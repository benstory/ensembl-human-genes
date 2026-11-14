#!/usr/bin/env Rscript
# =============================================================================
# run_hg19.R  -  MAIN SCRIPT: build the MANE Select pages on hg19 / GRCh37
#
# 1. Put this file next to mane_pages.R and mane_hg19.R.
# 2. Download the files listed in mane_hg19.R and set their paths below.
# 3. Run:  Rscript run_hg19.R        (or source("run_hg19.R") inside R)
#
# What it does:
#   MANE summary (ENST <-> NM_ <-> gene)  ┐
#   NCBI GRCh37 MANE mapping (GFF)        ┴─> hg19 MANE "GTF" ─┐
#   GRCh37 genome FASTA ───────────────────────────────────────┼─> one HTML page per
#   gnomAD text file (only MANE regions are read) ─────────────┘   MANE transcript + index.html
#   ...then checks every page's cDNA against MANE's own sequences.
# =============================================================================

## ---- CONFIG: edit these paths ---------------------------------------------
cfg <- list(
  mane_summary = "MANE.GRCh38.v1.5.summary.txt",
  # cDNA check; NULL to skip. Use the FULL ensembl_rna file: select_ensembl_rna has no
  # MANE Plus Clinical sequences, so those pages would go unchecked.
  mane_rna     = "MANE.GRCh38.v1.5.ensembl_rna.fna",
  ncbi_gff37   = "GCF_000001405.25_GRCh37.p13_genomic.gff",     # NCBI GRCh37 annotation (file ending _genomic.gff.gz)
  fasta        = "GRCh37.primary_assembly.genome.fa",    # uncompressed (or bgzipped)
  gnomad       = "~/Downloads/hg19_gnomad_genome.txt",   # gnomAD text file; NULL for no variants
  gnomad_cols  = c(chrom="#Chr", pos="Start", ref="Ref", alt="Alt", af="gnomAD_genome_ALL"),     # only if columns aren't recognised, e.g. c(chrom="CHR", pos="BP", ref="REF", alt="ALT", af="AF")
  # data source shown on the pages and in the footer. ANNOVAR's hg19_gnomad_genome table
  # (column gnomAD_genome_ALL) is gnomAD v2.0.1; v2.1.1 is the separate gnomad211_genome
  # table. If you switch files, change this label and gnomad_cols$af to match its header.
  gnomad_label = "gnomAD v2.0.1 genomes (GRCh37, via ANNOVAR)",
  genes        = NULL,     # e.g. c("BRCA1", "BRCA2") for a quick trial run; NULL = all MANE Select genes
  out_dir      = "mane_pages_hg19",
  intron_max   = 1000,
  variant_label = "gnomAD Genomes",
  min_af       = 1e-4,     # variants below 0.01% are left out
  af_default   = 0.001,    # AF filter the pages open with (0.1%)
  genome_build = "GRCh37 (hg19)",
  include_plus_clinical = TRUE,   # pages for MANE Plus Clinical transcripts too
  flags_gtf    = NULL     # optional: current Ensembl GRCh38 GTF for extra tags (basic, CCDS, APPRIS)
)


# optional: Rscript run_hg19.R my_config.R   (a file that changes entries of cfg)
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) for (f in args) source(f, local = TRUE)

## ---- load the functions ------------------------------------------------------
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else
    if (!is.null(sys.frames()[[1]]$ofile)) dirname(normalizePath(sys.frames()[[1]]$ofile)) else getwd()
})
source(file.path(here, "mane_pages.R"))
source(file.path(here, "mane_hg19.R"))

need <- c(cfg$mane_summary, cfg$fasta, cfg$ncbi_gff37, cfg$mane_rna, cfg$gnomad, cfg$flags_gtf)
miss <- need[!file.exists(need)]
if (length(miss)) stop("Missing file(s):\n  ", paste(miss, collapse = "\n  "), "\nEdit the CONFIG block in run_hg19.R.")
if (isTRUE(cfg$include_plus_clinical) && !is.null(cfg$mane_rna) && grepl("select_", basename(cfg$mane_rna), fixed = TRUE))
  warning("mane_rna looks like the Select-only file (", basename(cfg$mane_rna), "): MANE Plus Clinical pages ",
          "will not be checked. Use MANE.GRCh38.v1.x.ensembl_rna.fna instead.", call. = FALSE, immediate. = TRUE)
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
# remove the tables a previous run left behind: make_index.R reads whichever of
# them exist, so stale ones would put old QC / fallback flags on the new index
unlink(file.path(cfg$out_dir, c("pages.csv", "qc_cdna_vs_mane.csv", "mane_version_fallback.csv",
                                "mane_not_on_grch37.csv", "mane_vs_grch37_gtf.csv")))
t0 <- Sys.time()

## ---- 1. MANE annotation on hg19 ------------------------------------------------
message("\n[1/5] MANE summary + hg19 transcript models (NCBI GRCh37 mapping)")
summary <- read_mane_summary(cfg$mane_summary)
mane_status <- if (isTRUE(cfg$include_plus_clinical)) c("MANE Select", "MANE Plus Clinical") else "MANE Select"
mane_tags   <- gsub(" ", "_", mane_status)
gtf37 <- mane_gtf_from_ncbi_grch37(cfg$ncbi_gff37, summary, status = mane_status)
write.csv(attr(gtf37, "missing"), file.path(cfg$out_dir, "mane_not_on_grch37.csv"), row.names = FALSE)
write.csv(attr(gtf37, "version_fallback"), file.path(cfg$out_dir, "mane_version_fallback.csv"), row.names = FALSE)

if (!is.null(cfg$flags_gtf)) gtf37 <- add_transcript_flags(normalize_gtf(gtf37), cfg$flags_gtf)

## ---- 2. genome -------------------------------------------------------------------
message("\n[2/5] Genome FASTA")
if (requireNamespace("Rsamtools", quietly = TRUE)) {
  if (!file.exists(paste0(cfg$fasta, ".fai"))) {
    ok <- tryCatch({ Rsamtools::indexFa(cfg$fasta); TRUE }, error = function(e) FALSE)
    if (!ok) stop(cfg$fasta, " could not be indexed. If it is a plain .gz, unzip it first (gunzip) ",
                  "or recompress with bgzip, then run again.")
  }
  fasta <- Rsamtools::FaFile(cfg$fasta)        # reads only the bases each gene needs
} else {
  message("  Rsamtools not installed: loading the whole FASTA into memory (slow, ~3 GB)")
  fasta <- read_fasta_simple(cfg$fasta)
}

## ---- 3. gnomAD (only the regions the pages show) ---------------------------------
variants <- NULL
if (!is.null(cfg$gnomad)) {
  message("\n[3/5] gnomAD variants")
  # with cfg$genes set, read only those genes' regions (much faster for trial runs)
  reg_gtf <- if (is.null(cfg$genes)) gtf37 else
    gtf37[gtf37$gene_name %in% cfg$genes | sub("\\.[0-9]+$", "", gtf37$gene_id) %in% sub("\\.[0-9]+$", "", cfg$genes), ]
  if (!nrow(reg_gtf)) stop("none of cfg$genes (", paste(cfg$genes, collapse = ", "), ") is in the MANE annotation")
  variants <- read_variants_text(cfg$gnomad, regions = mane_regions(reg_gtf, tags = mane_tags, pad = cfg$intron_max + 100),
                                 cols = cfg$gnomad_cols)
} else message("\n[3/5] gnomAD: skipped")

## ---- 4. pages ---------------------------------------------------------------------
message("\n[4/5] Writing pages")
# provenance line shown in every page / RTF footer
ver <- function(x, rx) { m <- regmatches(basename(x), regexpr(rx, basename(x))); if (length(m)) m else basename(x) }
prov <- paste0("Generated ", format(Sys.Date()),
  " | MANE ", ver(cfg$mane_summary, "v[0-9][0-9.]*[0-9]"),
  " | NCBI GRCh37 MANE mapping ", ver(cfg$ncbi_gff37, "RS_[0-9_]*[0-9]"),
  " | genome ", basename(cfg$fasta),
  if (!is.null(cfg$gnomad)) paste0(" | variants ", cfg$gnomad_label, " (AF >= ", cfg$min_af * 100, "%)") else "",
  " | ", cfg$genome_build)
pages <- make_mane_pages(gtf37, fasta, variants = variants, out_dir = cfg$out_dir, genes = cfg$genes, tags = mane_tags,
                         intron_max = cfg$intron_max, variant_label = cfg$variant_label,
                         min_af = cfg$min_af, af_default = cfg$af_default,
                         genome_build = cfg$genome_build, provenance = prov,
                         variant_source = if (!is.null(cfg$gnomad)) paste0(cfg$gnomad_label, ", AF >= ", cfg$min_af * 100, "%"),
                         reference_cdna = cfg$mane_rna, verbose = FALSE)
write.csv(pages, file.path(cfg$out_dir, "pages.csv"), row.names = FALSE)

## ---- 5. QC: cDNA built from hg19 vs MANE's own transcript sequence -------------
if (!is.null(cfg$mane_rna)) {
  message("\n[5/5] cDNA check against MANE sequences")
  qc <- check_cdna_against_mane(pages, cfg$mane_rna)
  write.csv(qc, file.path(cfg$out_dir, "qc_cdna_vs_mane.csv"), row.names = FALSE)
} else message("\n[5/5] cDNA check: skipped")

# simple, fast index with risky pages flagged (replaces the table index)
source(file.path(here, "make_index.R"))
make_simple_index(cfg$out_dir)

message(sprintf("\nDone in %.1f min. Open %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
                file.path(normalizePath(cfg$out_dir), "index.html")))
