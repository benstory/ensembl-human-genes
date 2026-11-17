#!/usr/bin/env Rscript
# =============================================================================
# run_hg38.R  -  build MANE pages on hg38 / GRCh38
#
# MANE is native to GRCh38, so no NCBI remapping is needed: the MANE GTF is read
# directly and every page should match MANE exactly (QC is a sanity check here).
#
#   Rscript run_hg38.R            (or source("run_hg38.R") inside R)
#   Rscript run_hg38.R my_config.R   (optional: a file that changes entries of cfg)
#
# Downloads:
#   https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
#     MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz   transcript models (replaces the NCBI GRCh37 GFF)
#     MANE.GRCh38.v1.5.summary.txt.gz           NM_ accessions for file names and NCBI links
#     MANE.GRCh38.v1.5.ensembl_rna.fna.gz       QC (full file, not select_ensembl_rna)
#     MANE.GRCh38.v1.5.ensembl_protein.faa.gz   selenocysteine (U) positions
#   https://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/
#     Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz   -> gunzip, or re-bgzip
#     (UCSC hg38.fa works too; chr1 / 1 naming both work)
#   gnomAD (optional), e.g. ANNOVAR:
#     annotate_variation.pl -buildver hg38 -downdb -webfrom annovar gnomad41_genome humandb/
# =============================================================================

## ---- CONFIG: edit these ------------------------------------------------------
cfg <- list(
  mane_gtf     = "MANE.GRCh38.v1.5.ensembl_genomic.gtf.gz",          # <-- CHANGE: MANE GTF (hg38 models)
  mane_summary = "MANE.GRCh38.v1.5.summary.txt.gz",                  # <-- CHANGE: adds NM_ IDs (file names, NCBI links); NULL to skip
  mane_rna     = "MANE.GRCh38.v1.5.ensembl_rna.fna.gz",              # <-- CHANGE: QC; NULL to skip
  # selenocysteine positions (U) for selenoproteins such as SELENOP, used if the GTF has
  # no Selenocysteine rows; without it their Sec codons show as stops (*). NULL to skip
  mane_protein = "MANE.GRCh38.v1.5.ensembl_protein.faa.gz",          # <-- CHANGE
  ref_check_min = 0.95,   # stop if fewer variant REF alleles than this match the genome (wrong-build file); NULL = no check
  fasta        = "Homo_sapiens.GRCh38.dna.primary_assembly.fa",      # <-- CHANGE: hg38 genome (uncompressed or bgzipped)
  gnomad       = "humandb/hg38_gnomad41_genome.txt",                 # <-- CHANGE: hg38 variants; NULL for none
  # ANNOVAR's AF column name differs per table: check with  head -1 humandb/hg38_gnomad41_genome.txt
  gnomad_cols  = c(chrom = "#Chr", pos = "Start", ref = "Ref", alt = "Alt", af = "gnomad41_genome_AF"),
  gnomad_label = "gnomAD v4.1 genomes (GRCh38, via ANNOVAR)",
  genes        = NULL,                # <-- CHANGE: NULL = all MANE genes; trial run, e.g.
  # genes = c("NF1", "BRCA1", "BRCA2", "TP53", "CFTR", "DMD", "TTN", "KRAS", "SELENOP", "SHOX", "HBB"),
  out_dir      = "mane_pages_hg38",
  intron_max   = 1000,
  variant_label = "gnomAD Genomes",
  min_af       = 1e-4,
  af_default   = 0.001,
  genome_build = "GRCh38 (hg38)",     # drives UCSC hg38 / www.ensembl.org / gnomAD v4 links on the pages
  include_plus_clinical = TRUE,
  flags_gtf    = NULL                 # optional Ensembl GRCh38 GTF for basic / CCDS / APPRIS tags
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) for (f in args) source(f, local = TRUE)

## ---- load the functions ------------------------------------------------------
here <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(normalizePath(sub("^--file=", "", a[1]))) else
    if (!is.null(sys.frames()[[1]]$ofile)) dirname(normalizePath(sys.frames()[[1]]$ofile)) else getwd()
})
source(file.path(here, "mane_pages.R"))
source(file.path(here, "mane_hg19.R"))   # only for check_cdna_against_mane() / add_transcript_flags()

## Sec codon rows from the MANE protein FASTA: the U positions of each protein, placed on
## the first base (transcript orientation) of their codon in the CDS rows
sec_rows_from_protein <- function(gtf, faa, summary = NULL) {
  p <- read_fasta_simple(faa)
  names(p) <- sub("\\.[0-9]+$", "", sub("\\|.*$", "", names(p)))
  p <- p[grepl("U", p, fixed = TRUE)]
  if (!length(p)) return(NULL)
  cds <- gtf[gtf$feature == "CDS", , drop = FALSE]
  cds$pkey <- sub("\\.[0-9]+$", "", cds$protein_id)
  if (!is.null(summary) && any(is.na(cds$pkey)))      # GTF without protein_id: ENST -> ENSP via the summary
    cds$pkey[is.na(cds$pkey)] <- sub("\\.[0-9]+$", "", summary$Ensembl_prot[match(
      sub("\\.[0-9]+$", "", cds$transcript_id[is.na(cds$pkey)]), summary$enst)])
  cds <- cds[!is.na(cds$pkey) & cds$pkey %in% names(p), , drop = FALSE]
  if (!nrow(cds)) return(NULL)
  do.call(rbind, lapply(split(cds, cds$transcript_id), function(r) {
    minus <- r$strand[1] == "-"
    r <- r[order(r$start, decreasing = minus), ]
    pos <- unlist(Map(function(s, e) if (minus) e:s else s:e, r$start, r$end))
    k <- as.integer(gregexpr("U", p[[r$pkey[1]]], fixed = TRUE)[[1]])
    g <- pos[(k - 1L) * 3L + 1L]
    d <- r[rep(1L, length(g)), setdiff(names(r), "pkey")]
    d$feature <- "Selenocysteine"; d$start <- d$end <- g; d$frame <- NA; d
  }))
}

## variant REF alleles vs the genome: catches a variant file from the wrong build
check_variant_ref <- function(v, get_seq, min_match = 0.95, n = 2000) {
  v <- v[grepl("^[ACGT]+$", toupper(v$REF)), , drop = FALSE]      # skip ANNOVAR insertions ("-")
  if (nrow(v) < 20) return(invisible(NA))
  v <- v[unique(round(seq(1, nrow(v), length.out = min(n, nrow(v))))), , drop = FALSE]
  g <- mapply(function(c, p, r) toupper(tryCatch(get_seq(c, p, p + nchar(r) - 1L), error = function(e) "")),
              v$CHROM, as.integer(v$POS), v$REF, USE.NAMES = FALSE)
  ok <- mean(g == toupper(v$REF))
  message(sprintf("  REF matches the genome for %.1f%% of %s checked variants", 100 * ok, fmt_int(nrow(v))))
  if (ok < min_match) stop(sprintf(paste0("Only %.0f%% of variant REF alleles match %s (about 25%% = random). ",
    "The variant file is probably for another genome build: check cfg$gnomad."), 100 * ok, basename(cfg$fasta)), call. = FALSE)
  invisible(ok)
}

need <- c(cfg$mane_gtf, cfg$mane_summary, cfg$mane_protein, cfg$fasta, cfg$mane_rna, cfg$gnomad, cfg$flags_gtf)
miss <- need[!file.exists(need)]
if (length(miss)) stop("Missing file(s):\n  ", paste(miss, collapse = "\n  "), "\nEdit the CONFIG block in run_hg38.R.")
dir.create(cfg$out_dir, showWarnings = FALSE, recursive = TRUE)
unlink(file.path(cfg$out_dir, c("pages.csv", "qc_cdna_vs_mane.csv", "mane_version_fallback.csv",
                                "mane_not_on_grch37.csv", "mane_vs_grch37_gtf.csv")))
t0 <- Sys.time()

## ---- 1. MANE annotation (native GRCh38) -----------------------------------------
message("\n[1/5] MANE GRCh38 GTF")
mane_status <- if (isTRUE(cfg$include_plus_clinical)) c("MANE Select", "MANE Plus Clinical") else "MANE Select"
mane_tags   <- gsub(" ", "_", mane_status)
gtf38 <- read_gtf_simple(cfg$mane_gtf, keep = mane_tags)
# PAR genes (e.g. SHOX): keep the chrX copy only, so a chrY copy is never merged into the
# same model or written as a second page (GENCODE-style IDs carry a _PAR_Y suffix)
chr_key <- sub("^chr", "", gtf38$seqname)
tx_base <- sub("_PAR_Y$", "", gtf38$transcript_id)
gtf38 <- gtf38[!(chr_key == "Y" & tx_base %in% tx_base[chr_key == "X"]), ]
# the MANE GTF has no NM_ accession; the page builder takes it from transcript_name
summary <- NULL
if (!is.null(cfg$mane_summary)) {
  summary <- read_mane_summary(cfg$mane_summary)
  nm <- summary$RefSeq_nuc[match(sub("\\.[0-9]+$", "", gtf38$transcript_id), summary$enst)]
  gtf38$transcript_name <- ifelse(is.na(nm), gtf38$transcript_name, nm)
}
# selenocysteine: Sec rows may lack the MANE tag (so the read above drops them); read them
# separately, else take the U positions from the MANE protein FASTA
has_sec <- unique(gtf38$transcript_id[gtf38$feature == "Selenocysteine"])
sec <- tryCatch(read_gtf_simple(cfg$mane_gtf, keep = "Selenocysteine"), error = function(e) NULL)
if (!is.null(sec)) sec <- sec[sec$feature == "Selenocysteine" & sec$transcript_id %in% gtf38$transcript_id, , drop = FALSE]
if ((is.null(sec) || !nrow(sec)) && !is.null(cfg$mane_protein)) sec <- sec_rows_from_protein(gtf38, cfg$mane_protein, summary)
if (!is.null(sec) && nrow(sec <- sec[!(sec$transcript_id %in% has_sec), , drop = FALSE])) {
  gtf38 <- rbind(gtf38, sec[, names(gtf38)])
  message("  selenocysteine codons added for ", length(unique(sec$transcript_id)), " transcript(s)")
}
if (!is.null(cfg$genes)) {
  found <- cfg$genes[cfg$genes %in% gtf38$gene_name | cfg$genes %in% sub("\\.[0-9]+$", "", gtf38$gene_id)]
  if (length(lost <- setdiff(cfg$genes, found))) warning("not in the MANE GTF: ", paste(lost, collapse = ", "), call. = FALSE)
}
if (!is.null(cfg$flags_gtf)) gtf38 <- add_transcript_flags(gtf38, cfg$flags_gtf)

## ---- 2. genome -------------------------------------------------------------------
message("\n[2/5] Genome FASTA")
if (requireNamespace("Rsamtools", quietly = TRUE)) {
  if (!file.exists(paste0(cfg$fasta, ".fai"))) {
    ok <- tryCatch({ Rsamtools::indexFa(cfg$fasta); TRUE }, error = function(e) FALSE)
    if (!ok) stop(cfg$fasta, " could not be indexed. If it is a plain .gz, gunzip it or recompress with bgzip.")
  }
  fasta <- Rsamtools::FaFile(cfg$fasta)
} else {
  message("  Rsamtools not installed: loading the whole FASTA into memory (slow, ~3 GB)")
  fasta <- read_fasta_simple(cfg$fasta)
}

## ---- 3. gnomAD (only the regions the pages show) ---------------------------------
variants <- NULL
if (!is.null(cfg$gnomad)) {
  message("\n[3/5] gnomAD variants")
  reg_gtf <- if (is.null(cfg$genes)) gtf38 else
    gtf38[gtf38$gene_name %in% cfg$genes | sub("\\.[0-9]+$", "", gtf38$gene_id) %in% sub("\\.[0-9]+$", "", cfg$genes), ]
  if (!nrow(reg_gtf)) stop("none of cfg$genes (", paste(cfg$genes, collapse = ", "), ") is in the MANE annotation")
  variants <- read_variants_text(cfg$gnomad, regions = mane_regions(reg_gtf, tags = mane_tags, pad = cfg$intron_max + 100),
                                 cols = cfg$gnomad_cols)
  if (!is.null(cfg$ref_check_min)) check_variant_ref(variants, make_seq_getter(fasta), cfg$ref_check_min)
} else message("\n[3/5] gnomAD: skipped")

## ---- 4. pages ---------------------------------------------------------------------
message("\n[4/5] Writing pages")
ver <- function(x, rx) { m <- regmatches(basename(x), regexpr(rx, basename(x))); if (length(m)) m else basename(x) }
prov <- paste0("Generated ", format(Sys.Date()),
  " | MANE ", ver(cfg$mane_gtf, "v[0-9][0-9.]*[0-9]"), " GRCh38 GTF",
  " | genome ", basename(cfg$fasta),
  if (!is.null(cfg$gnomad)) paste0(" | variants ", cfg$gnomad_label, " (AF >= ", cfg$min_af * 100, "%)") else "",
  " | ", cfg$genome_build)
pages <- make_mane_pages(gtf38, fasta, variants = variants, out_dir = cfg$out_dir, genes = cfg$genes, tags = mane_tags,
                         intron_max = cfg$intron_max, variant_label = cfg$variant_label,
                         min_af = cfg$min_af, af_default = cfg$af_default,
                         genome_build = cfg$genome_build, provenance = prov,
                         variant_source = if (!is.null(cfg$gnomad)) paste0(cfg$gnomad_label, ", AF >= ", cfg$min_af * 100, "%"),
                         reference_cdna = cfg$mane_rna, verbose = FALSE)
write.csv(pages, file.path(cfg$out_dir, "pages.csv"), row.names = FALSE)
stops <- vapply(attr(pages, "models"), function(m) grepl("\\*.", m$protein), TRUE)
if (any(stops)) warning("internal stop (*) in the protein of: ", paste(names(stops)[stops], collapse = ", "),
                        " - a selenoprotein without Sec annotation? Set cfg$mane_protein.", call. = FALSE)

## ---- 5. QC: should be 100% identical on hg38 ----------------------------------------
if (!is.null(cfg$mane_rna)) {
  message("\n[5/5] cDNA check against MANE sequences")
  qc <- check_cdna_against_mane(pages, cfg$mane_rna)
  write.csv(qc, file.path(cfg$out_dir, "qc_cdna_vs_mane.csv"), row.names = FALSE)
  bad <- qc[!is.na(qc$identical) & !qc$identical, ]
  if (nrow(bad)) warning(nrow(bad), " hg38 page(s) differ from MANE - unexpected on GRCh38, check the FASTA: ",
                         paste(bad$gene, bad$transcript, collapse = ", "), call. = FALSE)
} else message("\n[5/5] cDNA check: skipped")

source(file.path(here, "make_index.R"))
make_simple_index(cfg$out_dir, title = "MANE transcripts (hg38)")

message(sprintf("\nDone in %.1f min. Open %s", as.numeric(difftime(Sys.time(), t0, units = "mins")),
                file.path(normalizePath(cfg$out_dir), "index.html")))
