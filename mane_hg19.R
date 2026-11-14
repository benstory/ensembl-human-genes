# =============================================================================
# mane_hg19.R : MANE Select transcripts on hg19 / GRCh37 for mane_pages.R
#
# MANE is defined on GRCh38 only. NCBI's GRCh37.p13 annotation release maps
# the MANE RefSeq transcripts (NM_) onto GRCh37. Because each MANE NM_ and its
# ENST are identical (same exons, UTRs and CDS), the NM_ mapping is the MANE
# transcript's structure on GRCh37. We relabel it with the ENST / ENSG from the
# MANE summary file:
#     mane_gtf <- mane_gtf_from_ncbi_grch37(gff, summary)
# then check the built cDNAs against MANE's own transcript sequences:
#     check_cdna_against_mane(pages, "MANE.GRCh38.v1.5.ensembl_rna.fna.gz")
#
# Downloads (MANE v1.5 is current at the time of writing):
#   MANE (NCBI): https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/
#     MANE.GRCh38.v1.5.summary.txt.gz        ENST <-> NM_ <-> gene, status
#     MANE.GRCh38.v1.5.ensembl_rna.fna.gz    MANE transcript sequences (QC). Use this
#                                            full file, not select_ensembl_rna: the
#                                            Select-only file has no MANE Plus Clinical.
#   NCBI GRCh37.p13 annotation with MANE mapped:
#     https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/annotation_releases/GCF_000001405.25-RS_2024_09/
#     -> the file ending in _genomic.gff.gz
#   GRCh37 genome FASTA (Ensembl):
#     https://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
#
# Needs only base R (plus mane_pages.R for the page step).
# =============================================================================

.open_text <- function(path) if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")

#' Read the MANE summary file. Adds `enst` / `nm` columns without versions.
read_mane_summary <- function(path) {
  s <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, quote = "", comment.char = "")
  names(s) <- sub("^#", "", names(s))
  need <- c("symbol", "Ensembl_Gene", "Ensembl_nuc", "RefSeq_nuc", "MANE_status")
  miss <- setdiff(need, names(s)); if (length(miss)) stop("summary file lacks: ", paste(miss, collapse = ", "))
  s$enst <- sub("\\.[0-9]+$", "", s$Ensembl_nuc)
  s$nm   <- sub("\\.[0-9]+$", "", s$RefSeq_nuc)
  s
}

# GFF3 key=value. Not URL-decoded: the values used here (accessions,
# chromosome names) never contain "%" - and URLdecode() is slow and turns NA
# into the string "NA".
.attr1 <- function(a, key) {
  m <- regexpr(paste0("(^|;)", key, "=[^;]*"), a)
  out <- rep(NA_character_, length(a))
  out[m > 0] <- sub(paste0("^;?", key, "="), "", regmatches(a, m))
  out
}
.gtf_attr1 <- function(a, key) {                  # GTF key "value"
  m <- regexpr(paste0(key, ' "[^"]*"'), a)
  out <- rep(NA_character_, length(a))
  out[m > 0] <- sub(paste0(key, ' "([^"]*)"'), "\\1", regmatches(a, m))
  out
}

#' Build a GTF-like data.frame of MANE transcripts on GRCh37 from NCBI's
#' GRCh37.p13 annotation GFF3, labelled with Ensembl IDs. Pass the result
#' straight to make_mane_pages().
#'
#' @param gff      path to the NCBI GRCh37.p13 *_genomic.gff.gz
#' @param summary  data.frame from read_mane_summary()
#' @param status   which MANE sets to keep
#' @param chunk    lines read at a time (memory control)
mane_gtf_from_ncbi_grch37 <- function(gff, summary, status = c("MANE Select", "MANE Plus Clinical"),
                                      chunk = 500000L) {
  want <- summary[summary$MANE_status %in% status, ]
  if (!nrow(want)) stop("no rows with MANE_status in: ", paste(status, collapse = ", "))
  nm_set <- want$nm                                          # NM_/NR_ without version (fallback-aware)
  con <- .open_text(gff); on.exit(close(con))
  regions <- character(0); keep <- character(0)
  repeat {
    l <- readLines(con, n = chunk, warn = FALSE)
    if (!length(l)) break
    l <- l[!startsWith(l, "#")]
    regions <- c(regions, l[grepl("\tregion\t", l, fixed = TRUE)])
    # cheap text filters first (exon/CDS rows of NM_/NR_ transcripts), then
    # pull the accession only from the survivors
    cand <- l[grepl("\t(exon|CDS)\t", l)]
    cand <- cand[grepl("Parent=rna-N[MR]_", cand)]
    acc <- sub("-[0-9]+$", "", sub("^rna-", "", .attr1(cand, "Parent")))
    keep <- c(keep, cand[sub("\\.[0-9]+$", "", acc) %in% nm_set])
  }
  if (!length(keep)) stop("no exon/CDS rows for the MANE accessions were found in ", gff,
                          " (check the file is an NCBI GRCh37.p13 annotation GFF)")
  # NC_ accession -> chromosome name ("13", "X", "MT") from the region rows
  rf <- do.call(rbind, strsplit(regions, "\t", fixed = TRUE))
  chr_of <- stats::setNames(.attr1(rf[, 9], "chromosome"), rf[, 1])
  mt <- grepl("genome=mitochondrion", rf[, 9]); chr_of[mt] <- "MT"
  chr_of <- chr_of[!is.na(chr_of)]

  f <- do.call(rbind, strsplit(keep, "\t", fixed = TRUE))
  df <- data.frame(seqid = f[, 1], feature = f[, 3], start = as.integer(f[, 4]), end = as.integer(f[, 5]),
                   strand = f[, 7], frame = suppressWarnings(as.integer(f[, 8])),
                   parent = sub("^rna-", "", .attr1(f[, 9], "Parent")), stringsAsFactors = FALSE)
  df$acc <- sub("-[0-9]+$", "", df$parent)
  # Selenocysteine: NCBI records Sec codons as transl_except on the CDS rows,
  # e.g. transl_except=(pos:complement(123..125)%2Caa:Sec). Turn them into
  # Selenocysteine rows so those codons translate to U instead of a stop.
  # A codon split by an intron is written pos:complement(join(1..2%2C10..10)); the row
  # then spans the whole range, and the page only uses its 5'-most base (start on +,
  # end on -), which is the codon's first base either way.
  sec_i <- which(df$feature == "CDS" & grepl("aa:Sec", f[, 9], fixed = TRUE))
  if (length(sec_i)) {
    te <- gsub("%2C", ",", f[sec_i, 9], fixed = TRUE)
    hits <- regmatches(te, gregexpr("pos:[cjoinmplet0-9.,()]+aa:Sec", te))
    sec <- do.call(rbind, lapply(seq_along(sec_i), function(j) {
      if (!length(hits[[j]])) return(NULL)
      xy <- do.call(rbind, lapply(hits[[j]], function(h) {
        v <- as.integer(regmatches(h, gregexpr("[0-9]+", h))[[1]]); c(min(v), max(v)) }))
      d <- df[rep(sec_i[j], nrow(xy)), ]; d$feature <- "Selenocysteine"; d$start <- xy[, 1]; d$end <- xy[, 2]; d$frame <- NA; d
    }))
    if (!is.null(sec)) df <- unique(rbind(df, sec))
  }
  # primary chromosomes only (NC_), one alignment per accession (avoid chrY PAR copy)
  df <- df[startsWith(df$seqid, "NC_") & df$seqid %in% names(chr_of), ]
  df$chrom <- unname(chr_of[df$seqid])
  # Version matching: use the exact MANE version when NCBI placed it on GRCh37;
  # otherwise fall back to the version that is placed (the highest, if several).
  # Pages built from a fallback carry a version note, and the cDNA check against
  # MANE's own sequence (reference_cdna) shows whether anything differs.
  df$base <- sub("\\.[0-9]+$", "", df$acc)
  avail <- lapply(split(df$acc, df$base), unique)
  vnum <- function(a) suppressWarnings(as.integer(sub("^.*\\.", "", a)))
  chosen <- vapply(seq_len(nrow(want)), function(k) {
    a <- avail[[want$nm[k]]]
    if (is.null(a)) NA_character_ else if (want$RefSeq_nuc[k] %in% a) want$RefSeq_nuc[k] else a[order(-vnum(a))][1]
  }, "")
  df <- df[df$acc %in% chosen[!is.na(chosen)], ]
  pick <- tapply(seq_len(nrow(df)), df$acc, function(i) {
    al <- unique(df$parent[i]); al_chr <- vapply(al, function(a) df$chrom[i][df$parent[i] == a][1], "")
    al[order(al_chr == "Y", al)][1]
  })
  df <- df[df$parent == pick[df$acc], ]

  m <- want[match(df$base, want$nm), ]
  # A fallback page is the RefSeq version that IS on GRCh37 (e.g. NM_000267.3), a real
  # transcript - not MANE. It is labelled with its own accession; the MANE IDs
  # (NM_000267.4 / ENST00000356175.8) are kept only for the warning and the sequence
  # comparison. Ensembl canonical / GENCODE Primary belong to the MANE ENST, so they are not set.
  isfb <- df$acc != m$RefSeq_nuc
  out <- data.frame(
    seqname = df$chrom, source = "NCBI_GRCh37_MANE", feature = df$feature,
    start = df$start, end = df$end, score = ".", strand = df$strand, frame = df$frame,
    gene_id = m$Ensembl_Gene, gene_name = m$symbol, gene_biotype = NA_character_,
    transcript_id = ifelse(isfb, df$acc, m$Ensembl_nuc), transcript_name = df$acc,
    transcript_biotype = ifelse(grepl("^NR_", m$RefSeq_nuc), "non_coding", "protein_coding"),
    # By definition (GENCODE): the Ensembl canonical transcript is the MANE Select
    # one where it exists, and GENCODE Primary contains all MANE Select and MANE
    # Plus Clinical transcripts. Other tags can be added with add_transcript_flags().
    tag = ifelse(isfb, paste0(gsub(" ", "_", m$MANE_status), ",MANE_version_fallback"),
          ifelse(m$MANE_status == "MANE Select", "MANE_Select,Ensembl_canonical,GENCODE_Primary",
                 paste0(gsub(" ", "_", m$MANE_status), ",GENCODE_Primary"))),
    refseq_id = df$acc,
    mane_refseq = ifelse(isfb, m$RefSeq_nuc, NA_character_),     # set only on fallback pages
    mane_enst   = ifelse(isfb, m$Ensembl_nuc, NA_character_),
    exon_id = NA_character_, stringsAsFactors = FALSE)
  # exon numbers in transcript order (one vectorised pass)
  ex <- which(out$feature == "exon")
  o <- ex[order(out$transcript_id[ex], ifelse(out$strand[ex] == "-", -out$start[ex], out$start[ex]))]
  out$exon_number <- NA_character_
  out$exon_number[o] <- as.character(stats::ave(seq_along(o), out$transcript_id[o], FUN = seq_along))
  placed <- unique(m$RefSeq_nuc)
  attr(out, "missing") <- want[!want$RefSeq_nuc %in% placed, c("symbol", "Ensembl_nuc", "RefSeq_nuc", "MANE_status")]
  fb <- !is.na(chosen) & chosen != want$RefSeq_nuc & want$RefSeq_nuc %in% placed
  attr(out, "version_fallback") <- data.frame(symbol = want$symbol[fb], MANE_status = want$MANE_status[fb],
    Ensembl_nuc = want$Ensembl_nuc[fb], MANE_RefSeq = want$RefSeq_nuc[fb], GRCh37_RefSeq = chosen[fb],
    stringsAsFactors = FALSE)
  message(sprintf("MANE on GRCh37 (NCBI mapping): %d transcripts placed (%d using an older/other RefSeq version), %d not found",
                  length(placed), sum(fb), nrow(attr(out, "missing"))))
  out
}

# ---- optional: extra transcript tags from a GTF ----------------------------------

#' Read only the rows of the given transcripts (unversioned ENST) from a GTF.
read_gtf_transcripts <- function(path, enst, chunk = 500000L) {
  con <- .open_text(path); on.exit(close(con))
  out <- character(0)
  repeat {
    l <- readLines(con, n = chunk, warn = FALSE)
    if (!length(l)) break
    l <- l[!startsWith(l, "#")]
    id <- sub("\\.[0-9]+$", "", .gtf_attr1(l, "transcript_id"))
    out <- c(out, l[!is.na(id) & id %in% enst])
  }
  if (!length(out)) return(data.frame(seqname = character(0), source = character(0), feature = character(0),
    start = integer(0), end = integer(0), score = character(0), strand = character(0), frame = character(0),
    attribute = character(0), stringsAsFactors = FALSE))
  f <- do.call(rbind, strsplit(out, "\t", fixed = TRUE))
  data.frame(seqname = f[, 1], source = f[, 2], feature = f[, 3], start = as.integer(f[, 4]),
             end = as.integer(f[, 5]), score = f[, 6], strand = f[, 7], frame = f[, 8],
             attribute = f[, 9], stringsAsFactors = FALSE)
}

#' Add tags (e.g. basic, CCDS, appris_principal_1) from a GTF that has them -
#' typically the current Ensembl GRCh38 GTF; tags are per transcript, so the
#' assembly does not matter. Only tags of the SAME ENST version are used (versions
#' are compared whether they are part of transcript_id or a transcript_version column).
add_transcript_flags <- function(gtf, flags_gtf, keep = c("basic", "CCDS", "Ensembl_canonical", "GENCODE_Primary",
                                                         "appris_principal_1", "appris_principal_2", "appris_principal_3",
                                                         "appris_principal_4", "appris_principal_5")) {
  ids <- unique(sub("\\.[0-9]+$", "", gtf$transcript_id))
  f <- read_gtf_transcripts(flags_gtf, ids)
  f <- f[f$feature == "transcript", ]
  if (!nrow(f)) { message("add_transcript_flags: no matching transcripts in ", flags_gtf); return(gtf) }
  tid <- .gtf_attr1(f$attribute, "transcript_id"); tv <- .gtf_attr1(f$attribute, "transcript_version")
  vid <- ifelse(is.na(tv) | grepl("\\.[0-9]+$", tid), tid, paste0(tid, ".", tv))
  tg <- vapply(regmatches(f$attribute, gregexpr('tag "[^"]*"', f$attribute)),
               function(v) paste(intersect(sub('^tag "(.*)"$', "\\1", v), keep), collapse = ","), "")
  ccds <- .gtf_attr1(f$attribute, "ccds_id")
  add <- stats::setNames(tg, vid); cc <- stats::setNames(ccds, vid)
  tv_g <- if ("transcript_version" %in% names(gtf)) gtf$transcript_version else NA_character_
  vid_g <- ifelse(is.na(tv_g) | grepl("\\.[0-9]+$", gtf$transcript_id), gtf$transcript_id,
                  paste0(gtf$transcript_id, ".", tv_g))
  hit <- !is.na(vid_g) & vid_g %in% names(add)
  extra <- add[vid_g[hit]]
  # merge once per distinct (existing tags, new tags) pair - there are only a handful -
  # instead of once per row
  old <- gtf$tag[hit]; old[is.na(old)] <- ""
  pair <- paste(old, unname(extra), sep = "\r"); up <- unique(pair)
  merged <- vapply(strsplit(up, "\r", fixed = TRUE), function(p)
    paste(unique(unlist(strsplit(p, ",", fixed = TRUE))), collapse = ","), "")
  gtf$tag[hit] <- merged[match(pair, up)]
  if (!"ccds_id" %in% names(gtf)) gtf$ccds_id <- NA_character_
  gtf$ccds_id[hit] <- ifelse(is.na(gtf$ccds_id[hit]), unname(cc[vid_g[hit]]), gtf$ccds_id[hit])
  message(sprintf("add_transcript_flags: tags added for %d of %d transcripts (same ENST version)",
                  length(unique(gtf$transcript_id[hit])), length(unique(gtf$transcript_id))))
  gtf
}

# ---- QC: built cDNA vs MANE transcript sequences -------------------------------

#' Compare each page's cDNA (built from the GRCh37 genome) with MANE's own
#' transcript sequence. Differences mean the GRCh37 reference differs from
#' the MANE transcript at those bases (or the model is not the MANE one).
#'
#' @param pages  return value of make_mane_pages() (it carries the models)
#' @param rna    path to MANE.GRCh38.v1.5.ensembl_rna.fna.gz. Must be the Ensembl
#'               (ENST-keyed) file: pages are matched by ENST, so a refseq_rna file
#'               (NM_-keyed, may include poly-A tails) would match nothing.
check_cdna_against_mane <- function(pages, rna) {
  models <- attr(pages, "models")
  ids <- sub("\\.[0-9]+$", "", vapply(models, function(m)
    if (!is.null(m$mane_enst) && !is.na(m$mane_enst)) m$mane_enst else m$transcript_id, ""))
  con <- .open_text(rna); on.exit(close(con))
  l <- readLines(con, warn = FALSE)
  h <- which(startsWith(l, ">"))
  nm <- sub("\\.[0-9]+$", "", sub("^>(\\S+).*", "\\1", l[h]))
  # headers may be "ENST...|..." style; take the first ENST/NM-like token
  nm <- sub("^([A-Z]+_?[0-9]+).*", "\\1", nm)
  keep <- which(nm %in% ids)
  ends <- c(h[-1] - 1, length(l))
  ref <- stats::setNames(vapply(keep, function(k) toupper(paste(l[(h[k] + 1):ends[k]], collapse = "")), ""), nm[keep])
  res <- data.frame(gene = vapply(models, `[[`, "", "gene_name"),
                    transcript = vapply(models, `[[`, "", "transcript_id"),
                    compared_with = vapply(models, function(m) if (!is.null(m$mane_refseq) && !is.na(m$mane_refseq))
                      paste0("MANE ", m$mane_refseq, " / ", m$mane_enst, " (page is an older/other version)") else "MANE", ""),
                    built_bp = vapply(models, function(m) nchar(m$cdna), 0L),
                    mane_bp = NA_integer_, mismatches = NA_integer_, identical = NA,
                    stringsAsFactors = FALSE, row.names = NULL)
  for (k in seq_along(models)) {
    r <- ref[ids[k]]; if (is.na(r)) next
    a <- strsplit(models[[k]]$cdna, "")[[1]]; b <- strsplit(r, "")[[1]]
    res$mane_bp[k] <- length(b)
    res$mismatches[k] <- if (length(a) == length(b)) sum(a != b) else NA
    res$identical[k] <- identical(models[[k]]$cdna, unname(r))
  }
  miss <- is.na(res$mane_bp)
  if (any(miss)) {
    pc <- vapply(models, function(m) "MANE_Plus_Clinical" %in% m$tags, TRUE)
    warning(sprintf("%d page(s) not found in %s and therefore NOT checked%s", sum(miss), basename(rna),
      if (any(miss & pc)) sprintf(" (%d MANE Plus Clinical - use the full MANE.GRCh38.v1.x.ensembl_rna.fna.gz, not a Select-only file)",
                                  sum(miss & pc)) else ""), call. = FALSE)
  }
  message(sprintf("cDNA identical to MANE: %d / %d  (length differs: %d, base mismatches: %d)",
                  sum(res$identical, na.rm = TRUE), sum(!is.na(res$identical)),
                  sum(res$built_bp != res$mane_bp, na.rm = TRUE), sum(res$mismatches > 0, na.rm = TRUE)))
  res
}
