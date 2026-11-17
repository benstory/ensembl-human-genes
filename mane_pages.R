# =============================================================================
# mane_pages.R
#
# Build one self-contained HTML page per MANE Select transcript, named after
# the gene symbol (e.g. BRCA2.html). Each page has two views modelled on the
# Ensembl transcript pages:
#
#   * cDNA sequence : cDNA / CDS / translation lines, alternating exon colours,
#                     UTR highlight, alternating codons
#   * Exons         : exon / intron table with coordinates, phases, lengths and
#                     sequence (flanks, UTR, CDS, intron ends)
#
# plus an index.html listing every page.
#
# Expected inputs (already loaded in your R session):
#
#   gtf   : any of
#             - data.frame from rtracklayer::import("x.gtf") |> as.data.frame()
#             - GRanges from rtracklayer::import("x.gtf")
#             - data.frame with the 9 raw GTF columns (attribute column unparsed)
#             - read_gtf_simple("x.gtf[.gz]") below (base R, no packages)
#   fasta : any of
#             - Biostrings::DNAStringSet (readDNAStringSet("genome.fa"))
#             - Rsamtools::FaFile("genome.fa")  (indexed; best for whole genomes)
#             - named character vector (read_fasta_simple("x.fa") below)
#           Record names may be "13", "chr13", "13 dna:chromosome ...", or a
#           region extract such as "13:32314508-32401268" (samtools faidx
#           style); region records are mapped back onto chromosome coordinates.
#   variants (optional) : gnomAD / VCF data with allele frequencies, any of
#             - data.frame with CHROM, POS, REF, ALT, AF (ID optional; column
#               names are case-insensitive; or an INFO column holding AF=)
#             - VariantAnnotation::readVcf() object (uses info(vcf)$AF)
#             - vcfR::read.vcfR() object
#
# Usage:
#   source("mane_pages.R")
#   make_mane_pages(gtf, fasta, variants = gnomad, out_dir = "mane_pages")
#
# Page features: cDNA + Exons views, click / Shift+click / drag to select
# bases with a copy box (in a collapsible bottom drawer), variants
# coloured by AF, and a dropdown for how much intron sequence is shown.
#
# Only base R is required. Biostrings / Rsamtools / GenomicRanges are used
# automatically if your objects come from them.
# =============================================================================

# ---- small helpers ----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub("\"", "&quot;", x, fixed = TRUE)
}

fmt_int <- function(x) formatC(x, format = "d", big.mark = ",")

revcomp <- function(x) {
  if (!nzchar(x)) return(x)
  x <- chartr("ACGTRYKMBVDHNacgtrykmbvdhn", "TGCAYRMKVBHDNtgcayrmkvbhdn", x)
  intToUtf8(rev(utf8ToInt(x)))
}

.CODONS <- local({
  b <- c("T", "C", "A", "G")
  cod <- paste0(rep(b, each = 16), rep(rep(b, each = 4), 4), rep(b, 16))
  aa  <- strsplit("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "")[[1]]
  stats::setNames(aa, cod)
})

translate_dna <- function(s) {
  n <- nchar(s) %/% 3
  if (n == 0) return("")
  cod <- substring(s, seq(1, by = 3, length.out = n), seq(3, by = 3, length.out = n))
  aa <- .CODONS[toupper(cod)]
  aa[is.na(aa)] <- "X"
  paste(aa, collapse = "")
}

json_str <- function(x) {
  x <- gsub("\\\\", "\\\\\\\\", x)
  x <- gsub("\"", "\\\\\"", x)
  x <- gsub("</", "<\\/", x, fixed = TRUE)          # JSON "<\/" (keeps "</script" out of the page)
  # raw control characters are invalid inside JSON strings (JSON.parse would throw)
  x <- gsub("\n", "\\n", x, fixed = TRUE); x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE); x <- gsub("[[:cntrl:]]", " ", x)
  paste0("\"", x, "\"")
}

safe_filename <- function(x) gsub("[^A-Za-z0-9._-]", "_", x)

# ---- optional base-R readers (handy for small/test files) -------------------

#' Read a GTF without Bioconductor. `keep` holds fixed strings; a line is kept
#' if it contains any of them (default MANE Select + Plus Clinical), which keeps memory small on a
#' full genome GTF. Use keep = NULL to read everything. Handles .gz.
read_gtf_simple <- function(path, keep = c("MANE_Select", "MANE_Plus_Clinical"), chunk = 500000L) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  out <- list()
  repeat {
    l <- readLines(con, n = chunk, warn = FALSE)
    if (!length(l)) break
    l <- l[!startsWith(l, "#")]
    if (!is.null(keep)) l <- l[Reduce(`|`, lapply(keep, function(k) grepl(k, l, fixed = TRUE)))]
    if (length(l)) out[[length(out) + 1]] <- l
  }
  l <- unlist(out, use.names = FALSE)
  if (!length(l)) stop("No GTF lines kept from ", path)
  f <- do.call(rbind, strsplit(l, "\t", fixed = TRUE))
  df <- data.frame(seqname = f[, 1], source = f[, 2], feature = f[, 3],
                   start = as.integer(f[, 4]), end = as.integer(f[, 5]),
                   score = f[, 6], strand = f[, 7], frame = f[, 8],
                   attribute = f[, 9], stringsAsFactors = FALSE)
  normalize_gtf(df)
}

#' Read a (small) FASTA into a named character vector. Names are the first
#' word of each header line.
read_fasta_simple <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  l <- readLines(con, warn = FALSE)
  h <- which(startsWith(l, ">"))
  nm <- sub("\\s.*$", "", sub("^>", "", l[h]))
  ends <- c(h[-1] - 1, length(l))
  seqs <- mapply(function(a, b) if (b > a) paste(l[(a + 1):b], collapse = "") else "",
                 h, ends, USE.NAMES = FALSE)
  stats::setNames(seqs, nm)
}

# ---- GTF normalisation -------------------------------------------------------

.GTF_ATTRS <- c("gene_id", "gene_version", "gene_name", "gene_biotype", "gene_type",
                "transcript_id", "transcript_version", "transcript_name",
                "transcript_biotype", "transcript_type", "exon_number", "exon_id",
                "exon_version", "protein_id", "protein_version", "ccds_id", "tag")

parse_gtf_attributes <- function(attr, keys = .GTF_ATTRS) {
  res <- list()
  for (k in keys) {
    if (k == "tag") {
      m <- regmatches(attr, gregexpr('tag "[^"]*"', attr))
      res[[k]] <- vapply(m, function(v) paste(sub('^tag "(.*)"$', "\\1", v), collapse = ","), "")
    } else {
      pat <- paste0("(^|;)\\s*", k, ' "([^"]*)"')
      m <- regexpr(pat, attr, perl = TRUE)
      v <- rep(NA_character_, length(attr))
      v[m > 0] <- sub(paste0('.*', k, ' "([^"]*)"$'), "\\1", regmatches(attr, m))
      res[[k]] <- v
    }
  }
  as.data.frame(res, stringsAsFactors = FALSE)
}

#' Coerce the many shapes a GTF can take in R into one data.frame with
#' columns seqname, feature, start, end, strand, frame + attribute columns.
normalize_gtf <- function(gtf) {
  df <- if (methods::is(gtf, "GRanges") || methods::is(gtf, "GenomicRanges")) {
    as.data.frame(gtf)
  } else as.data.frame(gtf, stringsAsFactors = FALSE)

  if (ncol(df) == 9 && all(names(df) == paste0("V", 1:9))) {
    names(df) <- c("seqname", "source", "feature", "start", "end", "score",
                   "strand", "frame", "attribute")
  }
  ren <- c(seqnames = "seqname", seqid = "seqname", chr = "seqname", chrom = "seqname",
           type = "feature", phase = "frame", attributes = "attribute")
  for (o in names(ren)) if (o %in% names(df) && !(ren[[o]] %in% names(df))) {
    names(df)[names(df) == o] <- ren[[o]]
  }
  need <- c("seqname", "feature", "start", "end", "strand")
  miss <- setdiff(need, names(df))
  if (length(miss)) stop("GTF is missing column(s): ", paste(miss, collapse = ", "))

  # unparsed 9th column -> parse only what we need
  if ("attribute" %in% names(df) && !("transcript_id" %in% names(df))) {
    df <- cbind(df, parse_gtf_attributes(df$attribute))
  }
  # list columns (e.g. CharacterList `tag` from rtracklayer GFF3) -> "a,b,c"
  for (k in names(df)) if (is.list(df[[k]]) || methods::is(df[[k]], "List")) {
    df[[k]] <- vapply(as.list(df[[k]]), function(v) paste(as.character(v), collapse = ","), "")
  }
  for (k in setdiff(.GTF_ATTRS, names(df))) df[[k]] <- NA_character_
  if (!("frame" %in% names(df))) df$frame <- NA
  df$seqname <- as.character(df$seqname)
  df$feature <- as.character(df$feature)
  df$strand  <- as.character(df$strand)
  df$start   <- as.integer(df$start)
  df$end     <- as.integer(df$end)
  df$frame   <- suppressWarnings(as.integer(as.character(df$frame)))
  df$tag     <- as.character(df$tag)
  df
}

# ---- FASTA access ------------------------------------------------------------

#' Returns function(chr, start, end) -> uppercase sequence of exactly
#' end-start+1 characters (positions outside the available sequence are "N").
make_seq_getter <- function(fasta) {
  kind <- if (methods::is(fasta, "FaFile")) "fafile"
          else if (methods::is(fasta, "XStringSet")) "xss"
          else if (is.list(fasta)) "chr"
          else if (is.character(fasta)) "chr"
          else stop("Unsupported fasta object of class ", class(fasta)[1])
  if (kind == "fafile") {
    idx  <- Rsamtools::scanFaIndex(fasta)
    raw  <- as.character(GenomicRanges::seqnames(idx))
    lens <- as.integer(GenomicRanges::end(idx))
  } else if (kind == "xss") {
    raw  <- names(fasta); lens <- as.integer(Biostrings::width(fasta))
  } else {
    fasta <- unlist(fasta)
    raw  <- names(fasta); lens <- nchar(fasta)
  }
  if (is.null(raw) || any(!nzchar(raw))) stop("FASTA records need names")
  first <- sub("\\s.*$", "", raw)
  rx <- "^(.+):([0-9,]+)-([0-9,]+)$"            # samtools-style region record
  is_reg <- grepl(rx, first)
  contig <- ifelse(is_reg, sub(rx, "\\1", first), first)
  offset <- rep(0, length(first))
  offset[is_reg] <- as.numeric(gsub(",", "", sub(rx, "\\2", first[is_reg]))) - 1
  recs <- data.frame(i = seq_along(first), raw = raw, contig = contig,
                     offset = offset, len = lens, stringsAsFactors = FALSE)
  norm <- function(x) { x <- sub("^chr", "", x, ignore.case = TRUE); ifelse(x %in% c("M", "MT"), "MT", x) }
  recs$key <- norm(recs$contig)

  fetch <- function(r, s, e) {             # local 1-based coordinates
    if (kind == "xss") return(as.character(Biostrings::subseq(fasta[[r$i]], s, e)))
    if (kind == "fafile") {
      gr <- GenomicRanges::GRanges(r$raw, IRanges::IRanges(s, e))
      return(as.character(Rsamtools::scanFa(fasta, gr))[[1]])
    }
    substr(fasta[[r$i]], s, e)
  }

  function(chr, start, end) {
    cand <- recs[recs$key == norm(chr), , drop = FALSE]
    if (!nrow(cand)) stop("Sequence '", chr, "' not found in FASTA (have: ",
                          paste(utils::head(unique(recs$contig), 10), collapse = ", "), ")")
    out <- rep("N", end - start + 1)
    for (k in seq_len(nrow(cand))) {
      r <- cand[k, ]
      s <- max(start, r$offset + 1); e <- min(end, r$offset + r$len)
      if (s > e) next
      seg <- strsplit(toupper(fetch(r, s - r$offset, e - r$offset)), "")[[1]]
      out[(s - start + 1):(e - start + 1)] <- seg
    }
    paste(out, collapse = "")
  }
}

# ---- gnomAD / VCF variants ---------------------------------------------------

#' Allele-frequency bands and their colours. The seven colours are the ones
#' Ensembl uses for variant highlighting (with Ensembl's black/white text);
#' here they are assigned to frequency bands, most common first. Override with
#' make_mane_pages(..., af_bins = my_bins) using the same columns.
.AF_BINS <- data.frame(
  min   = c(0.05, 0.01, 0.001, 1e-4, 1e-5, 1e-6, 0),
  bg    = c("#ff0000", "#ff7f50", "#ffd700", "#76ee00", "#7ac5cd", "#a2b5cd", "#02599c"),
  fg    = c("#ffffff", "#000000", "#000000", "#000000", "#000000", "#000000", "#ffffff"),
  label = c("&ge;5%", "1&ndash;5%", "0.1&ndash;1%", "0.01&ndash;0.1%", "0.001&ndash;0.01%",
            "0.0001&ndash;0.001%", "&lt;0.0001%"),
  stringsAsFactors = FALSE)

# per-band highlight + legend CSS
af_bin_css <- function(bins) {
  i <- seq_len(nrow(bins)) - 1
  paste(c(sprintf("::highlight(v%d){background-color:%s;color:%s;text-decoration:underline dotted 2px %s}",
                  i, bins$bg, bins$fg, ifelse(bins$fg == "#ffffff", "rgba(255,255,255,.9)", "rgba(0,0,0,.75)")),
          sprintf(".vsw.vb%d{background:%s;color:%s}", i, bins$bg, bins$fg)), collapse = "\n")
}

#' Coerce variant data into data.frame(chrom, pos, ref, alt, af, id).
#' Accepts a data.frame (CHROM/POS/REF/ALT/AF[/ID], any case; or an INFO
#' column containing AF=), a VariantAnnotation VCF object, or a vcfR object.
#' Multi-allelic rows ("A,G" with AF "0.1,0.02") are split.
normalize_variants <- function(v) {
  if (is.null(v)) return(NULL)
  if (methods::is(v, "VCF")) {                       # VariantAnnotation
    rr <- SummarizedExperiment::rowRanges(v)
    alt <- VariantAnnotation::alt(v)
    af <- VariantAnnotation::info(v)$AF
    df <- data.frame(
      chrom = as.character(GenomicRanges::seqnames(rr)), pos = GenomicRanges::start(rr),
      ref = as.character(VariantAnnotation::ref(v)),
      alt = vapply(as.list(alt), function(a) paste(as.character(a), collapse = ","), ""),
      af  = if (is.null(af)) NA else vapply(as.list(af), function(a) paste(a, collapse = ","), ""),
      id  = names(rr) %||% NA_character_, stringsAsFactors = FALSE)
  } else if (methods::is(v, "vcfR")) {                # vcfR
    fx <- v@fix
    df <- data.frame(chrom = fx[, "CHROM"], pos = as.integer(fx[, "POS"]), ref = fx[, "REF"],
                     alt = fx[, "ALT"], af = sub(".*(^|;)AF=([^;]+).*", "\\2", fx[, "INFO"]),
                     id = fx[, "ID"], stringsAsFactors = FALSE)
    df$af[!grepl("(^|;)AF=", fx[, "INFO"])] <- NA
  } else {
    df <- as.data.frame(v, stringsAsFactors = FALSE)
    nm <- .norm_colname(names(df))
    pick <- function(opts) { i <- match(opts, nm); i <- i[!is.na(i)][1]; if (is.na(i)) NULL else df[[i]] }
    chrom <- pick(.VAR_ALIASES$chrom)
    pos   <- pick(.VAR_ALIASES$pos)
    ref   <- pick(.VAR_ALIASES$ref)
    alt   <- pick(.VAR_ALIASES$alt)
    af    <- pick(.VAR_ALIASES$af)
    info  <- pick(.VAR_ALIASES$info)
    id    <- pick(.VAR_ALIASES$id)
    if (is.null(af) && !is.null(info)) {
      af <- ifelse(grepl("(^|;)AF=", info), sub(".*(^|;)AF=([^;]+).*", "\\2", info), NA)
    }
    if (is.null(chrom) || is.null(pos) || is.null(ref) || is.null(alt) || is.null(af)) {
      stop("variants need CHROM, POS, REF, ALT and AF columns (or an INFO column with AF=)")
    }
    df <- data.frame(chrom = as.character(chrom), pos = as.integer(pos), ref = as.character(ref),
                     alt = as.character(alt), af = as.character(af),
                     id = if (is.null(id)) NA_character_ else as.character(id), stringsAsFactors = FALSE)
  }
  # split multi-allelic rows
  multi <- grepl(",", df$alt)
  if (any(multi)) {
    a <- strsplit(df$alt[multi], ","); f <- strsplit(as.character(df$af[multi]), ",")
    k <- lengths(a)
    ex <- df[rep(which(multi), k), ]
    ex$alt <- unlist(a)
    ex$af  <- unlist(mapply(function(fv, n) rep_len(fv, n), f, k, SIMPLIFY = FALSE))
    df <- rbind(df[!multi, ], ex)
  }
  df$af <- suppressWarnings(as.numeric(df$af))
  df <- df[!is.na(df$pos) & !is.na(df$af), ]
  df$id[df$id %in% c(".", "")] <- NA
  df$key <- .chr_key(df$chrom)
  df <- df[order(df$key, df$pos), ]
  split(df, df$key)
}

.chr_key <- function(x) { x <- sub("^chr", "", x, ignore.case = TRUE); ifelse(x %in% c("M", "MT"), "MT", x) }

# column names accepted for variant tables (after .norm_colname: lower case,
# no leading "#", spaces/dots/dashes -> "_"). First match wins.
.VAR_ALIASES <- list(
  chrom = c("chrom", "chr", "chromosome", "seqnames", "seqname", "contig", "chrom_hg19", "chr_hg19"),
  pos   = c("pos", "position", "start", "bp", "pos_hg19", "position_grch37"),
  ref   = c("ref", "reference", "ref_allele", "reference_allele"),
  alt   = c("alt", "alternate", "alt_allele", "alternate_allele"),
  af    = c("af", "allele_frequency", "gnomad_af", "af_joint", "af_exomes", "af_genomes",
            "gnomad_exomes_af", "gnomad_genomes_af", "af_popmax", "maf", "freq"),
  info  = c("info"),
  id    = c("id", "rsid", "rsids", "rs", "variant_id", "gnomad_id"))
.norm_colname <- function(x) gsub("[ .\\-]+", "_", tolower(trimws(gsub('"', "", sub("^#+", "", x)))))

#' Read a gnomAD (or other) variant TEXT file, keeping only rows inside the
#' given regions so whole-genome files never have to fit in memory.
#' Handles .gz, tab or comma separated, a plain header row (any of the
#' .VAR_ALIASES names, e.g. a gnomAD browser export with Chromosome /
#' Position / Reference / Alternate / Allele Frequency) or a VCF written as
#' text (##meta lines, #CHROM header, AF= inside INFO).
#'
#' @param path     text file
#' @param regions  data.frame(chrom, start, end) to keep (see mane_regions());
#'                 NULL keeps everything
#' @param cols     optional manual mapping, e.g.
#'                 c(chrom = "CHR", pos = "BP", ref = "A1", alt = "A2", af = "AF_nfe")
#' @param chunk    lines per chunk
#' If the file is bgzipped and has a tabix index (path.tbi) and Rsamtools is
#' installed, only the regions are read through the index - use this for
#' genome-wide gnomAD files (bgzip file.vcf; tabix -p vcf file.vcf.gz).
#' @return data.frame(CHROM, POS, ID, REF, ALT, AF), ready for make_mane_pages(variants = ...)
read_variants_text <- function(path, regions = NULL, cols = NULL, chunk = 500000L) {
  con <- if (grepl("\\.b?gz$", path)) gzfile(path, "r") else file(path, "r")
  on.exit(close(con))
  # header: skip VCF "##" meta lines; header is the first remaining line
  repeat {
    h <- readLines(con, n = 1, warn = FALSE)
    if (!length(h)) stop("no header line found in ", path)
    if (!startsWith(h, "##") && nzchar(trimws(h))) break
  }
  sep <- if (grepl("\t", h, fixed = TRUE)) "\t" else if (grepl(",", h, fixed = TRUE)) "," else ""
  split_line <- function(x) {
    if (sep == ",") return(utils::read.table(text = x, sep = ",", quote = "\"", header = FALSE, colClasses = "character",
                                             comment.char = "", fill = TRUE, na.strings = character(0)))
    f <- if (sep == "\t") strsplit(x, "\t", fixed = TRUE) else strsplit(trimws(x), "[ ]+")
    n <- max(lengths(f)); f <- lapply(f, `length<-`, n)
    as.data.frame(do.call(rbind, f), stringsAsFactors = FALSE)
  }
  hdr <- unlist(split_line(h)[1, ]); nm <- .norm_colname(hdr)
  find <- function(key) {
    if (!is.null(cols) && !is.na(cols[key])) {
      i <- match(.norm_colname(cols[[key]]), nm)
      if (is.na(i)) stop("column '", cols[[key]], "' not found in ", path, " (have: ", paste(hdr, collapse = ", "), ")")
      return(i)
    }
    i <- match(.VAR_ALIASES[[key]], nm); i[!is.na(i)][1]
  }
  ix <- vapply(c("chrom", "pos", "ref", "alt", "af", "info", "id"), find, 1L)
  if (is.na(ix["af"]) && is.na(ix["info"])) stop("no allele-frequency column found in ", path,
    " (have: ", paste(hdr, collapse = ", "), "); pass cols = c(af = \"...\")")
  if (any(is.na(ix[c("chrom", "pos", "ref", "alt")]))) stop("need chromosome, position, ref and alt columns in ", path,
    " (have: ", paste(hdr, collapse = ", "), "); pass cols = c(chrom = ..., pos = ..., ref = ..., alt = ...)")
  message(sprintf("Reading %s: chrom=%s pos=%s ref=%s alt=%s af=%s", basename(path), hdr[ix["chrom"]], hdr[ix["pos"]],
                  hdr[ix["ref"]], hdr[ix["alt"]], if (is.na(ix["af"])) "INFO:AF" else hdr[ix["af"]]))

  reg <- NULL
  if (!is.null(regions)) {                       # merged intervals per chromosome
    regions$key <- .chr_key(regions$chrom)
    reg <- lapply(split(regions, regions$key), function(r) {
      r <- r[order(r$start), ]; s <- r$start; e <- r$end; os <- s[1]; oe <- e[1]; S <- c(); E <- c()
      for (k in seq_along(s)[-1]) if (s[k] <= oe + 1) oe <- max(oe, e[k]) else { S <- c(S, os); E <- c(E, oe); os <- s[k]; oe <- e[k] }
      list(s = c(S, os), e = c(E, oe))
    })
  }
  in_reg <- function(chrom, pos) {               # TRUE for positions inside the regions
    keep <- !is.na(pos)
    if (is.null(reg)) return(keep)
    key <- .chr_key(chrom); keep <- keep & key %in% names(reg)
    for (k in unique(key[keep])) {
      i <- which(keep & key == k); r <- reg[[k]]
      j <- findInterval(pos[i], r$s)
      keep[i] <- j > 0 & pos[i] <= r$e[pmax(j, 1)]
    }
    keep
  }
  # one field of every line by regex (C-level), so lines outside the regions
  # are dropped before the expensive split into columns
  field <- function(l, k) {
    if (sep == "\t") sub(sprintf("^(?:[^\t]*\t){%d}([^\t]*).*$", k - 1), "\\1", l, perl = TRUE)
    else sub(sprintf("^\\s*(?:\\S+\\s+){%d}(\\S+).*$", k - 1), "\\1", l, perl = TRUE)
  }
  process <- function(l) {
    l <- l[nzchar(l) & !startsWith(l, "#")]
    if (!length(l)) return(NULL)
    if (sep != "," && !is.null(reg)) {
      pos0 <- suppressWarnings(as.integer(field(l, ix["pos"])))
      l <- l[in_reg(field(l, ix["chrom"]), pos0)]
      if (!length(l)) return(NULL)
    }
    f <- split_line(l)
    pos <- suppressWarnings(as.integer(gsub(",", "", f[[ix["pos"]]])))
    keep <- in_reg(f[[ix["chrom"]]], pos)
    if (!any(keep)) return(NULL)
    f <- f[keep, , drop = FALSE]
    af <- if (!is.na(ix["af"])) f[[ix["af"]]] else {
      inf <- f[[ix["info"]]]; ifelse(grepl("(^|;)AF=", inf), sub(".*(^|;)AF=([^;]+).*", "\\2", inf), NA) }
    data.frame(CHROM = f[[ix["chrom"]]], POS = pos[keep],
      ID = if (is.na(ix["id"])) NA_character_ else f[[ix["id"]]],
      REF = f[[ix["ref"]]], ALT = f[[ix["alt"]]], AF = af, stringsAsFactors = FALSE)
  }
  out <- list(); n_read <- 0
  tbi <- paste0(path, ".tbi")
  if (!is.null(reg) && file.exists(tbi) && requireNamespace("Rsamtools", quietly = TRUE)) {
    # bgzipped + tabix-indexed: read only the regions (fast for genome-wide files)
    idx_names <- Rsamtools::headerTabix(path)$seqnames
    gr <- do.call(rbind, lapply(names(reg), function(k) {
      nm <- idx_names[.chr_key(idx_names) == k][1]
      if (is.na(nm)) NULL else data.frame(chr = nm, s = pmax(1, reg[[k]]$s), e = reg[[k]]$e)
    }))
    if (!is.null(gr) && nrow(gr)) {
      message("  using the tabix index: ", fmt_int(nrow(gr)), " regions")
      lines <- unlist(Rsamtools::scanTabix(path, param = GenomicRanges::GRanges(gr$chr, IRanges::IRanges(gr$s, gr$e))),
                      use.names = FALSE)
      n_read <- length(lines)
      for (st in seq(1, max(1, length(lines)), by = chunk))
        out[[length(out) + 1]] <- process(lines[st:min(st + chunk - 1, length(lines))])
    }
  } else {
    repeat {
      l <- readLines(con, n = chunk, warn = FALSE)
      if (!length(l)) break
      n_read <- n_read + length(l)
      out[[length(out) + 1]] <- process(l)
    }
  }
  out <- Filter(Negate(is.null), out)
  res <- if (length(out)) do.call(rbind, out) else
    data.frame(CHROM = character(0), POS = integer(0), ID = character(0), REF = character(0), ALT = character(0), AF = character(0))
  message(sprintf("  %s lines read, %s variants kept", fmt_int(n_read), fmt_int(nrow(res))))
  res
}

#' Regions covered by the pages (each transcript +/- `pad` bp), for
#' read_variants_text(regions = ...). Use pad >= intron_max of make_mane_pages().
mane_regions <- function(gtf, tags = c("MANE_Select", "MANE_Plus_Clinical"), pad = 1100) {
  df <- normalize_gtf(gtf)
  tag_rx <- paste0("(^|[,;| ])(", paste(tags, collapse = "|"), ")($|[,;| ])")
  df <- df[df$feature == "exon" & !is.na(df$tag) & grepl(tag_rx, df$tag), ]
  s <- tapply(df$start, df$transcript_id, min); e <- tapply(df$end, df$transcript_id, max)
  ch <- tapply(df$seqname, df$transcript_id, `[`, 1)
  data.frame(chrom = as.character(ch), start = as.numeric(s) - pad, end = as.numeric(e) + pad, stringsAsFactors = FALSE)
}

#' Variants overlapping a transcript locus (+ embedded flanks), in the
#' transcript-oriented "u" coordinate (u = 1 at the transcript 5' end).
variants_for_model <- function(vby, m, get_seq = NULL) {
  empty <- data.frame(u0 = integer(0), len = integer(0), label = character(0), af = numeric(0), id = character(0))
  if (is.null(vby)) return(empty)
  d <- vby[[.chr_key(m$chr)]]
  if (is.null(d)) return(empty)
  lo <- m$tx_start - m$embed - 50; hi <- m$tx_end + m$embed       # 50: allow long REF alleles
  a <- findInterval(lo, d$pos, left.open = TRUE) + 1L; b <- findInterval(hi, d$pos)   # lo <= pos <= hi
  if (a > b) return(empty)
  d <- d[a:b, ]
  L <- pmax(nchar(d$ref), 1L)
  u0 <- if (m$strand == "-") m$tx_end - (d$pos + L - 1L) + 1L else d$pos - m$tx_start + 1L
  # ANNOVAR writes indels with "-" (deletion: pos = first deleted base; insertion:
  # pos = base before the insert). Labels/gnomAD IDs need VCF style with the
  # preceding reference base: 13-100-AT-A, 13-100-A-AT. Highlighting keeps pos/L.
  vpos <- d$pos; vref <- d$ref; valt <- d$alt
  del <- valt == "-" & vref != "-"; ins <- vref == "-" & valt != "-"
  if (!is.null(get_seq) && any(del | ins)) {
    anc <- ifelse(del, vpos - 1L, vpos)
    rs <- min(anc[del | ins]); re <- max(anc[del | ins])
    sq <- toupper(tryCatch(get_seq(m$chr, rs, re), error = function(e) ""))
    if (nzchar(sq)) {
      b <- substring(sq, anc - rs + 1L, anc - rs + 1L)
      vref[del] <- paste0(b[del], vref[del]); valt[del] <- b[del]; vpos[del] <- anc[del]
      valt[ins] <- paste0(b[ins], valt[ins]); vref[ins] <- b[ins]
    }
  }
  data.frame(u0 = as.integer(u0), len = as.integer(L),
             label = paste(sub("^chr", "", d$chrom, ignore.case = TRUE), vpos, vref, valt, sep = "-"),
             af = d$af, id = d$id, stringsAsFactors = FALSE)
}

# ---- transcript model --------------------------------------------------------

#' Build everything needed for one page from the GTF rows of one transcript.
#' `embed` = bp of intron end / flank sequence embedded in the page (the
#' dropdown on the Exons view can show up to this much).
build_transcript_model <- function(rows, get_seq, flank = 50, embed = 1000) {
  embed <- max(embed, flank)
  ex <- rows[rows$feature == "exon", , drop = FALSE]
  if (!nrow(ex)) stop("no exon rows")
  strand <- ex$strand[1]; chr <- ex$seqname[1]; minus <- strand == "-"
  ex <- ex[order(if (minus) -ex$start else ex$start), , drop = FALSE]
  n_ex <- nrow(ex)
  len <- ex$end - ex$start + 1L
  cum <- c(0L, cumsum(len))
  tx_s <- min(ex$start); tx_e <- max(ex$end)
  to_u <- function(g) if (minus) tx_e - g + 1L else g - tx_s + 1L

  # fetch the whole locus once, then cut locally
  loc_s <- tx_s - embed; loc_e <- tx_e + embed
  locus <- get_seq(chr, loc_s, loc_e)
  gs <- function(s, e) {
    if (e < s) return("")
    x <- substr(locus, s - loc_s + 1, e - loc_s + 1)
    if (minus) revcomp(x) else x
  }
  g2c <- function(g) {
    i <- which(ex$start <= g & ex$end >= g)[1]
    if (is.na(i)) return(NA_integer_)
    as.integer(if (minus) cum[i] + ex$end[i] - g + 1 else cum[i] + g - ex$start[i] + 1)
  }

  ex_seq <- vapply(seq_len(n_ex), function(i) gs(ex$start[i], ex$end[i]), "")
  cdna <- paste(ex_seq, collapse = "")
  n <- nchar(cdna)

  # CDS in cDNA coordinates (stop codon included). Ensembl/GENCODE CDS rows
  # exclude the stop codon, so stop_codon rows are merged in.
  cr <- rows[rows$feature %in% c("CDS", "start_codon", "stop_codon"), , drop = FALSE]
  cds_s <- cds_e <- NA_integer_; protein <- ""; protein_id <- NA
  if (nrow(cr)) {
    gmin <- min(cr$start); gmax <- max(cr$end)
    cds_s <- g2c(if (minus) gmax else gmin)
    cds_e <- g2c(if (minus) gmin else gmax)
    cds_rows <- rows[rows$feature == "CDS", , drop = FALSE]
    if (nrow(cds_rows)) {
      first <- cds_rows[order(if (minus) -cds_rows$end else cds_rows$start), ][1, ]
      if (!is.na(first$frame) && first$frame > 0 &&
          (if (minus) first$end == gmax else first$start == gmin)) cds_s <- cds_s + first$frame
      pid <- cds_rows$protein_id[!is.na(cds_rows$protein_id)][1]
      pv  <- cds_rows$protein_version[!is.na(cds_rows$protein_version)][1]
      if (!is.na(pid)) protein_id <- if (!is.na(pv) && !grepl("\\.[0-9]+$", pid)) paste0(pid, ".", pv) else pid
    }
    cds <- substr(cdna, cds_s, cds_e)
    # GTF without stop_codon rows: extend over a following stop codon
    if (nchar(cds) %% 3 == 0 && !grepl("^(TAA|TAG|TGA)$", substr(cds, nchar(cds) - 2, nchar(cds))) &&
        cds_e + 3 <= n && grepl("^(TAA|TAG|TGA)$", substr(cdna, cds_e + 1, cds_e + 3))) {
      cds_e <- cds_e + 3L; cds <- substr(cdna, cds_s, cds_e)
    }
    protein <- translate_dna(cds)
    sec <- rows[rows$feature == "Selenocysteine", , drop = FALSE]
    if (nrow(sec)) for (g in (if (minus) sec$end else sec$start)) {
      k <- (g2c(g) - cds_s) %/% 3 + 1
      if (!is.na(k) && k <= nchar(protein)) substr(protein, k, k) <- "U"
    }
  }
  has_cds <- !is.na(cds_s)

  # exon phases (Ensembl convention, -1 = non-coding boundary)
  ph <- function(cpos, end_side) {
    if (!has_cds || cpos < cds_s || cpos > cds_e) return(-1L)
    if (end_side) as.integer((cpos - cds_s + 1) %% 3) else as.integer((cpos - cds_s) %% 3)
  }
  exon_ids <- ifelse(is.na(ex$exon_id), paste0("Exon ", seq_len(n_ex)),
                     ifelse(!is.na(ex$exon_version) & !grepl("\\.[0-9]+$", ex$exon_id),
                            paste0(ex$exon_id, ".", ex$exon_version), ex$exon_id))
  exons <- data.frame(
    no = seq_len(n_ex), id = exon_ids,
    gstart = ex$start, gend = ex$end,
    start5 = if (minus) ex$end else ex$start,       # transcript-oriented
    end3   = if (minus) ex$start else ex$end,
    cstart = cum[-length(cum)] + 1L, cend = cum[-1],
    ustart = to_u(if (minus) ex$end else ex$start), uend = to_u(if (minus) ex$start else ex$end),
    len = len, stringsAsFactors = FALSE)
  exons$sphase <- vapply(exons$cstart, ph, 0L, end_side = FALSE)
  exons$ephase <- vapply(exons$cend,   ph, 0L, end_side = TRUE)

  introns <- NULL
  if (n_ex > 1) {
    lo <- if (minus) ex$end[-1] + 1L else ex$end[-n_ex] + 1L
    hi <- if (minus) ex$start[-n_ex] - 1L else ex$start[-1] - 1L
    introns <- data.frame(no = seq_len(n_ex - 1), lo = lo, hi = hi,
                          start5 = if (minus) hi else lo, end3 = if (minus) lo else hi,
                          ustart = to_u(if (minus) hi else lo), uend = to_u(if (minus) lo else hi),
                          len = hi - lo + 1L, stringsAsFactors = FALSE)
    introns$seq <- vapply(seq_len(nrow(introns)), function(i) gs(introns$lo[i], introns$hi[i]), "")
  }
  up <- if (minus) gs(tx_e + 1, tx_e + embed) else gs(tx_s - embed, tx_s - 1)
  dn <- if (minus) gs(tx_s - embed, tx_s - 1) else gs(tx_e + 1, tx_e + embed)

  first_non_na <- function(x) x[!is.na(x) & nzchar(x)][1] %||% NA
  tid <- first_non_na(rows$transcript_id); tv <- first_non_na(rows$transcript_version)
  gid <- first_non_na(rows$gene_id);       gv <- first_non_na(rows$gene_version)
  addv <- function(id, v) if (!is.na(v) && !grepl("\\.[0-9]+$", id)) paste0(id, ".", v) else id
  tags <- unique(unlist(strsplit(paste(rows$tag[!is.na(rows$tag)], collapse = ","), "[,;|]")))
  tags <- trimws(tags[nzchar(tags)])

  list(
    transcript_id = addv(tid, tv), transcript_name = first_non_na(rows$transcript_name),
    gene_id = addv(gid, gv), gene_name = first_non_na(rows$gene_name) %||% gid,
    biotype = first_non_na(c(rows$transcript_biotype, rows$transcript_type)),
    gene_biotype = first_non_na(c(rows$gene_biotype, rows$gene_type)),
    ccds = first_non_na(rows$ccds_id), protein_id = protein_id, tags = tags,
    # Route A version fallback: this page is NOT the MANE transcript; these are the MANE IDs
    mane_refseq = first_non_na(rows$mane_refseq), mane_enst = first_non_na(rows$mane_enst),
    chr = chr, strand = strand, tx_start = tx_s, tx_end = tx_e,
    exons = exons, exon_seq = ex_seq, introns = introns, up = up, dn = dn,
    cdna = cdna, cds_s = cds_s, cds_e = cds_e, protein = protein, has_cds = has_cds,
    flank = flank, embed = embed)
}

# ---- HTML rendering ----------------------------------------------------------

# runs of identical class -> <span>s; class "" -> bare text
spans <- function(chars, cls) {
  if (!length(chars)) return("")
  r <- rle(cls); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  txt <- vapply(seq_along(ends), function(k) paste(chars[starts[k]:ends[k]], collapse = ""), "")
  txt <- html_escape(txt)
  paste(ifelse(nzchar(r$values), sprintf('<span class="%s">%s</span>', r$values, txt), txt),
        collapse = "")
}

render_cdna <- function(m, width = 60) {
  b <- strsplit(m$cdna, "")[[1]]; n <- length(b)
  exidx <- rep(seq_len(nrow(m$exons)), m$exons$len)
  pos <- seq_len(n)
  in_cds <- if (m$has_cds) pos >= m$cds_s & pos <= m$cds_e else rep(FALSE, n)
  codon <- ifelse(in_cds, (pos - m$cds_s) %/% 3, NA)
  cls <- paste0("x", (exidx - 1) %% 2,
                ifelse(in_cds, paste0(" c", codon %% 2), if (m$has_cds) " u" else ""))
  aa <- if (m$has_cds) strsplit(m$protein, "")[[1]] else character(0)
  is_res <- in_cds & (pos - m$cds_s) %% 3 == 1 & grepl("[A-Za-z]", aa[codon + 1])
  n_res <- cumsum(is_res)
  nw <- max(nchar(n), 4)
  num <- function(x) if (is.na(x)) strrep(" ", nw) else formatC(x, width = nw, format = "d")
  # line numbers (and their separating space) sit in .lnum spans so the
  # "Line numbers" markup box can hide them
  ln_l <- function(x) paste0('<span class="lnum">', num(x), ' </span>')
  ln_r <- function(x) paste0('<span class="lnum"> ', num(x), '</span>')
  out <- character(0)
  for (s in seq(1, n, by = width)) {
    e <- min(s + width - 1, n); idx <- s:e; padn <- width - length(idx)
    top <- sprintf('<span class="ln">%s<span class="bases" data-space="c" data-p="%d" data-n="%d">%s</span>%s%s\n</span>',
                   ln_l(s), s, length(idx), spans(b[idx], cls[idx]), strrep(" ", padn), ln_r(e))
    if (m$has_cds && any(in_cds[idx])) {
      ci <- idx[in_cds[idx]]
      cds_chars <- ifelse(in_cds[idx], b[idx], ".")
      cds_cls   <- ifelse(in_cds[idx], paste0("c", codon[idx] %% 2), "")
      off <- (idx - m$cds_s) %% 3
      pep_chars <- ifelse(!in_cds[idx], ".",
                          ifelse(off == 1, aa[codon[idx] + 1] %||% "-", "-"))
      pep_chars[is.na(pep_chars)] <- "-"
      l_cds <- if (in_cds[s]) s - m$cds_s + 1 else NA
      mid <- sprintf('<span class="ln cdsl">%s%s%s%s\n</span>', ln_l(l_cds),
                     spans(cds_chars, cds_cls), strrep(" ", padn), ln_r(max(ci) - m$cds_s + 1))
      # amino-acid letters on this line form a "p" (protein) container so
      # variant / selection highlights can be drawn on the translation too
      pep_txt <- html_escape(paste(pep_chars, collapse = ""))
      k_here <- codon[idx][in_cds[idx] & off == 1] + 1
      n_let <- sum(grepl("[A-Za-z*]", pep_chars))
      l_pep <- if (in_cds[s] && any(is_res[idx])) codon[s] + 1 else NA
      r_pep <- if (n_res[e] > 0) n_res[e] else NA
      pep_cont <- if (n_let > 0) sprintf('<span class="pep bases" data-space="p" data-p="%d" data-n="%d">%s</span>',
                                         min(k_here), n_let, pep_txt) else sprintf('<span class="pep">%s</span>', pep_txt)
      pep <- sprintf('<span class="ln pepl">%s%s%s%s\n</span>',
                     ln_l(l_pep), pep_cont, strrep(" ", padn), ln_r(r_pep))
    } else {
      dots <- paste0(ln_l(NA), strrep(".", length(idx)), "\n")
      mid <- sprintf('<span class="ln cdsl">%s</span>', dots)
      pep <- sprintf('<span class="ln pepl">%s</span>', dots)
    }
    out <- c(out, paste0(top, mid, pep, '<span class="ln gap">\n</span>'))
  }
  paste0('<pre class="seq cdna">', paste(out, collapse = ""), "</pre>")
}

# Protein view (Ensembl style): residues coloured by exon, residues whose
# codon spans an exon junction in red; numbered lines of `width` residues.
render_protein <- function(m, width = 60) {
  if (!m$has_cds) return('<p class="note">Non-coding transcript: no protein.</p>')
  aa <- strsplit(sub("\\*$", "", m$protein), "")[[1]]; n <- length(aa)
  if (n == 0) return('<p class="note">The CDS translates to a stop codon only: no protein to show.</p>')
  exidx <- rep(seq_len(nrow(m$exons)), m$exons$len)
  c0 <- m$cds_s + 3 * (seq_len(n) - 1)
  e1 <- exidx[c0]; e3 <- exidx[pmin(c0 + 2, length(exidx))]
  cls <- ifelse(e1 != e3, "so", paste0("x", (e1 - 1) %% 2))
  nw <- max(nchar(n), 4)
  num <- function(x) formatC(x, width = nw, format = "d")
  ln_l <- function(x) paste0('<span class="lnum">', num(x), ' </span>')
  ln_r <- function(x) paste0('<span class="lnum"> ', num(x), '</span>')
  out <- vapply(seq(1, n, by = width), function(s) {
    e <- min(s + width - 1, n); idx <- s:e
    sprintf('<span class="ln">%s<span class="bases" data-space="p" data-p="%d" data-n="%d">%s</span>%s%s\n</span>',
            ln_l(s), s, length(idx), spans(aa[idx], cls[idx]), strrep(" ", width - length(idx)), ln_r(e))
  }, "")
  paste0('<pre class="seq prot">', paste(out, collapse = ""), "</pre>")
}

wrap_spans <- function(chars, cls, width = 60) {
  if (!length(chars)) return("")
  st <- seq(1, length(chars), by = width)
  paste(vapply(st, function(s) { i <- s:min(s + width - 1, length(chars)); spans(chars[i], cls[i]) }, ""),
        collapse = "\n")
}

# Exon cells are rendered here; intron and flank cells are filled in by the
# page script so the "intron sequence shown" dropdown can change them.
# Line numbers follow Ensembl's exon export: the coding (CDS) position of the
# first and last coding base on each line; UTR-only, intron and flank lines
# get blank padding so the sequence stays aligned.
exon_lnw <- function(m) max(nchar(if (m$has_cds) m$cds_e - m$cds_s + 1 else 0), 4)
render_exons <- function(m, width = 60) {
  ph <- function(x) ifelse(x < 0, "-", as.character(x))
  nw <- exon_lnw(m)
  lnum <- function(x, left) paste0('<span class="lnum">', if (left) "" else " ",
                                   if (is.na(x)) strrep(" ", nw) else formatC(x, width = nw, format = "d"), if (left) " " else "", '</span>')
  exon_lines <- function(ch, kind, cp) {
    cdsp <- if (m$has_cds) ifelse(cp >= m$cds_s & cp <= m$cds_e, cp - m$cds_s + 1, NA) else rep(NA, length(cp))
    st <- seq(1, length(ch), by = width)
    paste(vapply(st, function(s) {
      i <- s:min(s + width - 1, length(ch)); k <- cdsp[i][!is.na(cdsp[i])]
      paste0(lnum(if (length(k)) min(k) else NA, TRUE), spans(ch[i], kind[i]),
             if (length(k)) paste0(strrep(" ", width - length(i)), lnum(max(k), FALSE)) else "")
    }, ""), collapse = "\n")
  }
  rows <- '<tr class="flank"><td></td><td>5\' upstream sequence</td><td></td><td></td><td></td><td></td><td></td><td class="sq" data-flank="up"></td></tr>'
  for (i in seq_len(nrow(m$exons))) {
    e <- m$exons[i, ]
    ch <- strsplit(m$exon_seq[i], "")[[1]]
    cp <- e$cstart:e$cend
    kind <- if (m$has_cds) ifelse(cp >= m$cds_s & cp <= m$cds_e, "cx", "ux") else rep("nx", length(cp))
    rows <- c(rows, sprintf(
      '<tr class="exon"><td>%d</td><td class="id">%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><pre class="seq bases" data-space="g" data-p="%d" data-n="%d">%s</pre></td></tr>',
      e$no, html_escape(e$id), fmt_int(e$start5), fmt_int(e$end3), ph(e$sphase), ph(e$ephase),
      fmt_int(e$len), e$ustart, e$len, exon_lines(ch, kind, cp)))
    if (!is.null(m$introns) && i <= nrow(m$introns)) {
      it <- m$introns[i, ]
      rows <- c(rows, sprintf(
        '<tr class="intron"><td></td><td>Intron %d-%d</td><td>%s</td><td>%s</td><td></td><td></td><td>%s</td><td class="sq" data-intron="%d"></td></tr>',
        i, i + 1, fmt_int(it$start5), fmt_int(it$end3), fmt_int(it$len), i - 1))
    }
  }
  rows <- c(rows, '<tr class="flank"><td></td><td>3\' downstream sequence</td><td></td><td></td><td></td><td></td><td></td><td class="sq" data-flank="dn"></td></tr>')
  paste0('<table class="exons"><thead><tr><th>No.</th><th>Exon / Intron</th><th>Start</th><th>End</th>',
         '<th>Start phase</th><th>End phase</th><th>Length</th><th>Sequence</th></tr></thead><tbody>',
         paste(rows, collapse = "\n"), "</tbody></table>")
}

# Genome build -> UCSC database, short label and matching Ensembl site
build_info <- function(build) {
  b <- tolower(build %||% "")
  if (grepl("hg19|grch37", b)) list(db = "hg19", short = "hg19", ens = "https://grch37.ensembl.org")
  else if (grepl("hg38|grch38", b)) list(db = "hg38", short = "hg38", ens = "https://www.ensembl.org")
  else list(db = NULL, short = NULL, ens = "https://www.ensembl.org")
}
ext_link <- function(url, text, title) sprintf('<a class="ext" href="%s" target="_blank" rel="noopener" title="%s">%s</a>',
                                               html_escape(url), html_escape(title), text)

# HGVS-style c./n. label for cDNA positions
c_dot <- function(m, pos) {
  if (!m$has_cds) return(paste0("n.", pos))
  ifelse(pos < m$cds_s, paste0("c.-", m$cds_s - pos),
         ifelse(pos > m$cds_e, paste0("c.*", pos - m$cds_e), paste0("c.", pos - m$cds_s + 1)))
}

#' Compare a page's cDNA (built from the genome FASTA) with the transcript's
#' own sequence (e.g. MANE ensembl_rna). Returns NULL when identical.
.is_fb <- function(m) isTRUE(!is.null(m$mane_refseq) && !is.na(m$mane_refseq))
.fb_word <- function(m) {   # "older"/"newer" version than MANE
  v <- suppressWarnings(as.integer(sub("^.*\\.", "", c(m$transcript_id, m$mane_refseq))))
  if (anyNA(v) || v[1] == v[2]) "other" else if (v[1] < v[2]) "older" else "newer"
}
.ref_key <- function(m) if (.is_fb(m)) m$mane_enst else m$transcript_id

ref_diff <- function(m, ref) {
  if (is.null(ref)) return(NULL)
  r <- ref[sub("\\.[0-9]+$", "", .ref_key(m))]
  if (!length(r) || is.na(r)) return(NULL)
  r <- toupper(unname(r))
  if (nchar(r) != nchar(m$cdna))
    return(list(pos = integer(0), len_diff = TRUE,
                txt = sprintf("length differs (transcript %s bp, genome %s bp)",
                              fmt_int(nchar(r)), fmt_int(nchar(m$cdna)))))
  a <- strsplit(m$cdna, "")[[1]]; b <- strsplit(r, "")[[1]]; d <- which(a != b)
  if (!length(d)) return(NULL)
  lab <- sprintf("%s (transcript %s, genome %s)", c_dot(m, d), b[d], a[d])
  list(pos = d, len_diff = FALSE,
       txt = paste0(length(d), " position", if (length(d) > 1) "s" else "", ": ",
                    paste(utils::head(lab, 10), collapse = "; "), if (length(d) > 10) "; ..." else ""))
}

# Page data for the script: sequences, coordinates, intron ends, variants.
page_json <- function(m, vars, bins = .AF_BINS, build = NULL, prov = NULL, refdiff = NULL) {
  bins_js <- paste(sprintf("[%s,%s,%s,%s]", format(bins$min, scientific = FALSE, trim = TRUE, drop0trailing = TRUE),
                           vapply(bins$bg, json_str, ""), vapply(bins$fg, json_str, ""),
                           vapply(bins$label, json_str, "")), collapse = ",")
  num <- function(x) ifelse(is.na(x), "null", format(x, scientific = FALSE, trim = TRUE))
  e <- m$exons
  ex_js <- paste(sprintf("[%s,%s,%d,%d,%d,%d,%d,%d]", vapply(e$id, json_str, ""),
                         vapply(m$exon_seq, json_str, ""), e$gstart, e$gend,
                         e$cstart, e$cend, e$ustart, e$uend), collapse = ",")
  in_js <- ""
  if (!is.null(m$introns)) {
    it <- m$introns; k <- m$embed
    head <- ifelse(it$len <= 2 * k, it$seq, substr(it$seq, 1, k))
    tail <- ifelse(it$len <= 2 * k, "", substr(it$seq, it$len - k + 1, it$len))
    in_js <- paste(sprintf("[%d,%d,%d,%s,%s]", it$ustart, it$uend, it$len,
                           vapply(head, json_str, ""), vapply(tail, json_str, "")), collapse = ",")
  }
  v_js <- ""
  if (nrow(vars)) {
    v_js <- paste(sprintf("[%d,%d,%s,%s,%s]", vars$u0, vars$len, vapply(vars$label, json_str, ""),
                          sprintf("%.6g", vars$af),
                          ifelse(is.na(vars$id), "null", vapply(ifelse(is.na(vars$id), "", vars$id), json_str, ""))),
                  collapse = ",")
  }
  sprintf(paste0('{"transcript_id":%s,"gene_name":%s,"location":%s,"chr":%s,"strand":%s,',
                 '"tx_s":%d,"tx_e":%d,"cds_s":%s,"cds_e":%s,"flank":%d,"embed":%d,"lnw":%d,',
                 '"cdna":%s,"cds":%s,"protein":%s,"up":%s,"dn":%s,"has_variants":%s,',
                 '"build":%s,"prov":%s,"tname":%s,"mane":%s,"refdiff":[%s],"lendiff":%s,',
                 '"bins":[%s],"exons":[%s],"introns":[%s],"variants":[%s]}'),
    json_str(m$transcript_id), json_str(m$gene_name),
    json_str(sprintf("%s:%d-%d:%s", m$chr, m$tx_start, m$tx_end, if (m$strand == "-") -1 else 1)),
    json_str(m$chr), json_str(m$strand), m$tx_start, m$tx_end, num(m$cds_s), num(m$cds_e),
    as.integer(m$flank), as.integer(m$embed), as.integer(exon_lnw(m)),
    json_str(m$cdna), json_str(if (m$has_cds) substr(m$cdna, m$cds_s, m$cds_e) else ""),
    json_str(m$protein), json_str(m$up), json_str(m$dn),
    if (isTRUE(attr(vars, "supplied"))) "true" else "false",
    json_str(build %||% ""), json_str(prov %||% ""), json_str(if (is.na(m$transcript_name) || m$transcript_name == m$transcript_id) "" else m$transcript_name),
    json_str(if (.is_fb(m)) m$mane_refseq else ""),
    paste(refdiff$pos, collapse = ","), if (isTRUE(refdiff$len_diff)) "true" else "false",
    bins_js, ex_js, in_js, v_js)
}

.PAGE_CSS <- r"---(
:root{--ink:#1d2733;--muted:#5d6b7a;--line:#d9dee4;--panel:#f4f6f8;--head:#24435f;--stripe:#eef0f2;
--exon2:#1044ee;--utr:#ffdf33;--codon:#fff9af;--pep:#dd0000;--fl:#008000;--ux:#cd6839;--cx:#1044ee;--in:#8a8a8a}
*{box-sizing:border-box}
body{margin:0;font:14px/1.45 system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;color:var(--ink);background:#fff}
header.top{background:var(--head);color:#fff;padding:10px 20px;display:flex;gap:16px;align-items:baseline;flex-wrap:wrap}
header.top a{color:#cfe0f2;text-decoration:none;font-size:13px}
header.top h1{font-size:18px;margin:0;font-weight:600}
main{padding:16px 20px 40px;max-width:1400px}
.summary{display:grid;grid-template-columns:max-content 1fr;gap:6px 18px;background:var(--panel);border:1px solid var(--line);padding:12px 16px;border-radius:4px}
.summary dt{font-weight:600;color:var(--muted)}
.summary dd{margin:0}
a.ext{color:#1864ab;text-decoration:none;border-bottom:1px dotted #1864ab}
a.ext:hover{border-bottom-style:solid}
a.ext::after{content:"\2197";font-size:.8em;margin-left:2px;opacity:.7}
.badge{display:inline-block;padding:1px 8px;margin:0 4px 2px 0;border-radius:3px;font-size:12px;background:#e3e8ee;color:#33414f}
.badge.mane{background:#24435f;color:#fff}
.varbar{display:flex;flex-wrap:wrap;gap:6px 16px;align-items:center;margin-top:10px;font-size:13px;padding:8px 14px;border:1px solid var(--line);border-radius:4px}
.varbar b{color:var(--muted)}
.varbar select{font:inherit;font-size:13px;padding:2px 4px}
.vsw{font-family:Menlo,Consolas,"Courier New",monospace;font-size:11px;padding:1px 5px;border-radius:2px;text-decoration:underline dotted 1.5px;white-space:nowrap}
.vlegend{display:inline-flex;flex-wrap:wrap;gap:3px;align-items:center}
.varbar.off .vlegend,.varbar.off .vfilter{opacity:.45}
nav.tabs{margin:18px 0 0;display:flex;gap:4px;border-bottom:2px solid var(--head)}
nav.tabs button{font:inherit;border:1px solid var(--line);border-bottom:none;background:var(--panel);padding:7px 16px;cursor:pointer;border-radius:4px 4px 0 0;color:var(--ink)}
nav.tabs button[aria-selected=true]{background:var(--head);color:#fff;border-color:var(--head)}
section.view{display:none;padding-top:12px}
section.view.on{display:block}
.toolbar{display:flex;flex-wrap:wrap;gap:8px 18px;align-items:center;margin:4px 0 12px}
.toolbar button{font:inherit;font-size:13px;padding:4px 10px;border:1px solid #9aa8b6;background:#fff;border-radius:3px;cursor:pointer}
.toolbar .dl button{margin-right:4px}
.toolbar button:hover{background:var(--panel)}
.key{display:grid;grid-template-columns:max-content 1fr;gap:6px 14px;font-size:13px;background:var(--panel);border:1px solid var(--line);padding:10px 14px;border-radius:4px;margin-bottom:12px;max-width:820px}
.key b{color:var(--muted);font-weight:600}
.key label{margin-right:14px;cursor:pointer;white-space:nowrap}
.key select{font:inherit;font-size:13px;padding:1px 4px}
.sw{display:inline-block;padding:0 6px;font-family:Menlo,Consolas,"Courier New",monospace;font-size:12px;border:1px solid var(--line)}
pre.seq{font:12.5px/1.35 Menlo,Consolas,"DejaVu Sans Mono","Courier New",monospace;margin:0;white-space:pre}
.cdna{overflow-x:auto;padding:4px 0}
.x1{color:var(--exon2)}
.u{background:var(--utr)}
.c1{background:var(--codon)}
.pep{color:var(--pep)}
.cdna.no-exons .x1{color:inherit}
.cdna.no-utr .u{background:none}
.cdna.no-codons .c1{background:none}
.cdna.no-trans .cdsl,.cdna.no-trans .pepl{display:none}
.prot{overflow-x:auto;padding:4px 0}
.so{color:#e00000}
.prot.no-exons .x1{color:inherit}
.prot.no-so .so{color:inherit}
.no-lnum .lnum{display:none}
.build{display:inline-block;background:#1864ab;color:#fff;font-weight:700;font-size:12px;border-radius:3px;padding:1px 8px;margin-left:6px;letter-spacing:.02em;vertical-align:1px}
header.top .build{background:#fff;color:var(--head)}
.strand-badge{display:inline-block;background:#fff3bf;color:#7a5b00;border:1px solid #f0d36b;font-size:12px;font-weight:600;border-radius:3px;padding:0 7px;margin-left:4px}
.nochr{color:#8a1c1c;font-size:12.5px;font-weight:600}
.verwarn{margin-top:10px;padding:10px 14px;border:2px solid #e8b10c;background:#fff9db;color:#5c4300;border-radius:4px;font-size:13.5px}
.badge.old{background:#fff3bf;color:#6b4e00;border:1px solid #e8b10c;font-weight:600}
.refwarn{margin-top:10px;padding:8px 14px;border:1px solid #f1aeb5;background:#fff5f5;color:#8a1c1c;border-radius:4px;font-size:13px}
footer.prov{margin-top:28px;padding-top:8px;border-top:1px solid var(--line);color:var(--muted);font-size:11.5px}
.info-i{display:inline-block;width:14px;height:14px;line-height:14px;text-align:center;border-radius:50%;background:#9aa8b6;color:#fff;font-size:10px;font-weight:700;font-style:italic;cursor:help}
form.goto{display:inline-flex;gap:4px;align-items:center}
form.goto input{font:13px Menlo,Consolas,"Courier New",monospace;padding:4px 8px;border:1px solid #9aa8b6;border-radius:3px;width:330px}
form.goto .gmsg{font-size:12px;color:#c92a2a;margin-left:4px}
#pickbar .warns{margin-top:4px;font-size:12.5px;color:#c92a2a;gap:4px 14px}
#pickbar .warns:empty{display:none}
#pickbar .warns .n{color:var(--muted)}
.chip{display:inline-flex;align-items:stretch;border:1px solid #c9d1da;background:#fff;border-radius:4px;margin:1px 4px 1px 0;font-size:12.5px;white-space:nowrap;overflow:hidden;vertical-align:middle}
.chip .k{background:#e6ebf0;color:#4a5866;font-size:11.5px;font-weight:600;padding:1px 7px;border-right:1px solid #c9d1da}
.chip .v{padding:1px 8px;font-variant-numeric:tabular-nums}
.chip .sub{color:var(--muted);font-size:11.5px}
.sep{color:#9aa8b6;margin:0 .35em}
/* Some browsers leave highlight pixels (underline/background edge) behind when a
   highlight is removed. Toggling the block background forces the whole block to
   repaint; the transparent underline makes the text paint area include underlines. */
main.hlr .seq,main.hlr table.exons{background-color:rgba(255,255,255,0)}
.bases{text-decoration:underline 3px transparent;text-underline-offset:auto}
.toolbar button.rtf{font-weight:600}
table.exons{border-collapse:collapse;font-size:13px}
table.exons th{position:sticky;top:0;background:var(--head);color:#fff;text-align:left;padding:6px 10px;font-weight:600;white-space:nowrap;z-index:1}
table.exons td{border-bottom:1px solid var(--line);padding:4px 10px;vertical-align:top;white-space:nowrap;background:#fff}
table.exons tr.alt td{background:var(--stripe)}
table.exons td:nth-child(3),table.exons td:nth-child(4),table.exons td:nth-child(7){text-align:right;font-variant-numeric:tabular-nums}
table.exons td:nth-child(5),table.exons td:nth-child(6){text-align:center}
table.exons tr.intron td{color:var(--muted)}
table.exons tr.flank td{color:var(--muted);font-style:italic}
table.exons tr.flank pre{font-style:normal}
.fl{color:var(--fl)}.ux{color:var(--ux)}.cx{color:var(--cx)}.nx{color:var(--ink)}.in{color:var(--in)}
table.exons.no-utr .ux{color:var(--cx)}
table.exons.no-introns tr.intron{display:none}
table.exons.no-flanks tr.flank{display:none}
.note{color:var(--muted);font-size:12px;margin-top:10px}
.hint{font-size:12px;color:var(--muted)}
.bases{cursor:pointer}
::highlight(sel){background-color:#a5d8ff;color:#000}
::highlight(anc){background-color:#1c7ed6;color:#fff}
::highlight(vsel){color:#e00000;background-color:#ffe066;text-decoration:underline solid 3px #e00000}
::highlight(vchk){color:#0000e0;background-color:#ffe066;text-decoration:underline solid 3px #0000e0}
#tip .sq{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px;border:1px solid rgba(255,255,255,.5)}
#tip{position:fixed;z-index:20;pointer-events:none;background:#1d2733;color:#fff;font-size:12px;line-height:1.4;padding:6px 9px;border-radius:4px;max-width:360px;display:none;box-shadow:0 2px 8px rgba(0,0,0,.25)}
#tip .af{color:#ffd43b}
#pickbar{position:fixed;left:0;right:0;bottom:0;background:#fff;border-top:2px solid var(--head);box-shadow:0 -4px 14px rgba(0,0,0,.12);display:none;z-index:10}
body.has-pick #pickbar{display:block}
#pickbar .dtab{position:absolute;left:20px;bottom:100%;display:flex;align-items:center;gap:10px;background:var(--head);color:#fff;
  padding:4px 12px 4px 6px;border-radius:6px 6px 0 0;font-size:12px;cursor:ns-resize;user-select:none;touch-action:none;max-width:calc(100vw - 40px)}
#pickbar .dtoggle{background:rgba(255,255,255,.15);border:none;color:#fff;border-radius:3px;width:24px;height:20px;padding:0;font-size:11px;cursor:pointer;line-height:20px}
#pickbar .dtoggle:hover{background:rgba(255,255,255,.3)}
#pickbar .dtitle{font-weight:600}
#pickbar .dsum{opacity:.85;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;font-variant-numeric:tabular-nums}
#pickbar .grip{opacity:.6;letter-spacing:-1px}
#pickbar .dbody{padding:10px 20px;max-height:80vh;overflow:auto}
#pickbar.collapsed .dbody{display:none}
#pickbar.collapsed{border-top-width:2px;box-shadow:none}
#pickbar .row{display:flex;gap:10px;align-items:center;flex-wrap:wrap;max-width:1400px}
#pickbar .info{font-size:13px;color:var(--ink);font-variant-numeric:tabular-nums}
#pickbar .info b{color:var(--head)}
#pickbar textarea{flex:1 1 420px;min-width:220px;height:48px;resize:vertical;font:12.5px/1.35 Menlo,Consolas,"Courier New",monospace;padding:6px 8px;border:1px solid #9aa8b6;border-radius:3px;word-break:break-all}
#pickbar button{font:inherit;font-size:13px;padding:5px 12px;border:1px solid #9aa8b6;background:#fff;border-radius:3px;cursor:pointer}
#pickbar button.primary{background:var(--head);border-color:var(--head);color:#fff}
#pickbar button.mark{color:#c92a2a;border-color:#e8a0a0;font-weight:600}
.markvar a{color:#c92a2a}
.mvb{color:#0000e0}.mvb a{color:#3399ff}
#pickbar label{font-size:13px;white-space:nowrap;cursor:pointer}
)---"

.PAGE_JS <- r"---(
(function(){
"use strict";
var D=JSON.parse(document.getElementById("seqdata").textContent);
var $=function(s,r){return (r||document).querySelector(s)};
var $$=function(s,r){return [].slice.call((r||document).querySelectorAll(s))};
var HL=!!(window.CSS&&CSS.highlights&&window.Highlight);
var fmt=function(n){return n.toLocaleString("en-US")};
var isBase=function(ch){return /[A-Za-z*]/.test(ch)};   // "*" = stop / untranslated Sec residue
function revcomp(s){var m={A:"T",C:"G",G:"C",T:"A",N:"N",a:"t",c:"g",g:"c",t:"a",n:"n"};
  return s.split("").reverse().map(function(x){return m[x]||x}).join("")}
function esc(s){return String(s).replace(/[&<>"]/g,function(c){return{"&":"&amp;","<":"&lt;",">":"&gt;","\"":"&quot;"}[c]})}

/* ---------- coordinates ----------
   c = cDNA position, u = position along the unspliced transcript
   (u=1 at the transcript 5' end, <1 upstream, >L downstream) */
var L=D.tx_e-D.tx_s+1, NC=D.cdna.length;
var EX=D.exons.map(function(e){return{id:e[0],seq:e[1],gs:e[2],ge:e[3],cs:e[4],ce:e[5],us:e[6],ue:e[7]}});
var IN=D.introns.map(function(i){return{us:i[0],ue:i[1],len:i[2],head:i[3],tail:i[4]}});
function exonU(u){for(var i=0;i<EX.length;i++)if(u>=EX[i].us&&u<=EX[i].ue)return i;return -1}
function exonC(c){for(var i=0;i<EX.length;i++)if(c>=EX[i].cs&&c<=EX[i].ce)return i;return -1}
function u2c(u){var i=exonU(u);return i<0?null:EX[i].cs+(u-EX[i].us)}
function c2u(c){var i=exonC(c);return i<0?null:EX[i].us+(c-EX[i].cs)}
function u2g(u){return D.strand==="-"?D.tx_e-u+1:D.tx_s+u-1}
function cdotC(c){if(D.cds_s===null)return "n."+c;
  if(c<D.cds_s)return "c.-"+(D.cds_s-c); if(c>D.cds_e)return "c.*"+(c-D.cds_e); return "c."+(c-D.cds_s+1)}
function cdotU(u){
  var c=u2c(u); if(c!==null)return cdotC(c);
  if(u<1)return D.cds_s===null?"n.-"+(1-u):"c.-"+(D.cds_s-u);
  if(u>L)return D.cds_s===null?"n.*"+(u-L):"c.*"+(NC+(u-L)-D.cds_e);
  for(var i=0;i<EX.length-1;i++)if(u>EX[i].ue&&u<EX[i+1].us){
    var d1=u-EX[i].ue,d2=EX[i+1].us-u;
    return d1<=d2?cdotC(EX[i].ce)+"+"+d1:cdotC(EX[i+1].cs)+"-"+d2}
  return "?";
}
function locU(u){var i=exonU(u);if(i>=0)return "exon "+(i+1);if(u<1)return "5′ flank";if(u>L)return "3′ flank";
  for(var k=0;k<EX.length-1;k++)if(u>EX[k].ue&&u<EX[k+1].us)return "intron "+(k+1);return "?"}
function gtxt(ua,ub){var a=u2g(ua),b=u2g(ub);return D.chr+":"+(a===b?fmt(a):fmt(Math.min(a,b))+"-"+fmt(Math.max(a,b)))+" ("+D.strand+")"}

/* unspliced sequence from the embedded segments (null if not embedded) */
var SEG=[{u0:1-D.up.length,s:D.up}];
EX.forEach(function(e){SEG.push({u0:e.us,s:e.seq})});
IN.forEach(function(i){SEG.push({u0:i.us,s:i.head});if(i.tail)SEG.push({u0:i.ue-i.tail.length+1,s:i.tail})});
SEG.push({u0:L+1,s:D.dn});
SEG.sort(function(a,b){return a.u0-b.u0});   // 5'->3' order, so intron+exon selections read correctly
function useq(a,b){var out="",got=0;
  SEG.forEach(function(g){var s=Math.max(a,g.u0),e=Math.min(b,g.u0+g.s.length-1);
    if(s<=e){out+=g.s.slice(s-g.u0,e-g.u0+1);got+=e-s+1}});
  return got===b-a+1?out.toUpperCase():null}

/* ---------- variants ---------- */
var VAR=D.variants.map(function(v){return{u0:v[0],u1:v[0]+v[1]-1,label:v[2],af:v[3],id:v[4]}})
  .sort(function(a,b){return a.u0-b.u0});
/* AF bands (from R): [[min AF, background, text colour, label], ...] highest first */
var BINS=D.bins;
var vmin=$("#vmin")?+$("#vmin").value:0, vshow=false;                    // off until the variant checkbox is ticked
function vbin(af){for(var i=0;i<BINS.length;i++)if(af>=BINS[i][0])return i;return BINS.length-1}
function vvis(v){return vshow&&v.af>=vmin}
function varsIn(ua,ub){return VAR.filter(function(v){return vvis(v)&&v.u1>=ua&&v.u0<=ub})}
function afTxt(af){return af>=0.01?(af*100).toFixed(2)+"%":af>=0.0001?(af*100).toFixed(3)+"%":af.toExponential(2)}

/* ---------- tabs ---------- */
function show(id){
  $$("section.view").forEach(function(s){s.classList.toggle("on",s.id===id)});
  $$("nav.tabs button").forEach(function(b){b.setAttribute("aria-selected",b.dataset.view===id)});
}
$$("nav.tabs button").forEach(function(b){
  b.addEventListener("click",function(){
    if(b.dataset.view!==currentView())clearVariantColours();
    show(b.dataset.view);history.replaceState(null,"","#"+b.dataset.view);
    if(typeof paintSel==="function")paintSel()})});
/* switching tabs unticks the variant box and clears its colours; some browsers
   left the colours on the wrong elements after a tab change */
function clearVariantColours(){var c=document.getElementById("vshow");if(!c||!c.checked)return;
  c.checked=false;c.dispatchEvent(new Event("change"))}
function fromHash(){var h=location.hash.slice(1);show((h==="exons"||h==="protein")&&document.getElementById(h)?h:"cdna")}
function currentView(){var s=$("section.view.on");return s?s.id:"cdna"}
fromHash(); window.addEventListener("hashchange",fromHash);

/* ---------- striping (white / grey, recomputed when rows are hidden) ---------- */
function restripe(){var k=0;$$("table.exons tbody tr").forEach(function(r){
  if(r.offsetParent===null&&getComputedStyle(r).display==="none")return; r.classList.toggle("alt",(k++)%2===1)})}

/* ---------- markup toggles ---------- */
$$("input[data-toggle]").forEach(function(c){
  var t=$(c.dataset.target);
  t.classList.toggle(c.dataset.toggle,!c.checked);   // browsers may restore an unticked box on reload / Back
  c.addEventListener("change",function(){t.classList.toggle(c.dataset.toggle,!c.checked);restripe()})});

/* ---------- helpers ---------- */
function wrap(s){return s.replace(/(.{60})/g,"$1\n").replace(/\n$/,"")}
function download(name,text,type){var a=document.createElement("a");
  a.href=URL.createObjectURL(new Blob([text],{type:type||"text/plain"}));a.download=name;
  document.body.appendChild(a);a.click();a.remove()}

/* ---------- Exons view: intron / flank cells (dropdown) ---------- */
var LPAD='<span class="lnum">'+" ".repeat((D.lnw||4)+1)+"</span>";
function padLines(html){return LPAD+html.replace(/\n/g,"\n"+LPAD)}
function box(seq,p,cls){return '<span class="bases '+cls+'" data-space="g" data-p="'+p+'" data-n="'+seq.length+'">'+wrap(seq.toLowerCase())+"</span>"}
var DOTS="..........";
function renderFlex(n){
  var m=Math.min(Math.max(D.flank,n),D.up.length);
  var up=$('td[data-flank="up"]'),dn=$('td[data-flank="dn"]');
  var su=D.up.slice(D.up.length-m), sd=D.dn.slice(0,Math.min(m,D.dn.length));
  if(up)up.innerHTML='<pre class="seq">'+padLines(DOTS+(m>50?"\n":"")+box(su,1-su.length,"fl"))+"</pre>";
  if(dn)dn.innerHTML='<pre class="seq">'+padLines(box(sd,L+1,"fl")+(m>50?"\n":"")+DOTS)+"</pre>";
  $$("td[data-intron]").forEach(function(td){
    var it=IN[+td.dataset.intron], html;
    if(it.len<=2*n){html=box(useq(it.us,it.ue)||"",it.us,"in")}
    else{var h=useq(it.us,it.us+n-1),t=useq(it.ue-n+1,it.ue),br=n>25?"\n":"";
      html=box(h,it.us,"in")+br+'<span class="in">'+DOTS+"</span>"+br+box(t,it.ue-n+1,"in")}
    td.innerHTML='<pre class="seq">'+padLines(html)+"</pre>";
  });
  paintVariants();paintSel();
}
var isel=$("#intronsize");
if(isel)isel.addEventListener("change",function(){renderFlex(+isel.value)});

/* ---------- highlighting helpers ---------- */
var LC=new WeakMap();
function letters(cont){var a=LC.get(cont);if(a)return a;a=[];
  var w=document.createTreeWalker(cont,NodeFilter.SHOW_TEXT),n;
  while((n=w.nextNode()))for(var i=0;i<n.data.length;i++)if(isBase(n.data[i]))a.push([n,i]);
  LC.set(cont,a);return a}
function addRange(hl,cont,k0,k1){var a=letters(cont);if(!a[k0]||!a[k1])return;
  var r=document.createRange();r.setStart(a[k0][0],a[k0][1]);r.setEnd(a[k1][0],a[k1][1]+1);hl.add(r)}
/* intervals {c:[[a,b]...],u:[[a,b]...],p:[[a,b]...]} painted onto every
   container of the matching space (c = cDNA, g/u = genomic, p = protein) */
function spaceKey(sp){return sp==="g"?"u":sp}
function paint(hl,iv){
  $$(".bases").forEach(function(cont){
    var list=iv[spaceKey(cont.dataset.space)];if(!list||!list.length)return;
    var p=+cont.dataset.p,n=+cont.dataset.n,q=p+n-1;
    for(var j=0;j<list.length;j++){var a=Math.max(list[j][0],p),b=Math.min(list[j][1],q);
      if(a<=b)addRange(hl,cont,a-p,b-p)}
  });
}
function cFromU(ua,ub){var o=[];EX.forEach(function(e){var a=Math.max(ua,e.us),b=Math.min(ub,e.ue);
  if(a<=b)o.push([e.cs+a-e.us,e.cs+b-e.us])});return o}
function uFromC(ca,cb){var o=[];EX.forEach(function(e){var a=Math.max(ca,e.cs),b=Math.min(cb,e.ce);
  if(a<=b)o.push([e.us+a-e.cs,e.us+b-e.cs])});return o}
/* residues (protein positions) whose codons overlap cDNA interval [ca,cb] */
var NAA=D.cds_s===null?0:D.protein.replace(/\*$/,"").length;
var NP=D.cds_s===null?0:D.protein.length;      // residues incl. the final stop "*" (shown on the cDNA tab only)
function pFromC(ca,cb){if(D.cds_s===null)return [];
  var a=Math.max(ca,D.cds_s),b=Math.min(cb,D.cds_s+3*NP-1);if(a>b)return [];
  return [[Math.floor((a-D.cds_s)/3)+1,Math.floor((b-D.cds_s)/3)+1]]}
function pFromCs(cs){var o=[];cs.forEach(function(x){pFromC(x[0],x[1]).forEach(function(y){o.push(y)})});return o}
function cOfP(a,b){return [D.cds_s+3*(a-1),D.cds_s+3*b-1]}
function ivOf(s){
  if(s.space==="c")return{c:[[s.a,s.b]],u:uFromC(s.a,s.b),p:pFromC(s.a,s.b)};
  if(s.space==="p"){var c=cOfP(s.a,s.b);return{p:[[s.a,s.b]],c:[c],u:uFromC(c[0],c[1])}}
  var cs=cFromU(s.a,s.b);return{u:[[s.a,s.b]],c:cs,p:pFromCs(cs)}}
/* u-range a mouse position stands for (a residue stands for its codon) */
function uRangeOf(p){if(p.space==="c"){var u=c2u(p.pos);return [u,u]}
  if(p.space==="p"){var c=cOfP(p.pos,p.pos),us=uFromC(c[0],c[1]);return [us[0][0],us[us.length-1][1]]}
  return [p.pos,p.pos]}
/* Highlights are kept as one persistent object per name and cleared/refilled;
   an empty one is removed from the registry. Some browsers do not repaint
   when a highlight is merely replaced, which left colours behind. */
var REG={};
function setHL(name,hl,prio){if(!HL)return;
  var h=REG[name];if(!h){h=REG[name]=new Highlight();h.priority=prio}
  h.clear();hl.forEach(function(r){h.add(r)});
  if(h.size){if(CSS.highlights.get(name)!==h)CSS.highlights.set(name,h)}else CSS.highlights.delete(name)}
function nudgeRepaint(){var m=$("main");m.classList.toggle("hlr");
  requestAnimationFrame(function(){m.classList.toggle("hlr")})}

function paintVariants(){
  if(HL){
    var hs=BINS.map(function(){return new Highlight()}),iv=BINS.map(function(){return{c:[],u:[],p:[]}});
    VAR.forEach(function(v){if(!vvis(v))return;var b=vbin(v.af),cs=cFromU(v.u0,v.u1);
      iv[b].u.push([v.u0,v.u1]);cs.forEach(function(x){iv[b].c.push(x)});pFromCs(cs).forEach(function(x){iv[b].p.push(x)})});
    hs.forEach(function(h,i){paint(h,iv[i]);setHL("v"+i,h,3)});
    nudgeRepaint();
  }
  var cnt=$("#vcount");if(cnt)cnt.textContent=vshow?VAR.filter(vvis).length+" shown of "+VAR.length:VAR.length+" in this locus";
}

/* ---------- variant controls + tooltip ---------- */
var vs=$("#vshow"),vm=$("#vmin");
if(vs)vs.addEventListener("change",function(){vshow=vs.checked;$(".varbar").classList.toggle("off",!vshow);paintVariants();updateBar()});
if(vm)vm.addEventListener("change",function(){vmin=+vm.value;paintVariants();updateBar()});
var tip=$("#tip"),raf=0;
document.addEventListener("mousemove",function(e){
  if(raf)return;raf=requestAnimationFrame(function(){raf=0;
    var p=posAt(e.clientX,e.clientY),html="";
    if(p&&VAR.length){var ur=uRangeOf(p),u=ur[0];
      var vv=varsIn(ur[0],ur[1]);
      if(vv.length)html=vv.map(function(v){var b=BINS[vbin(v.af)];
        return "<span class=\"sq\" style=\"background:"+b[1]+"\"></span><b>"+esc(v.label)+"</b>"+(v.id?" &middot; "+esc(v.id):"")+
        "<br>AF <span class=\"af\">"+afTxt(v.af)+"</span> ("+b[3]+") &middot; "+cdotU(v.u0)+" &middot; "+locU(v.u0)+
        (pFromCs(cFromU(v.u0,v.u1)).length?" &middot; aa "+pFromCs(cFromU(v.u0,v.u1))[0][0]:"")}).join("<hr style=\"border:0;border-top:1px solid #556;margin:4px 0\">")}
    if(html){tip.innerHTML=html;tip.style.display="block";
      var x=e.clientX+14,y=e.clientY+16;if(x+370>innerWidth)x=e.clientX-370;tip.style.left=x+"px";tip.style.top=y+"px"}
    else tip.style.display="none";
  })});

/* ---------- mouse -> position ---------- */
function caret(x,y){
  if(document.caretPositionFromPoint){var p=document.caretPositionFromPoint(x,y);return p&&{node:p.offsetNode,offset:p.offset}}
  if(document.caretRangeFromPoint){var r=document.caretRangeFromPoint(x,y);return r&&{node:r.startContainer,offset:r.startOffset}}
  return null}
function posAt(x,y){
  var c=caret(x,y);if(!c||c.node.nodeType!==3)return null;
  var cont=c.node.parentElement&&c.node.parentElement.closest(".bases");if(!cont)return null;
  var d=c.node.data,off=c.offset;
  var hit=function(o){if(o<0||o>=d.length||!isBase(d[o]))return false;
    var r=document.createRange();r.setStart(c.node,o);r.setEnd(c.node,o+1);
    var b=r.getBoundingClientRect();return x>=b.left-0.5&&x<=b.right+0.5&&y>=b.top-1&&y<=b.bottom+1};
  if(!hit(off)){if(hit(off-1))off--;else return null}
  var a=letters(cont);for(var k=0;k<a.length;k++)if(a[k][0]===c.node&&a[k][1]===off)
    return{space:cont.dataset.space,pos:+cont.dataset.p+k};
  return null}
/* number of bases before a boundary point (text node or element node + offset) */
function lettersBefore(cont,node,off){var a=letters(cont),r=document.createRange();
  try{r.setStart(node,off);r.collapse(true)}catch(e){return a.length}
  for(var k=0;k<a.length;k++)if(r.comparePoint(a[k][0],a[k][1])>=0)return k;
  return a.length}
/* base at a selection boundary: the first base at/after a start, the last base
   before an end, only in boxes of `space` if given. The boundary can be an element
   (a cell, the table, the whole section) when the mouse is let go past the end of a
   line, so boxes are compared with the boundary point itself, not its node. */
function boundary(node,off,isStart,view,space){
  var el=node.nodeType===3?node.parentElement:node,cont=el&&el.closest(".bases");
  if(cont&&(!space||cont.dataset.space===space)){var k=lettersBefore(cont,node,off)-(isStart?0:1);k=Math.max(0,Math.min(k,+cont.dataset.n-1));
    return{space:cont.dataset.space,pos:+cont.dataset.p+k}}
  var P=document.createRange();try{P.setStart(node,off);P.collapse(true)}catch(e){return null}
  var list=$$(".bases",view).filter(function(c){return !space||c.dataset.space===space}),i,c,k;
  if(isStart){for(i=0;i<list.length;i++){c=list[i];
      if(P.comparePoint(c,c.childNodes.length)<0)continue;          // box ends before the point
      k=lettersBefore(c,node,off);if(k<+c.dataset.n)return{space:c.dataset.space,pos:+c.dataset.p+k}}}
  else{for(i=list.length-1;i>=0;i--){c=list[i];
      if(P.comparePoint(c,0)>0)continue;                            // box starts after the point
      k=lettersBefore(c,node,off);if(k>0)return{space:c.dataset.space,pos:+c.dataset.p+k-1}}}
  return null}
/* base on the text line at height y next to x: the last one at or left of x
   (leftOf) or the first one at or right of x; if x is past the line's end, the
   line's last base (leftOf) / the next line's first; null if y is on no line */
function baseNear(x,y,view,space,leftOf){var hit=null,edge=null;
  $$(".bases",view).forEach(function(c){if(c.dataset.space!==space)return;
    var cb=c.getBoundingClientRect();if(y<cb.top-1||y>cb.bottom+1)return;
    var a=letters(c),on=false,k,b,r,P=function(k){return{space:space,pos:+c.dataset.p+k}};
    for(k=0;k<a.length;k++){r=document.createRange();r.setStart(a[k][0],a[k][1]);r.setEnd(a[k][0],a[k][1]+1);b=r.getBoundingClientRect();
      if(y<b.top-1||y>b.bottom+1){if(on){if(!leftOf&&!hit&&!edge)edge=P(k);break}continue}
      if(!on){on=true;if(leftOf&&!hit&&!edge&&k>0)edge=P(k-1)}
      if(leftOf?b.left<=x:(!hit&&b.right>=x))hit=P(k)}});
  return hit||edge}

/* ---------- selection state ---------- */
var anchor=null,sel=null,selVar=null;     // selVar: "highlighted variant(s)" - set by clicking a coloured base
var markPos=null;                          // "marked base" - set with the Mark button (u coordinate), any base
/* gnomAD links: a variant page, or the region around one position.
   hg19 / GRCh37 pages open gnomAD v2.1.1 (dataset=gnomad_r2_1); any other build
   gets no dataset parameter, so gnomAD opens its default (GRCh38) dataset. */
var GNOMAD="https://gnomad.broadinstitute.org/";
var G37=/hg19|grch37/i.test(D.build||"");
var GDS=G37?"?dataset=gnomad_r2_1":"", GNAME=G37?"gnomAD v2.1.1 (GRCh37)":"gnomAD";
function gChr(){return String(D.chr).replace(/^chr/i,"")}
function gnomadPos(u){var g=u2g(u);
  return "<a href=\""+GNOMAD+"region/"+gChr()+"-"+Math.max(1,g-20)+"-"+(g+20)+GDS+"\" target=\"_blank\" rel=\"noopener\""+
    " title=\"Open this position in "+GNAME+", \u00b120 bp\">"+esc(D.chr+":"+fmt(g))+" \u2197</a>"}
function gnomadVar(v){var id=String(v.label).replace(/^chr/i,"");
  return "<a href=\""+GNOMAD+"variant/"+encodeURIComponent(id)+GDS+"\" target=\"_blank\" rel=\"noopener\""+
    " title=\"Open this variant in "+GNAME+"\">"+esc(v.label)+" \u2197</a>"}
function markIv(){return markPos===null?null:varIv([{u0:markPos,u1:markPos}])}
function region(space,a,b){var s={space:space,a:Math.min(a,b),b:Math.max(a,b)};
  s.len=s.b-s.a+1;
  if(space==="c"){s.seq=D.cdna.slice(s.a-1,s.b);var ua=c2u(s.a),ub=c2u(s.b);s.ua=ua;s.ub=ub;
    s.junction=(ub-ua)!==(s.b-s.a)}
  else if(space==="p"){s.a=Math.min(s.a,NP);s.b=Math.min(s.b,NP);s.len=s.b-s.a+1;
    s.seq=D.protein.slice(s.a-1,s.b);var c=cOfP(s.a,s.b);
    s.ua=c2u(c[0]);s.ub=c2u(c[1]);s.ca=c[0];s.cb=c[1];s.junction=false}
  else{s.seq=useq(s.a,s.b);s.ua=s.a;s.ub=s.b;s.junction=false}
  return s}
function describe(s){
  var where=locU(s.ua)===locU(s.ub)?locU(s.ua):locU(s.ua)+" → "+locU(s.ub);
  if(s.space==="p"){var aa=D.protein;
    return "aa "+(s.len===1?fmt(s.a):fmt(s.a)+"-"+fmt(s.b))+" (p."+aa[s.a-1]+s.a+(s.len>1?"_"+aa[s.b-1]+s.b:"")+") &middot; "+
      cdotC(s.ca)+"_"+cdotC(s.cb).replace(/^[cn]\./,"")+" &middot; "+where}
  var h=s.len===1?cdotU(s.ua):cdotU(s.ua)+"_"+cdotU(s.ub).replace(/^[cn]\./,"");
  var cd=s.space==="c"?"cDNA "+(s.len===1?fmt(s.a):fmt(s.a)+"-"+fmt(s.b))+" &middot; ":"";
  return cd+h+" &middot; "+where+" &middot; "+gtxt(s.ua,s.ub)+
    (s.junction?" <span class=\"hint\">(spans exon junction)</span>":"")}
function varIv(list){var o={u:[],c:[],p:[]};list.forEach(function(v){var cs=cFromU(v.u0,v.u1);
  o.u.push([v.u0,v.u1]);cs.forEach(function(x){o.c.push(x)});pFromCs(cs).forEach(function(x){o.p.push(x)})});return o}
function varTxt(v){return v.label+(v.id?" ("+v.id+")":"")+" AF "+afTxt(v.af)+" · "+cdotU(v.u0)}
/* variants (AF >= current filter) inside a selection; all=true also counts
   them while the variant checkbox is off. Uses the selection's own pieces,
   so a cDNA selection across a junction does not count intronic variants. */
function varsInSel(s,all){var iv=ivOf(s).u;return VAR.filter(function(v){if(v.af<vmin||(!all&&!vshow))return false;
  for(var i=0;i<iv.length;i++)if(v.u1>=iv[i][0]&&v.u0<=iv[i][1])return true;return false})}
function pctTxt(x){return String(parseFloat((x*100).toPrecision(3)))+"%"}
var VLABEL=$("#vlabel")?$("#vlabel").textContent.trim():"variants";
/* one-line header for "Add context" */
function ctxHeader(s,isRC){
  var kind=s.space==="c"?"cDNA":s.space==="p"?"protein":"genomic (unspliced)";
  if(D.strand==="-"&&s.space!=="p")kind+=", transcript orientation";
  if(isRC)kind+=", reverse complement";
  var pos=s.space==="p"?"p."+D.protein[s.a-1]+s.a+(s.len>1?"_"+D.protein[s.b-1]+s.b:""):
    (s.len===1?cdotU(s.ua):cdotU(s.ua)+"_"+cdotU(s.ub).replace(/^[cn]\./,""));
  var where=locU(s.ua)===locU(s.ub)?locU(s.ua):locU(s.ua)+" to "+locU(s.ub);
  return ">"+[D.gene_name,D.transcript_id+(D.tname?" ("+D.tname+")":"")+(D.mane?" (not MANE; MANE is "+D.mane+")":""),pos,where,gtxt(s.ua,s.ub),D.build,kind,
    fmt(s.len)+(s.space==="p"?" aa":" bp")].filter(function(x){return x}).join(" | ")}

function paintSel(){
  if(!HL)return;
  var hs=new Highlight(),ha=new Highlight(),hv=new Highlight(),hc=new Highlight();
  if(sel)paint(sel.len===1?ha:hs,ivOf(sel));
  // highlighted variant: blue on yellow, on screen only (the RTF keeps its gnomAD colour)
  if(selVar&&vshow)paint(hc,varIv(selVar));
  // marked base: red on yellow, drawn above the highlighted variant where they overlap
  if(markPos!==null)paint(hv,markIv());
  setHL("sel",hs,1);setHL("anc",ha,1);setHL("vchk",hc,5);setHL("vsel",hv,6);
  nudgeRepaint();
}

/* ---------- bottom drawer (tab: click = open/close, drag = resize) ---------- */
var bar=$("#pickbar"),dbody=$(".dbody",bar),dtab=$(".dtab",bar);
var store={get:function(k){try{return localStorage.getItem(k)}catch(e){return null}},
           set:function(k,v){try{localStorage.setItem(k,v)}catch(e){}}};
var dOpen=store.get("mp_open")!=="0",dH=+store.get("mp_h")||0;
function padMain(){requestAnimationFrame(function(){
  $("main").style.paddingBottom=(getComputedStyle(bar).display==="none"?40:bar.offsetHeight+24)+"px"})}
function applyDrawer(){
  bar.classList.toggle("collapsed",!dOpen);
  var t=$(".dtoggle",bar);t.textContent=dOpen?"▼":"▲";t.setAttribute("aria-expanded",dOpen);
  t.title=dOpen?"Collapse":"Expand";
  dbody.style.height=dH?dH+"px":"";padMain()}
function openDrawer(o){dOpen=o;store.set("mp_open",o?"1":"0");applyDrawer()}
dtab.addEventListener("pointerdown",function(e){
  if(e.button!==0)return;e.preventDefault();
  var y0=e.clientY,h0=dOpen?dbody.offsetHeight:0,moved=false;dtab.setPointerCapture(e.pointerId);
  function mv(ev){var dy=y0-ev.clientY;if(!moved&&Math.abs(dy)<5)return;moved=true;
    if(!dOpen){dOpen=true;bar.classList.remove("collapsed")}
    dH=Math.round(Math.max(70,Math.min(innerHeight*0.8,h0+dy)));dbody.style.height=dH+"px";padMain()}
  function up(){dtab.removeEventListener("pointermove",mv);dtab.removeEventListener("pointerup",up);
    if(moved){store.set("mp_h",String(dH));openDrawer(true)}else openDrawer(!dOpen)}
  dtab.addEventListener("pointermove",mv);dtab.addEventListener("pointerup",up)});
window.addEventListener("resize",padMain);

function copyText(text,btn){var done=function(){var o=btn.textContent;btn.textContent="Copied";setTimeout(function(){btn.textContent=o},1100)};
  if(navigator.clipboard&&window.isSecureContext)navigator.clipboard.writeText(text).then(done,function(){fallback()});else fallback();
  function fallback(){var t=document.createElement("textarea");t.value=text;document.body.appendChild(t);t.select();document.execCommand("copy");t.remove();done()}}

function updateBar(){
  document.body.classList.toggle("has-pick",!!sel||!!(selVar&&vshow)||markPos!==null);
  var mb=$("button.mark",bar),one=sel&&sel.len===1&&sel.space!=="p"&&sel.ua!==null;
  mb.hidden=!one&&markPos===null;
  mb.textContent=one&&markPos!==sel.ua?"\u2605 Mark this base":"\u2605 Unmark";
  var sp=$(".pane.seq",bar),rc=$("input.rc",sp);
  if(sel){var s=sel.seq;
    if(s===null){$("textarea",sp).value="";$(".info",sp).innerHTML="<b>Selection includes intron sequence that is not embedded in this page</b> (re-run with a larger intron_max)"}
    else{var isRC=rc.checked&&sel.space!=="p";if(isRC)s=revcomp(s);
      $("textarea",sp).value=($("input.ctx",sp).checked?ctxHeader(sel,isRC)+"\n":"")+s;
      var vv=varsInSel(sel,false);
      var unit=sel.space==="p"?" aa":" bp";
      $(".info",sp).innerHTML=(sel.len===1?"<b>"+s+"</b> &nbsp;":"<b>"+fmt(sel.len)+unit+"</b> &nbsp;")+describe(sel)+
        (rc.checked&&sel.space!=="p"?" &middot; reverse complement":"")+(vv.length?" &middot; <span style=\"color:#c92a2a\">"+vv.length+" variant"+(vv.length>1?"s":"")+"</span>":"")}}
  // warnings: hidden variants, reference differences, strand
  var W=[];
  if(sel&&sel.seq!==null){
    if(!vshow&&VAR.length){var hv=varsInSel(sel,true);
      if(hv.length)W.push("⚠ "+hv.length+" variant"+(hv.length>1?"s":"")+" ≥ "+pctTxt(vmin)+" in this selection (hidden – tick "+esc(VLABEL)+" to see them)")}
    if(D.refdiff.length){var ci=ivOf(sel).c,hit=D.refdiff.filter(function(c){return ci.some(function(x){return c>=x[0]&&c<=x[1]})});
      if(hit.length)W.push("⚠ includes "+hit.length+" base"+(hit.length>1?"s":"")+(D.mane?" where this transcript differs from MANE "+esc(D.mane):
        " where the "+esc(D.build||"genome")+" reference differs from the transcript")+" ("+hit.slice(0,5).map(cdotC).join(", ")+")")}
    if(D.lendiff)W.push("⚠ cDNA length differs from "+(D.mane?"MANE "+esc(D.mane):"the transcript")+": c. numbers and amino acids here may be shifted – use the genomic position");
    if(D.strand==="-"&&sel.space!=="p")W.push("<span class=\"n\">Gene on the reverse strand: sequence is 5′→3′ of the transcript (reverse complement of genomic +)</span>")}
  $(".warns",sp).innerHTML=W.map(function(x){return "<span>"+x+"</span>"}).join("");
  var mv=$(".markvar",sp);
  var ML=[];
  if(markPos!==null)ML.push("<b>\u2605 Marked base:</b> "+esc(cdotU(markPos))+" &middot; "+esc(locU(markPos))+" &middot; "+gnomadPos(markPos)+
    " <span class=\"hint\">(red in every tab and in the RTF)</span>");
  if(selVar&&vshow)ML.push("<span class=\"mvb\"><b>Highlighted variant"+(selVar.length>1?"s":"")+":</b> "+selVar.map(function(v){
    return gnomadVar(v)+(v.id?" ("+esc(v.id)+")":"")+" AF "+afTxt(v.af)+" &middot; "+esc(cdotU(v.u0))}).join(" ; ")+"</span>"+
    " <span class=\"hint\">(blue on screen only, for checking; the RTF keeps its gnomAD colour)</span>");
  mv.innerHTML=ML.join("<br>");
  // summary shown on the drawer tab (useful when collapsed)
  $(".dsum",bar).textContent=(sel?(sel.len===1?(sel.space==="p"?"1 aa":"1 base"):fmt(sel.len)+(sel.space==="p"?" aa":" bp"))+" · "+cdotU(sel.ua):"")+
    (markPos!==null?"  \u2605 "+cdotU(markPos):selVar&&vshow?"  variant "+selVar[0].label:"");
  paintSel();padMain();
}
$("button.copy",bar).addEventListener("click",function(){copyText($(".pane.seq textarea",bar).value,this)});
$("input.rc",bar).addEventListener("change",updateBar);
$("input.ctx",bar).addEventListener("change",function(){store.set("mp_ctx",this.checked?"1":"0");updateBar()});
$("button.clear",bar).addEventListener("click",function(){anchor=null;sel=null;selVar=null;markPos=null;updateBar()});
$("button.mark",bar).addEventListener("click",function(){
  if(sel&&sel.len===1&&sel.space!=="p"&&sel.ua!==null&&markPos!==sel.ua){
    markPos=sel.ua;var vv=varsIn(markPos,markPos);selVar=vv.length?vv:null}
  else{markPos=null;selVar=null}
  updateBar()});

/* function commit(space,a,b){sel=region(space,a,b);updateBar()} */
/* function commit(space,a,b){if(currentView()==="protein")return;sel=region(space,a,b);updateBar()} */
function commit(space,a,b){if(space==="p")return;sel=region(space,a,b);updateBar()}

/* ---------- mouse handling ---------- */
document.addEventListener("mousedown",function(e){
  if(e.shiftKey&&anchor&&e.target.closest&&e.target.closest(".bases"))e.preventDefault()});
document.addEventListener("mouseup",function(e){
  if(e.button!==0)return;
  // drag: judged by where the selection started, so letting go past the end of a
  // line, in the margin or over the drawer still counts
  var ns=window.getSelection(),an=ns&&ns.rangeCount&&!ns.isCollapsed?ns.anchorNode:null;
  var dv=an&&(an.nodeType===3?an.parentElement:an);dv=dv&&dv.closest&&dv.closest("section.view");
  if(dv){
    var r=ns.getRangeAt(0),a=boundary(r.startContainer,r.startOffset,true,dv),b=a&&boundary(r.endContainer,r.endOffset,false,dv,a.space);
    // let go outside the sequence (past a line's end, in the margin): browsers may stretch
    // the selection to the whole section, so take the end the mouse moved from the mouse
    if(a){var inSeq=function(n){var el=n.nodeType===3?n.parentElement:n;el=el&&el.closest(".bases");return el&&el.dataset.space===a.space},
        fwd=ns.focusNode===r.endContainer&&ns.focusOffset===r.endOffset,m;
      if(fwd&&!inSeq(r.endContainer)){m=baseNear(e.clientX,e.clientY,dv,a.space,true);if(m&&m.pos>=a.pos)b=m}
      else if(!fwd&&!inSeq(r.startContainer)){m=baseNear(e.clientX,e.clientY,dv,a.space,false);if(m&&b&&m.pos<=b.pos)a=m}}
    // never leave a browser selection over the sequence: pressing inside one starts
    // a text drag instead of a new selection, so the next try would do nothing
    if($$(".bases",dv).some(function(c){return r.intersectsNode(c)}))ns.removeAllRanges();
    if(a&&b&&b.pos>=a.pos){anchor={space:a.space,pos:a.pos};commit(a.space,a.pos,b.pos)}
    return}
  if(bar.contains(e.target))return;
  var view=e.target.closest&&e.target.closest("section.view");if(!view)return;
  var p=posAt(e.clientX,e.clientY);if(!p)return;
  if(e.shiftKey&&anchor&&anchor.space===p.space){commit(p.space,anchor.pos,p.pos)}
  else{
    // a single click on a coloured (visible) variant highlights it in every tab;
    // this is only for spot-checking, so the marked base stays marked
    var ur=uRangeOf(p),vv=varsIn(ur[0],ur[1]);if(vv.length)selVar=vv;
    anchor=p;commit(p.space,p.pos,p.pos)}
});
document.addEventListener("keydown",function(e){if(e.key!=="Escape")return;anchor=null;sel=null;selVar=null;markPos=null;updateBar()});

/* ---------- Go to: c.123, c.68-7, c.-20, c.*50, exon 3, intron 2, chr:pos, p.45 ---------- */
function parseGo(t){
  t=t.trim().replace(/,/g,"").replace(/\s+/g," ");if(!t)return null;var m,c,u;
  if((m=/^(?:p\.|aa ?)([A-Za-z]{0,3})(\d+)$/i.exec(t))){var k=+m[2];
    if(!NAA||k<1||k>NAA)return{err:NAA?"the protein has "+NAA+" aa":"non-coding transcript"};
    c=D.cds_s+3*(k-1);return{p:k,c:c,u:c2u(c),pref:"protein"}}
  if((m=/^ex(?:on)? ?(\d+)$/i.exec(t))){var e=EX[+m[1]-1];if(!e)return{err:"there are "+EX.length+" exons"};return{u:e.us,c:e.cs}}
  if((m=/^in(?:tron)? ?(\d+)$/i.exec(t))){var it=IN[+m[1]-1];if(!it)return{err:"there are "+IN.length+" introns"};return{u:it.us,c:null,pref:"exons"}}
  if((m=/^([cn])\.([-*]?)(\d+)(?:([+-])(\d+))?$/i.exec(t))){var N=+m[3],cp;
    if(D.cds_s===null||m[1].toLowerCase()==="n")cp=m[2]==="-"?1-N:m[2]==="*"?NC+N:N;else cp=m[2]==="-"?D.cds_s-N:m[2]==="*"?D.cds_e+N:D.cds_s+N-1;
    var inx=cp>=1&&cp<=NC,base=inx?c2u(cp):(cp<1?cp:cp-NC+L);
    if(m[4]){if(!inx)return{err:"+/- offsets need a position inside the transcript"};
      u=base+(m[4]==="+"?+m[5]:-m[5]);return{u:u,c:u2c(u),pref:u2c(u)===null?"exons":null}}
    return{u:base,c:inx?cp:null}}
  var g=null;
  if((m=/^(?:chr)?([0-9XYMT]+):(\d+)$/i.exec(t))){
    var mt=function(x){x=String(x).toUpperCase().replace(/^CHR/,"");return x==="M"?"MT":x};
    if(mt(m[1])!==mt(D.chr))return{err:"this gene is on chromosome "+D.chr};g=+m[2]}
  else if(/^\d+$/.test(t)){var n=+t;
    if(n>=D.tx_s-D.embed&&n<=D.tx_e+D.embed)g=n;else return parseGo("c."+n)}
  if(g!==null){u=D.strand==="-"?D.tx_e-g+1:g-D.tx_s+1;return{u:u,c:u2c(u)}}
  return{err:"not recognised – try c.123, c.68-7, exon 3, "+D.chr+":"+D.tx_s+" or p.45"}}
function goTo(t,msg){
  var r=parseGo(t);msg.textContent="";if(!r)return;if(r.err){msg.textContent=r.err;return}
  if(r.u!=null&&(r.u<1-D.up.length||r.u>L+D.dn.length)){msg.textContent="outside the region on this page";return}
  var cur=currentView(),view,space,pos;
  if(r.pref==="protein"&&$("#protein")){view="protein";space="p";pos=r.p}
  else if(cur==="protein"&&r.pref!=="exons"&&r.c!=null&&pFromC(r.c,r.c).length&&pFromC(r.c,r.c)[0][0]<=NAA){view="protein";space="p";pos=pFromC(r.c,r.c)[0][0]}
  else if(cur==="cdna"&&r.pref!=="exons"&&r.c!=null){view="cdna";space="c";pos=r.c}
  else if(cur==="exons"||r.c==null||r.pref==="exons"){view="exons";space="g";pos=r.u}
  else{view="cdna";space="c";pos=r.c}
  if(view!==cur)clearVariantColours();
  show(view);history.replaceState(null,"","#"+view);
  if(view==="exons"){
    // widen the intron/flank dropdown if the position is not shown yet
    var need=0,ilen=0;
    if(pos<1)need=1-pos;else if(pos>L)need=pos-L;
    else if(exonU(pos)<0)for(var i=0;i<IN.length;i++)if(pos>=IN[i].us&&pos<=IN[i].ue){
      var d=Math.min(pos-IN[i].us+1,IN[i].ue-pos+1);ilen=IN[i].len;
      if(ilen>2*D.embed&&d>D.embed){msg.textContent=fmt(d)+" bp into the intron – beyond the "+fmt(D.embed)+" bp included on this page";return}
      need=d;break}
    if(need>0&&isel){var shown=function(v){return (pos<1||pos>L)?Math.max(D.flank,v)>=need:(v>=need||ilen<=2*v)};
      if(!shown(+isel.value)){var opt=[].slice.call(isel.options).map(function(o){return +o.value}).filter(shown);
        if(!opt.length){msg.textContent="beyond the sequence included on this page";return}
        isel.value=opt[0];renderFlex(opt[0])}}
    // switch back on rows hidden by the markup boxes
    var tg=(pos<1||pos>L)?"no-flanks":(exonU(pos)<0?"no-introns":null),cb=tg&&$('#exons input[data-toggle="'+tg+'"]');
    if(cb&&!cb.checked){cb.checked=true;cb.dispatchEvent(new Event("change"))}
  }
  var cont=$$("#"+view+" .bases").filter(function(c){return c.dataset.space===space&&pos>=+c.dataset.p&&pos<+c.dataset.p+ +c.dataset.n})[0];
  if(!cont){msg.textContent="not shown in this tab";return}
  var a=letters(cont)[pos-+cont.dataset.p],rg=document.createRange();rg.setStart(a[0],a[1]);rg.setEnd(a[0],a[1]+1);
  window.scrollBy({top:rg.getBoundingClientRect().top-innerHeight*0.35,behavior:"smooth"});
  anchor={space:space,pos:pos};commit(space,pos,pos);
}
$$("form.goto").forEach(function(f){var inp=$("input",f),msg=$(".gmsg",f);
  f.addEventListener("submit",function(e){e.preventDefault();goTo(inp.value,msg)});
  inp.addEventListener("input",function(){msg.textContent=""})});

/* ---------- RTF export of the current tab, with the colours shown ----------
   Text colours/backgrounds are read from the page as displayed (exons, UTR,
   codons, translation, intron/flank, splice residues), hidden rows/lines are
   left out, and gnomAD colours + the marked base are added on top. */
function rtfExport(){
  var view=currentView(),root=$("section.view#"+view),styleC=new Map(),hidC=new Map();
  function hex(c){var m=/rgba?\(([^)]+)\)/.exec(c||"");if(!m)return null;var p=m[1].split(",").map(parseFloat);
    if(p.length>3&&p[3]===0)return null;return "#"+p.slice(0,3).map(function(v){return ("0"+Math.round(v).toString(16)).slice(-2)}).join("")}
  function hidden(el){if(!el||el===root)return false;if(hidC.has(el))return hidC.get(el);
    var h=getComputedStyle(el).display==="none"||hidden(el.parentElement);hidC.set(el,h);return h}
  function baseStyle(el){var s=styleC.get(el);if(s)return s;var cs=getComputedStyle(el),bg=null,e=el;
    while(e&&e.tagName!=="PRE"&&e.tagName!=="TD"){bg=hex(getComputedStyle(e).backgroundColor);if(bg)break;e=e.parentElement}
    s={fg:hex(cs.color)||"#000000",bg:bg};styleC.set(el,s);return s}
  // overlay lookups: key = space + position
  var vmap={},smap={};
  if(vshow){VAR.forEach(function(v){if(!vvis(v))return;var b=vbin(v.af),iv=varIv([v]);
    ["u","c","p"].forEach(function(k){iv[k].forEach(function(x){for(var i=x[0];i<=x[1];i++){var key=k+i;if(!(key in vmap)||vmap[key]>b)vmap[key]=b}})})})}
  [markIv()].forEach(function(iv){if(!iv)return;
    ["u","c","p"].forEach(function(k){iv[k].forEach(function(x){for(var i=x[0];i<=x[1];i++)smap[k+i]=1})})});
  function overlay(space,pos,st){var key=spaceKey(space)+pos;
    if(smap[key])return{fg:"#e00000",bg:st.bg,size:40,b:1};
    if(key in vmap){var B=BINS[vmap[key]];return{fg:B[2],bg:B[1],uld:1}}
    return st}
  // walk a <pre>: returns lines of [char, style]
  function collect(pre,lines){var cur=[],cnt=new Map(),w=document.createTreeWalker(pre,NodeFilter.SHOW_TEXT),n;
    while((n=w.nextNode())){var el=n.parentElement;if(hidden(el))continue;
      var st=baseStyle(el),cont=el.closest(".bases");
      for(var i=0;i<n.data.length;i++){var ch=n.data[i];
        if(ch==="\n"){lines.push(cur);cur=[];continue}
        var s2=st;if(cont&&isBase(ch)){var k=cnt.get(cont)||0;cnt.set(cont,k+1);s2=overlay(cont.dataset.space,+cont.dataset.p+k,st)}
        cur.push([ch,s2])}}
    if(cur.length)lines.push(cur);return lines}
  var T=function(t,st){return t.split("").map(function(ch){return[ch,st]})};
  var HEAD={fg:"#24435f",f:1,size:24,b:1},SUB={fg:"#5d6b7a",f:1,size:18};
  var title={cdna:"cDNA sequence",exons:"Exons",protein:"Protein sequence"}[view];
  var LINES=[];
  LINES.push(T(D.gene_name+" \u2013 "+D.transcript_id+" \u2013 "+title,HEAD));
  LINES.push(T($(".summary dd:nth-of-type(2)").textContent.trim(),SUB));
  if(vshow&&BINS.length){var lg=T((($("#vlabel")&&$("#vlabel").textContent.trim())||"Variants")+" variants, AF bands: ",SUB);
    BINS.forEach(function(B){lg=lg.concat(T(" "+B[3].replace(/&ge;/g,"\u2265").replace(/&ndash;/g,"\u2013").replace(/&lt;/g,"<")+" ",{fg:B[2],bg:B[1],f:1,size:18,uld:1}),T("  ",SUB))});
    LINES.push(lg)}
  if(markPos!==null)LINES.push(T("Marked base (red, 20 pt): ",{fg:"#e00000",f:1,size:18,b:1}).concat(
    T(cdotU(markPos)+" \u00b7 "+locU(markPos)+" \u00b7 "+D.chr+":"+fmt(u2g(markPos))+" ("+(D.build||"genome")+")",SUB)));
  if(selVar&&vshow)LINES.push(T("Highlighted variant"+(selVar.length>1?"s":"")+": ",{fg:"#5d6b7a",f:1,size:18,b:1}).concat(T(selVar.map(varTxt).join(" ; "),SUB)));
  if($(".verwarn"))LINES.push(T($(".verwarn").textContent.trim(),{fg:"#6b4e00",f:1,size:18}));
  if($(".refwarn"))LINES.push(T($(".refwarn").textContent.trim(),{fg:"#8a1c1c",f:1,size:18}));
  if(D.prov)LINES.push(T(D.prov,SUB));
  LINES.push([]);
  if(view==="exons"){
    // sequence blocks only (like Ensembl's export), separated by a blank line
    $$("table.exons tbody tr",root).forEach(function(tr){if(hidden(tr))return;
      var td=tr.children,pre=$("pre",td[td.length-1]);if(pre){collect(pre,LINES);LINES.push([])}});
  } else {var pre=$("pre",root);if(pre)collect(pre,LINES)}
  // ---- write RTF
  var cols=[],ci={};function cidx(h){if(!(h in ci)){cols.push(h);ci[h]=cols.length}return ci[h]}
  function esc(t){return t.replace(/[\\{}]/g,function(c){return "\\"+c}).replace(/[^\x20-\x7e]/g,function(c){var k=c.charCodeAt(0);return "\\u"+(k>32767?k-65536:k)+"?"})}
  function sig(s){return [s.fg,s.bg,s.size,s.b,s.uld,s.f].join("|")}
  var body=LINES.map(function(line){
    if(!line.length)return "{\\pard\\fs18\\f0\\par}";
    var out="",run="",cur=null;
    function flush(){if(!run)return;var s=cur,ctl="";
      if(s.f)ctl+="\\f"+s.f;if(s.fg)ctl+="\\cf"+cidx(s.fg);
      if(s.bg){var j=cidx(s.bg);ctl+="\\chshdng0\\chcbpat"+j+"\\cb"+j}
      if(s.size)ctl+="\\fs"+s.size;if(s.b)ctl+="\\b";if(s.uld)ctl+="\\uld";
      out+="{"+ctl+" "+esc(run)+"}";run=""}
    line.forEach(function(x){if(!cur||sig(x[1])!==sig(cur)){flush();cur=x[1]}run+=x[0]});flush();
    return "{\\pard\\fs18\\f0 "+out+"\\par}"}).join("\n");
  var ct=cols.map(function(h){return "\\red"+parseInt(h.substr(1,2),16)+"\\green"+parseInt(h.substr(3,2),16)+"\\blue"+parseInt(h.substr(5,2),16)+";"}).join("");
  var doc="{\\rtf1\\ansi\\deff0{\\fonttbl{\\f0\\fmodern Courier New;}{\\f1\\fswiss Arial;}}\n{\\colortbl ;"+ct+"}\n"+
    "\\paperw16838\\paperh11906\\margl720\\margr720\\margt720\\margb720\\landscape\n"+body+"\n}";
  download(D.gene_name+"_"+D.transcript_id+"_"+view+".rtf",doc,"application/rtf");
  return doc}
$$("button.rtf").forEach(function(b){b.addEventListener("click",rtfExport)});
window.__rtfExport=rtfExport;   // used by the test suite

/* ---------- start ---------- */
$("input.ctx",bar).checked=store.get("mp_ctx")==="1";
// browsers may restore a ticked variant box on reload / Back: follow it
if(vs){vshow=vs.checked;$(".varbar").classList.toggle("off",!vshow)}
applyDrawer();updateBar();
renderFlex(isel?+isel.value:25);restripe();
})();
)---"

page_html <- function(m, vars, width = 60, intron_ends = 25, intron_opts = c(25, 50, 100, 250, 500, 1000),
                      index_href = "index.html", bins = .AF_BINS, variant_label = "lab-specific",
                      af_default = 0.001, build = NULL, prov = NULL, var_source = NULL, refdiff = NULL,
                      ref_checked = FALSE, var_nochr = FALSE) {
  intron_ends <- min(intron_ends, m$embed)    # more than is embedded would break the Exons view
  strand_txt <- if (m$strand == "-") "reverse strand" else "forward strand"
  loc <- sprintf("Chromosome %s: %s-%s %s.", m$chr, fmt_int(m$tx_start), fmt_int(m$tx_end), strand_txt)
  n_aa <- if (m$has_cds) nchar(sub("\\*$", "", m$protein)) else 0
  cds_len <- if (m$has_cds) m$cds_e - m$cds_s + 1 else 0
  chip <- function(label, value) sprintf('<span class="chip"><span class="k">%s</span><span class="v">%s</span></span>', label, value)
  bi <- build_info(build)
  nm <- if (!is.na(m$transcript_name) && grepl("^[NX][MR]_[0-9]+", m$transcript_name)) m$transcript_name else NA
  about <- paste0(
    if (!is.na(nm)) paste0(ext_link(paste0("https://www.ncbi.nlm.nih.gov/nuccore/", nm), paste0("<b>", html_escape(nm), "</b>"),
                                    "Open this RefSeq transcript at NCBI"), ' <span class="sep">&middot;</span> ') else "",
    sprintf("This transcript has %d exon%s ", nrow(m$exons), if (nrow(m$exons) == 1) "" else "s"),
    chip("cDNA", paste(fmt_int(nchar(m$cdna)), "bp")),
    if (m$has_cds) chip("CDS", sprintf("%s bp <span class=\"sub\">(cDNA %s&ndash;%s)</span>", fmt_int(cds_len),
                                       fmt_int(m$cds_s), fmt_int(m$cds_e))) else "",
    if (m$has_cds) chip("Protein", paste(fmt_int(n_aa), "aa")) else "",
    chip("Biotype", html_escape(m$biotype %||% "-")),
    if (!is.na(m$protein_id)) chip("Protein ID", html_escape(m$protein_id)) else "",
    if (!is.na(m$ccds)) chip("CCDS", html_escape(m$ccds)) else "")
  nice_tag <- function(t) switch(t, MANE_Select = "MANE Select", MANE_Plus_Clinical = "MANE Plus Clinical",
                                 Ensembl_canonical = "Ensembl Canonical", GENCODE_Primary = "GENCODE Primary",
                                 basic = "GENCODE basic", sub("^appris ", "APPRIS ", gsub("_", " ", t)))
  fb_status <- if ("MANE_Plus_Clinical" %in% m$tags) "MANE Plus Clinical" else "MANE Select"
  if (.is_fb(m)) m$tags <- character(0)
  badges <- paste(vapply(m$tags, function(t) sprintf('<span class="badge%s">%s</span>',
                    if (grepl("^MANE", t)) " mane" else "", html_escape(nice_tag(t))), ""), collapse = "")
  # pre-escaped; the dash is an HTML entity so no non-ASCII string enters paste()
  tname2 <- if (!is.na(m$transcript_name) && m$transcript_name != m$transcript_id) m$transcript_name else NA
  title <- paste(c(html_escape(m$gene_name), html_escape(m$transcript_id),
                   if (!is.na(tname2)) html_escape(tname2)), collapse = " &ndash; ")

  dl <- '<button class="rtf" title="Download this tab as an RTF file with the colours shown">&#8681; Download RTF</button>'
  hint <- '<span class="hint">Click a base &middot; Shift+click or drag for a range &middot; click a coloured variant to highlight it &middot; select a base and press &#9733; Mark to mark it &middot; Esc clears</span>'
  goto <- paste0('<form class="goto" autocomplete="off"><input type="text" spellcheck="false" ',
    'placeholder="Go to: c.5946, c.68-7, exon 11, ', m$chr, ':', m$tx_start, ', p.100" aria-label="Go to position">',
    '<button type="submit">Go</button><span class="gmsg"></span></form>')
  toolbar <- paste0('<div class="toolbar">', goto, dl, hint, '</div>')
  tog <- function(label, cls, target) sprintf(
    '<label><input type="checkbox" checked data-toggle="%s" data-target="%s"> %s</label>', cls, target, label)
  opts <- sort(unique(c(intron_opts[intron_opts <= m$embed], intron_ends)))
  isel <- paste0('<select id="intronsize" aria-label="Intron and flank sequence shown">',
                 paste(sprintf('<option value="%d"%s>%s bp</option>', opts, ifelse(opts == intron_ends, " selected", ""),
                               fmt_int(opts)), collapse = ""), '</select>')
  mins <- bins$min[bins$min > 0]
  pct <- function(x) paste0(format(x * 100, scientific = FALSE, drop0trailing = TRUE, trim = TRUE), "%")
  varbar <- if (isTRUE(attr(vars, "supplied"))) paste0(
    '<div class="varbar off"><b>Variant allele frequencies</b>',
    '<label', if (!is.null(var_source)) paste0(' title="', html_escape(var_source), '"') else "",
    '><input type="checkbox" id="vshow"> <b id="vlabel">', html_escape(variant_label), '</b>',
    if (!is.null(var_source)) ' <span class="info-i">i</span>' else "", '</label>',
    '<label class="vfilter">AF &ge; <select id="vmin">',
    paste(sprintf('<option value="%s"%s>%s</option>', format(rev(mins), scientific = FALSE, drop0trailing = TRUE, trim = TRUE),
                  ifelse(abs(rev(mins) - af_default) < 1e-12, " selected", ""), pct(rev(mins))), collapse = ""),
    '</select></label>',
    '<span class="vlegend">', paste(sprintf('<span class="vsw vb%d">%s</span>', seq_len(nrow(bins)) - 1, bins$label), collapse = ""),
    '</span><span class="hint" id="vcount"></span><span class="hint">Hover a coloured base for details</span>',
    if (var_nochr) paste0('<span class="nochr">&#9888; No variant data for chromosome ', html_escape(sub("^chr", "", m$chr)),
                          ' in this dataset &ndash; absence of colours does not mean absence of variants</span>') else "",
    '</div>')
  else '<div class="varbar"><b>Variant allele frequencies</b><span class="hint">No variant data supplied</span></div>'

  pickbar <- paste0(
'<div id="pickbar" role="region" aria-label="Selected sequence">
<div class="dtab" title="Click to open/close &middot; drag to resize"><button class="dtoggle" aria-label="Open or close" aria-expanded="true">&#9660;</button>
<span class="dtitle">Selection</span><span class="dsum"></span><span class="grip" aria-hidden="true">&#8942;&#8942;</span></div>
<div class="dbody">
<div class="pane seq"><div class="row"><span class="info"></span></div>
<div class="row markvar" style="color:#c92a2a;font-size:13px;margin-top:2px"></div><div class="row" style="margin-top:6px">
<textarea readonly spellcheck="false" aria-label="Selected sequence"></textarea>
<button class="copy primary">Copy</button>
<label><input type="checkbox" class="rc"> Reverse complement</label>
<label title="Adds a header line: gene, transcript, position, build, strand"><input type="checkbox" class="ctx"> Add context</label>
<button class="mark" hidden title="Mark the selected base in red in every tab and in the RTF (works on any base)">&#9733; Mark this base</button>
<button class="clear">Clear</button></div>
<div class="row warns"></div></div>
</div></div><div id="tip"></div>')

  paste0(
'<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>', title, '</title>
<style>', .PAGE_CSS, af_bin_css(bins), '</style></head><body>
<header class="top"><h1>Transcript: ', html_escape(m$transcript_id),
  if (!is.na(tname2)) paste0('<span class="sep">&ndash;</span>', html_escape(tname2)) else "",
  if (.is_fb(m)) ' <span class="badge old">NOT MANE</span>' else "",
'</h1><a href="', index_href, '">&larr; all MANE transcripts</a></header>
<main>
<dl class="summary">
<dt>Gene</dt><dd><b>', html_escape(m$gene_name), '</b><span class="sep">&ndash;</span>',
  if (grepl("^ENSG", m$gene_id)) ext_link(paste0(bi$ens, "/Homo_sapiens/Gene/Summary?g=", sub("\\.[0-9]+$", "", m$gene_id)),
                                          html_escape(m$gene_id), "Open this gene in Ensembl") else html_escape(m$gene_id),
  if (!is.na(m$gene_biotype)) paste0(" &middot; ", html_escape(m$gene_biotype)) else "", '</dd>
<dt>Location', if (!is.null(bi$short)) paste0(" (", bi$short, ")") else "", '</dt><dd>',
  if (!is.null(bi$db)) ext_link(sprintf("https://genome.ucsc.edu/cgi-bin/hgTracks?db=%s&position=%s:%d-%d", bi$db,
                                        if (grepl("^chr", m$chr)) m$chr else if (m$chr %in% c("MT", "M")) "chrM" else paste0("chr", m$chr),
                                        m$tx_start, m$tx_end),
                                html_escape(loc), paste("Open this region in the UCSC Genome Browser,", bi$db))
  else html_escape(loc),
  if (!is.null(build)) paste0(' <span class="build" title="Genome assembly of all coordinates and sequences on this page">',
                              html_escape(build), '</span>') else "",
  if (m$strand == "-") ' <span class="strand-badge" title="Sequences are shown 5&prime;&rarr;3&prime; of the transcript, i.e. the reverse complement of the genomic + strand">reverse strand</span>' else "",
  '</dd>
<dt>About this transcript</dt><dd>', about, '</dd>
<dt>Flags</dt><dd>', if (.is_fb(m)) paste0('<span class="badge old">', tools::toTitleCase(.fb_word(m)), ' version of the ', fb_status,
    ' transcript (', html_escape(m$mane_refseq), ')</span>') else badges, '</dd>
</dl>
', if (.is_fb(m))
  paste0('<div class="verwarn">&#9888; <b>This is not the MANE transcript.</b> The ', fb_status, ' transcript for ',
         html_escape(m$gene_name), ' is <b>', html_escape(m$mane_refseq), '</b> (', html_escape(m$mane_enst),
         '), which is not placed on ', html_escape(build %||% "this genome build"), '. This page shows <b>',
         html_escape(m$transcript_id), '</b>, the ', .fb_word(m), ' version that is placed &ndash; a real RefSeq transcript with its own ',
         'exon structure, sequence and c. numbering. Check variant nomenclature against ', html_escape(m$mane_refseq), ' before reporting. ',
         if (!ref_checked) 'Not compared with the MANE sequence.'
         else if (is.null(refdiff)) 'Compared with the MANE sequence: <b>identical</b>.'
         else paste0('Compared with the MANE sequence: <b>different</b> &ndash; ', html_escape(refdiff$txt), '.',
                     if (isTRUE(refdiff$len_diff)) ' c. numbers and the translation on this page may be shifted against it.' else ""),
         '</div>') else "",
  if (!is.null(refdiff) && !.is_fb(m)) {
    if (isTRUE(refdiff$len_diff)) paste0('<div class="refwarn">&#9888; <b>The cDNA on this page is not the same length as the transcript</b> ',
     html_escape(sub("^length differs ", "", refdiff$txt)), '. The ', html_escape(build %||% "genome"),
     ' exons do not add up to the transcript sequence, so c. numbers and the translation shown here may be shifted against the official numbering. ',
     'Find variants by genomic position.</div>')
    else paste0('<div class="refwarn">&#9888; <b>The ', html_escape(build %||% "genome"),
     ' reference differs from the transcript sequence</b> at ', html_escape(refdiff$txt),
     '. c. numbering is not affected. Primers over these bases may not match the transcript.</div>')
  } else "", '
', varbar, '
<nav class="tabs"><button data-view="cdna">cDNA sequence</button><button data-view="exons">Exons</button>',
  if (m$has_cds) '<button data-view="protein">Protein sequence</button>' else "", '</nav>

<section class="view" id="cdna">
', toolbar, '
<div class="key">
<b>Codons</b><div><span class="sw">ATG</span><span class="sw c1">CCT</span> alternating codons</div>
<b>Exons</b><div><span class="sw">An exon</span> <span class="sw x1">Another exon</span></div>
<b>Other</b><div><span class="sw u">UTR</span> &nbsp;<span class="sw pep">-M-</span> translation</div>
<b>Markup</b><div>', tog("Codons", "no-codons", "pre.cdna"), tog("Exons", "no-exons", "pre.cdna"),
  tog("UTR", "no-utr", "pre.cdna"), tog("CDS &amp; translation lines", "no-trans", "pre.cdna"),
  tog("Line numbers", "no-lnum", "pre.cdna"), '</div>
</div>
', render_cdna(m, width), '
</section>

<section class="view" id="exons">
', toolbar, '
<div class="key">
<b>Key</b><div><span class="sw fl">flanking</span> <span class="sw ux">UTR</span> <span class="sw cx">CDS</span> <span class="sw in">intron</span>',
  if (m$strand == "-") ' &nbsp;<span class="strand-badge">reverse strand</span> <span class="hint">Start &gt; End: coordinates run high &rarr; low; sequence is 5&prime;&rarr;3&prime; of the transcript</span>' else "", '</div>
<b>Intron sequence</b><div>', isel, ' <span class="hint">of intron shown next to each exon; the middle of longer introns appears as ..........</span></div>
<b>Markup</b><div>', tog("Introns", "no-introns", "table.exons"), tog("Flanking sequence", "no-flanks", "table.exons"),
  tog("UTR", "no-utr", "table.exons"), tog("Line numbers", "no-lnum", "table.exons"), '</div>
</div>
', render_exons(m, width), '
<p class="note">Start/End are given in transcript orientation (5&prime; &rarr; 3&prime;). Phase &minus; = exon boundary in UTR.
Intron and flanking sequence can be selected here; selections in this view are genomic (unspliced).</p>
</section>
', if (m$has_cds) paste0('
<section class="view" id="protein">
', toolbar, '
<div class="key">
<b>Exons</b><div><span class="sw">An exon</span> <span class="sw x1">Another exon</span> <span class="sw so">Residue overlaps splice site</span></div>
<b>Markup</b><div>', tog("Exons", "no-exons", "pre.prot"), tog("Splice-site residues", "no-so", "pre.prot"),
  tog("Line numbers", "no-lnum", "pre.prot"), '</div>
</div>
', render_protein(m, width), '
</section>') else "", '
', if (!is.null(prov)) paste0('<footer class="prov">', html_escape(prov), '</footer>') else "", '
</main>
', pickbar, '
<script type="application/json" id="seqdata">', page_json(m, vars, bins, build, prov, refdiff), '</script>
<script>', .PAGE_JS, '</script>
</body></html>')
}

index_html <- function(summ, title = "MANE transcripts", build = NULL, prov = NULL) {
  if (is.null(summ$refseq)) summ$refseq <- ""
  if (is.null(summ$status)) summ$status <- ""
  rows <- sprintf(
    '<tr><td><a href="%s"><b>%s</b></a></td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s:%s-%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td></tr>',
    html_escape(summ$file), html_escape(summ$gene_name),
    paste0(html_escape(summ$refseq), ifelse(is.null(summ$version_note) | summ$version_note == "", "",
           paste0(' <span class="badge old" title="', html_escape(summ$version_note), '">&#9888; NOT MANE</span>'))),
    ifelse(!is.null(summ$version_note) & summ$version_note != "",
           paste0('<span class="badge old" title="', html_escape(summ$version_note), '">', sub("^Not MANE \\(([a-z]+) version\\).*$", "\\1", summ$version_note), ' version of ', summ$status, '</span>'),
    ifelse(summ$status == "MANE Plus Clinical", '<span class="badge mane">MANE Plus Clinical</span>',
           ifelse(summ$status == "", "", '<span class="badge">MANE Select</span>'))),
    html_escape(summ$transcript_id),
    html_escape(summ$gene_id), html_escape(summ$chr), fmt_int(summ$start), fmt_int(summ$end),
    summ$strand, summ$exons, fmt_int(summ$cdna_bp), ifelse(is.na(summ$protein_aa), "-", fmt_int(summ$protein_aa)))
  paste0('<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>', title, '</title><style>',
.PAGE_CSS, '
table.idx{border-collapse:collapse;font-size:13px;margin-top:12px}
table.idx th{background:var(--head);color:#fff;text-align:left;padding:6px 12px;white-space:nowrap}
table.idx td{border-bottom:1px solid var(--line);padding:5px 12px;white-space:nowrap}
table.idx tr.alt td{background:var(--stripe)}
table.idx a{color:var(--head)}
input.q{font:inherit;padding:6px 10px;border:1px solid #9aa8b6;border-radius:3px;width:280px}
</style></head><body>
<header class="top"><h1>', title, '</h1>', if (!is.null(build)) paste0('<span class="build">', html_escape(build), '</span>') else "",
  '<span style="font-size:13px;opacity:.85">', nrow(summ), ' pages</span></header>
<main><input class="q" placeholder="Filter by gene, transcript or chromosome" autofocus>
<table class="idx"><thead><tr><th>Gene</th><th>RefSeq</th><th>MANE</th><th>Transcript</th><th>Gene ID</th><th>Location</th><th>Strand</th><th>Exons</th><th>cDNA bp</th><th>Protein aa</th></tr></thead><tbody>',
paste(rows, collapse = "\n"), '</tbody></table>',
  if (!is.null(prov)) paste0('<footer class="prov">', html_escape(prov), '</footer>') else "", '</main>
<script>
var q=document.querySelector("input.q"),rows=[].slice.call(document.querySelectorAll("table.idx tbody tr"));
function stripe(){var k=0;rows.forEach(function(r){if(r.style.display!=="none")r.classList.toggle("alt",(k++)%2===1)})}
q.addEventListener("input",function(){var v=q.value.toLowerCase();rows.forEach(function(r){r.style.display=r.textContent.toLowerCase().indexOf(v)>-1?"":"none"});stripe()});
stripe();
</script></body></html>')
}

# ---- main entry point ------------------------------------------------------

#' Write one HTML page per MANE Select transcript (named <GENE>.html) plus
#' index.html into `out_dir`. Returns a summary data.frame invisibly.
#'
#' @param variants    optional gnomAD/VCF variants (data.frame with CHROM, POS,
#'                    REF, ALT, AF; a VariantAnnotation VCF; or a vcfR object).
#' @param tags        GTF tag(s) that select transcripts. Add "MANE_Plus_Clinical"
#'                    to include those too (files then get the transcript ID
#'                    appended when a gene has more than one page).
#' @param genes       optional character vector of gene symbols / IDs to restrict to.
#' @param flank       minimum bp of up/downstream sequence on the Exons view.
#' @param intron_ends default bp shown at each end of an intron (dropdown default).
#' @param variant_label label of the checkbox that turns variant colouring on
#'                    (off by default).
#' @param min_af      variants with AF below this are dropped (default 1e-4 = 0.01%).
#' @param af_default  AF filter the pages open with (default 0.001 = 0.1%).
#' @param af_bins     data.frame(min, bg, fg, label) of AF bands; default
#'                    .AF_BINS uses Ensembl's seven variant colours.
#' @param genome_build  assembly label shown on every page, e.g. "GRCh37 (hg19)".
#' @param provenance  one line naming data sources / versions / date (page + RTF footer).
#' @param variant_source text shown when hovering the variant checkbox.
#' @param reference_cdna transcript sequences to check the built cDNA against
#'                    (path to MANE ensembl_rna FASTA, or named character vector keyed by
#'                    ENST; a refseq_rna file is keyed by NM_ and would match nothing);
#'                    differences are shown as a banner on the page.
#' @param intron_max  bp of each intron end embedded in the page = the largest
#'                    dropdown choice. Introns shorter than 2 x intron_max are
#'                    embedded whole. Bigger values = bigger files.
make_mane_pages <- function(gtf, fasta, variants = NULL, out_dir = "mane_pages",
                            tags = c("MANE_Select", "MANE_Plus_Clinical"),
                            genes = NULL, flank = 50, intron_ends = 25, intron_max = 1000, width = 60,
                            variant_label = "CHUV-specific", af_bins = .AF_BINS, min_af = 1e-4, af_default = 0.001,
                            genome_build = NULL, provenance = NULL, variant_source = NULL, reference_cdna = NULL,
                            verbose = TRUE) {
  if (is.character(reference_cdna) && length(reference_cdna) == 1 && is.null(names(reference_cdna))) {
    # an unnamed single string is a file path, never a sequence
    if (!file.exists(reference_cdna)) stop("reference_cdna file not found: ", reference_cdna)
    reference_cdna <- read_fasta_simple(reference_cdna)
  }
  if (!is.null(reference_cdna)) names(reference_cdna) <- sub("\\.[0-9]+$", "", sub("\\|.*$", "", names(reference_cdna)))
  af_bins <- af_bins[order(-af_bins$min), ]
  # only bands above the AF cut-off are kept (and drawn in the legend)
  af_bins <- af_bins[c(Inf, utils::head(af_bins$min, -1)) > min_af, , drop = FALSE]
  df <- normalize_gtf(gtf)
  tag_rx <- paste0("(^|[,;| ])(", paste(tags, collapse = "|"), ")($|[,;| ])")
  hit <- !is.na(df$tag) & grepl(tag_rx, df$tag)
  tx_ids <- unique(df$transcript_id[hit & !is.na(df$transcript_id)])
  if (!length(tx_ids)) {
    stop("No transcripts tagged ", paste(tags, collapse = "/"), " found in the GTF. ",
         "If you imported with rtracklayer, check that the `tag` column kept all tags ",
         "(some versions keep only one). read_gtf_simple() keeps all of them.")
  }
  feats <- c("exon", "CDS", "start_codon", "stop_codon", "Selenocysteine", "transcript")
  df <- df[df$transcript_id %in% tx_ids & df$feature %in% feats, , drop = FALSE]
  if (!is.null(genes)) {
    df <- df[df$gene_name %in% genes | df$gene_id %in% genes | sub("\\.[0-9]+$", "", df$gene_id) %in% genes, , drop = FALSE]
    tx_ids <- intersect(tx_ids, unique(df$transcript_id))
  }
  by_tx <- split(df, df$transcript_id)[tx_ids]
  get_seq <- make_seq_getter(fasta)
  vby <- normalize_variants(variants)
  # chromosomes the variant data covers (e.g. gnomAD v2 has no chrY); pages on
  # uncovered chromosomes say so in the variant bar and the tooltip names them
  var_chroms <- if (!is.null(vby)) names(vby)[vapply(vby, nrow, 0L) > 0] else character(0)
  if (!is.null(vby) && min_af > 0) vby <- lapply(vby, function(d) d[d$af >= min_af, , drop = FALSE])
  if (verbose && !is.null(vby)) message("Variants loaded (AF >= ", min_af, "): ", fmt_int(sum(vapply(vby, nrow, 0L))))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  models <- list()
  for (t in tx_ids) {
    m <- tryCatch(build_transcript_model(by_tx[[t]], get_seq, flank = flank, embed = intron_max),
                  error = function(e) { warning("Skipping ", t, ": ", conditionMessage(e), call. = FALSE); NULL })
    if (!is.null(m)) models[[t]] <- m
  }
  if (!length(models)) stop("No pages could be built")
  # only chromosomes that have a page are named: the variants may have been read for
  # the page regions only (trial run with `genes`), so other chromosomes say nothing
  no_chr <- setdiff(unique(.chr_key(vapply(models, `[[`, "", "chr"))), var_chroms)
  if (!is.null(vby) && length(no_chr))
    variant_source <- paste0(variant_source %||% "Variant data", " (no data for ", paste0("chr", no_chr, collapse = ", "), ")")

  sym <- vapply(models, function(m) m$gene_name, "")
  # file names: gene symbol + RefSeq accession without its version (BRCA2_NM_000059.html) when the
  # transcript has one; otherwise the gene symbol (+ ENST if a gene has two pages).
  # Only the file name drops the version - pages always show the full ID.
  refseq <- vapply(models, function(m) if (!is.na(m$transcript_name) && grepl("^[NX][MR]_[0-9]+", m$transcript_name))
                     m$transcript_name else NA_character_, "")
  status <- vapply(models, function(m) if ("MANE_Plus_Clinical" %in% m$tags) "MANE Plus Clinical"
                     else if ("MANE_Select" %in% m$tags) "MANE Select" else "", "")
  dup <- sym %in% sym[duplicated(sym)]
  files <- ifelse(!is.na(refseq), paste0(safe_filename(sym), "_", safe_filename(sub("\\.[0-9]+$", "", refseq)), ".html"),
           ifelse(dup, paste0(safe_filename(sym), "_", safe_filename(names(models)), ".html"),
                  paste0(safe_filename(sym), ".html")))
  files <- paste0(make.unique(sub("\\.html$", "", files), sep = "_"), ".html")

  summ <- data.frame(gene_name = sym, transcript_id = vapply(models, `[[`, "", "transcript_id"),
                     gene_id = vapply(models, `[[`, "", "gene_id"),
                     chr = vapply(models, `[[`, "", "chr"),
                     start = vapply(models, `[[`, 0, "tx_start"), end = vapply(models, `[[`, 0, "tx_end"),
                     strand = vapply(models, `[[`, "", "strand"),
                     exons = vapply(models, function(m) nrow(m$exons), 0L),
                     cdna_bp = vapply(models, function(m) nchar(m$cdna), 0L),
                     protein_aa = vapply(models, function(m) if (m$has_cds) nchar(sub("\\*$", "", m$protein)) else NA_integer_, 0L),
                     variants = 0L, ref_diff = "", version_note = "", refseq = ifelse(is.na(refseq), "", refseq), status = status,
                     file = files, stringsAsFactors = FALSE, row.names = NULL)

  for (k in seq_along(models)) {
    vars <- variants_for_model(vby, models[[k]], get_seq)
    attr(vars, "supplied") <- !is.null(vby)
    summ$variants[k] <- nrow(vars)
    rd <- ref_diff(models[[k]], reference_cdna)
    ref_ok <- !is.null(reference_cdna) && sub("\\.[0-9]+$", "", .ref_key(models[[k]])) %in% names(reference_cdna)
    summ$version_note[k] <- if (.is_fb(models[[k]]))
      paste0("Not MANE (", .fb_word(models[[k]]), " version): MANE is ", models[[k]]$mane_refseq, " / ", models[[k]]$mane_enst, ", not on this build") else ""
    summ$ref_diff[k] <- if (is.null(rd)) "" else rd$txt
    writeLines(page_html(models[[k]], vars, width, intron_ends, bins = af_bins, variant_label = variant_label,
                         af_default = af_default, build = genome_build, prov = provenance,
                         var_source = variant_source, refdiff = rd, ref_checked = ref_ok,
                         var_nochr = !is.null(vby) && !(.chr_key(models[[k]]$chr) %in% var_chroms)),
               file.path(out_dir, files[k]), useBytes = TRUE)
    if (verbose) message(sprintf("  %-12s %-20s %5d variants -> %s", sym[k], summ$transcript_id[k], nrow(vars), files[k]))
  }
  ord <- order(summ$gene_name)
  writeLines(index_html(summ[ord, ], build = genome_build, prov = provenance), file.path(out_dir, "index.html"), useBytes = TRUE)
  if (verbose) message(sprintf("Wrote %d page(s) + index.html to %s", nrow(summ), normalizePath(out_dir)))
  invisible(structure(summ, models = models))
}
