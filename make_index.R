#!/usr/bin/env Rscript
# Simple, fast index for the MANE pages, with risky pages flagged.
#
#   Rscript make_index.R mane_pages_hg19          (writes mane_pages_hg19/index.html)
#   or in R: source("make_index.R"); make_simple_index("mane_pages_hg19")
#
# Reads from the pages folder (all optional except the .html files):
#   pages.csv                  page list written by run_hg19.R
#   qc_cdna_vs_mane.csv        hg19 cDNA vs the MANE sequence
#   mane_version_fallback.csv  pages built from another RefSeq version (not MANE)
#   mane_not_on_grch37.csv     MANE transcripts with no page
#
# Flags
#   RED   length differs  hg19 cDNA is longer/shorter than the NM: c. numbers on the page
#                         may not match official NM numbering. Use genomic positions.
#   AMBER N bases differ  same length; some hg19 bases differ from the NM (c. numbering OK)
#   AMBER NOT MANE        page is another RefSeq version of the MANE transcript
#   GREY  not checked     no QC row for this page

make_simple_index <- function(dir, out = file.path(dir, "index.html"), title = "MANE transcripts (hg19)") {
  rd <- function(f) { p <- file.path(dir, f); if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE) else NULL }
  esc <- function(x) { x <- as.character(x); x[is.na(x)] <- ""
    x <- gsub("&", "&amp;", x, fixed = TRUE); x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE); gsub('"', "&quot;", x, fixed = TRUE) }
  col <- function(d, k, n = NROW(d)) if (!is.null(d) && k %in% names(d)) d[[k]] else rep(NA, n)

  pg <- rd("pages.csv")
  if (is.null(pg)) {                       # no pages.csv: list the files, parse GENE_NM_xxx.html
    f <- setdiff(list.files(dir, pattern = "\\.html$"), basename(out))
    f <- setdiff(f, "index.html")
    b <- sub("\\.html$", "", f)
    # gene = everything before the accession (symbols like HLA-A become HLA_A in file names)
    acc_rx <- "_([NX][MR]_[0-9]+|ENS[A-Z]*T[0-9]+).*$"
    pg <- data.frame(file = f, gene_name = sub(acc_rx, "", b),
                     refseq = ifelse(grepl("_[NX][MR]_[0-9]+", b), sub("^.*_([NX][MR]_[0-9]+).*$", "\\1", b), ""),
                     transcript_id = rep("", length(f)), status = rep("", length(f)), chr = rep("", length(f)),
                     version_note = rep("", length(f)), stringsAsFactors = FALSE)
    message("pages.csv not found - listing ", length(f), " html files (QC flags need pages.csv)")
  }
  pg <- pg[file.exists(file.path(dir, pg$file)), , drop = FALSE]
  n <- nrow(pg)

  qc <- rd("qc_cdna_vs_mane.csv"); fb <- rd("mane_version_fallback.csv"); miss <- rd("mane_not_on_grch37.csv")
  i  <- if (is.null(qc)) rep(NA_integer_, n) else match(pg$transcript_id, qc$transcript)
  built <- col(qc, "built_bp")[i]; mane <- col(qc, "mane_bp")[i]; mm <- col(qc, "mismatches")[i]

  len_diff <- !is.na(built) & !is.na(mane) & built != mane
  base_diff <- !len_diff & !is.na(mm) & mm > 0
  unchecked <- is.na(i) | is.na(mane)
  not_mane <- (!is.na(col(pg, "version_note")) & nzchar(col(pg, "version_note"))) |
              (!is.null(fb) & pg$refseq %in% col(fb, "GRCh37_RefSeq"))
  mane_id <- if (is.null(fb)) rep("", n) else ifelse(is.na(m <- match(pg$refseq, fb$GRCh37_RefSeq)), "", fb$MANE_RefSeq[m])

  flag <- function(on, cls, txt, tip) ifelse(on, sprintf('<span class="f %s" title="%s">%s</span>', cls, esc(tip), txt), "")
  flags <- paste0(
    flag(len_diff, "red", "length differs",
         sprintf("hg19 cDNA %s bp vs NM %s bp. c. numbers on this page may not match the official NM numbering - find variants by genomic position.", built, mane)),
    flag(base_diff, "amb", paste0(mm, ifelse(mm == 1, " base differs", " bases differ")),
         "Same length as the NM; these hg19 bases differ from it (listed on the page). c. numbering matches the NM."),
    flag(not_mane, "amb", "NOT MANE",
         paste0("Built from another RefSeq version", ifelse(nzchar(mane_id), paste0("; MANE is ", mane_id), ""), ".")),
    flag(unchecked & !is.null(qc), "gry", "not checked", "No QC comparison with the MANE sequence for this page."))
  risky <- len_diff | base_diff | not_mane
  o <- order(toupper(pg$gene_name), pg$file)

  # The list is kept as data and only the matching rows are drawn (hiding thousands of
  # <li> elements is what made searching slow).
  js <- function(x) paste0('"', gsub("</", "<\\/", gsub('"', '\\"', gsub("\\", "\\\\", x, fixed = TRUE), fixed = TRUE), fixed = TRUE), '"')
  chr <- col(pg, "chr"); chr <- ifelse(is.na(chr) | chr == "", "", paste0("chr", sub("^chr", "", chr)))
  keys <- tolower(paste(pg$gene_name, pg$refseq, pg$transcript_id, chr))[o]
  li <- sprintf('<li><a href="%s"><b>%s</b> %s</a> <span class="s">%s%s</span> %s</li>',
                esc(pg$file), esc(pg$gene_name), esc(pg$refseq),
                esc(col(pg, "status")), ifelse(nzchar(chr), paste0(" &middot; ", esc(chr)), ""),
                flags)[o]

  if (!is.null(miss) && nrow(miss)) miss <- miss[order(toupper(miss$symbol), miss$RefSeq_nuc), , drop = FALSE]
  nopage <- if (!is.null(miss) && nrow(miss))
    paste0('<details><summary>', nrow(miss), ' MANE transcript', if (nrow(miss) > 1) "s", ' with no page (not placed on GRCh37)</summary><ul class="np">',
           paste(sprintf("<li><b>%s</b> %s %s</li>", esc(miss$symbol), esc(miss$RefSeq_nuc), esc(miss$MANE_status)), collapse = ""),
           "</ul></details>") else ""

  html <- paste0('<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>', esc(title), '</title>
<style>
body{font:14px/1.45 system-ui,-apple-system,Segoe UI,Roboto,sans-serif;margin:0 auto;max-width:900px;padding:16px;color:#1d2a36;background:#fff}
h1{font-size:20px;margin:0 0 10px}
.bar{position:sticky;top:0;background:#fff;padding:8px 0;border-bottom:1px solid #dde3e9;display:flex;gap:14px;align-items:center;flex-wrap:wrap}
input[type=search]{font:inherit;padding:7px 10px;width:260px;max-width:100%;border:1px solid #9aa8b6;border-radius:4px}
#n{color:#5b6b7a;font-size:13px}
ul{list-style:none;margin:6px 0;padding:0}
li{padding:4px 2px;border-bottom:1px solid #eef1f4}
a{color:#1f4e79;text-decoration:none} a:hover{text-decoration:underline}
.s{color:#5b6b7a;font-size:12.5px}
.f{display:inline-block;font-size:11.5px;font-weight:600;border-radius:3px;padding:0 6px;margin-left:4px;cursor:help}
.red{background:#fde2e2;color:#8a1c1c;border:1px solid #f1aeb5}
.amb{background:#fff3bf;color:#6b4e00;border:1px solid #e8b10c}
.gry{background:#eef1f4;color:#5b6b7a;border:1px solid #cfd7df}
.key{font-size:12.5px;color:#5b6b7a;margin:8px 0}
details{margin-top:16px;font-size:13px} .np li{color:#5b6b7a}
</style></head><body>
<h1>', esc(title), '</h1>
<div class="key"><span class="f red">length differs</span> c. numbers may not match the NM &ndash; use genomic positions &nbsp;
<span class="f amb">bases differ</span> some hg19 bases differ, numbering OK &nbsp;
<span class="f amb">NOT MANE</span> another RefSeq version. Hover a flag for details.</div>
<div class="bar"><input type="search" id="q" placeholder="Search gene, NM, ENST, chr" autofocus>
<label><input type="checkbox" id="rk"> flagged only (', sum(risky), ')</label><span id="n"></span></div>
<ul id="L"></ul><div id="more" class="s"></div>', nopage, '
<script>
var K=[', paste(js(keys), collapse = ","), '];
var H=[', paste(js(li), collapse = ",\n"), '];
var RK=[', paste(as.integer(risky[o]), collapse = ","), '];
var MAX=300,L=document.getElementById("L"),q=document.getElementById("q"),rk=document.getElementById("rk"),
    N=document.getElementById("n"),M=document.getElementById("more"),t=null;
function run(){var v=q.value.trim().toLowerCase(),f=rk.checked,out=[],c=0;
  for(var i=0;i<K.length;i++){if((!f||RK[i])&&(!v||K[i].indexOf(v)>-1)){c++;if(out.length<MAX)out.push(H[i])}}
  L.innerHTML=out.join("");N.textContent=c+" of "+K.length;
  M.textContent=c>MAX?"Showing the first "+MAX+" of "+c+" - type more to narrow the list.":""}
q.addEventListener("input",function(){clearTimeout(t);t=setTimeout(run,60)});
rk.addEventListener("change",run);run();
</script></body></html>')
  writeLines(html, out, useBytes = TRUE)
  message(sprintf("Index: %d pages (%d length differs, %d bases differ, %d not MANE) -> %s",
                  n, sum(len_diff), sum(base_diff), sum(not_mane), out))
  invisible(data.frame(file = pg$file, len_diff, base_diff, not_mane, unchecked))
}

if (sys.nframe() == 0L) {
  a <- commandArgs(trailingOnly = TRUE)
  if (!length(a)) stop("usage: Rscript make_index.R <pages folder> [output.html]")
  make_simple_index(a[1], if (length(a) > 1) a[2] else file.path(a[1], "index.html"))
}
