#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 3.Create_Sample_HTML.R   (base-R edition, no htmltools)
#
# Per-sample HTML report. Uses ONLY packages already present in the project
# SIFs (R_Genomics_v1_blast.def, Viral_Genomics_v5.def):
#       data.table, optparse  (and base R for HTML emission)
#
# Charts are emitted as hand-built inline SVG, so no plotting package is
# required. FastQC archives are read with utils::unz(), also base R.
#
# Modes (--html_mode):
#   single  : one HTML file with external image references and lazy table
#   multi   : index + per-chromosome pages + circos + QC + table pages (low-RAM)
#   both    : both (default; matches the original request)
#
# Inputs:
#   --integrations  CSV/TSV, site-level. At least one of each pair below:
#                       chrom        | chromosome
#                       position     | integration_site
#                       sample_id    | sample
#                       clonal_id    | clone_id
#                       read_support | reads_at_site
#                   plus optional: timepoint, viral_orientation, gene_name,
#                   in_repeat_t2t, in_repeat_hg38, MATCH_TYPE, IPDA_INTACT
#   --master_tsv    OPTIONAL read-level parseSAM.pl output
#                   (<prefix>_MasterOfMasterFrame.tsv). Supplies READ,
#                   FLANK_TYPE, HIV_POS, HIV_LEN, HOST_SEQ_OVERLAP, HIV_SEQ.
#                   Read-level counts come from here when available, since
#                   --integrations is one row per SITE, not per read.
#   --fastqc_zip    OPTIONAL FastQC archive. If omitted and --project_dir is
#                   given, resolves
#                   <project_dir>/output/00_QualityControl/fastqc/<sample>_fastqc.zip
#   --project_dir   OPTIONAL project root used to locate the FastQC archive
#   --circos_png    PNG produced by make_circos.R. Embedded inline when it fits
#                   under --embed_max_kb, otherwise referenced externally.
#   --sample_id     Sample identifier (used in titles + file names)
#   --outdir        Where to write HTML(s)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--integrations",  type = "character"),
  make_option("--master_tsv",    type = "character", default = NA),
  make_option("--fastqc_zip",    type = "character", default = NA),
  make_option("--project_dir",   type = "character", default = NA),
  make_option("--circos_png",    type = "character", default = NA),
  make_option("--sample_id",     type = "character", default = "sample"),
  make_option("--html_mode",     type = "character", default = "both"),
  make_option("--embed_max_kb",  type = "double",    default = 2048),
  make_option("--outdir",        type = "character", default = ".")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dt <- fread(opt$integrations)

# ---- Column aliasing --------------------------------------------------------
# The pipeline's combined.csv uses upstream-friendly names. The HTML and
# circos scripts internally prefer the canonical names below. Alias whichever
# is present without dropping anything else.
.alias <- list(
  # Prefer the integer column (integration_pos / canonical_pos) over the
  # "chr1:40058490" string in `integration_site`, since downstream code
  # coerces position to integer.
  chrom        = c("chrom",        "chromosome",       "Chrom",      "Chromosome"),
  position     = c("position",     "integration_pos",  "canonical_pos",
                   "Position",     "integration_site"),
  sample_id    = c("sample_id",    "sample",           "Sample",     "SampleID"),
  clonal_id    = c("clonal_id",    "clone_id",         "Clone_ID",   "ClonalID",
                   "CLONE_ID"),
  read_support = c("read_support", "reads_at_site",    "ReadSupport","Reads_At_Site")
)

apply_aliases <- function(d) {
  for (canon in names(.alias)) {
    if (canon %in% colnames(d)) next
    candidates <- intersect(.alias[[canon]], colnames(d))
    if (length(candidates) > 0) {
      setnames(d, candidates[1], canon)
    }
  }
  d
}
dt <- apply_aliases(dt)

# Coerce types where possible
if ("position" %in% colnames(dt)) {
  suppressWarnings(dt[, position := as.integer(position)])
}
if ("chrom" %in% colnames(dt)) {
  dt[, chrom := as.character(chrom)]
}

# ---- Sequence-bearing columns ----------------------------------------------
# HIV_SEQ and HOST_SEQ_OVERLAP can each be multi-kb. They must never be pasted
# into a table cell; instead derive compact summaries and offer the raw
# sequences as a FASTA sidecar.
.SEQ_COLS <- c("HIV_SEQ", "HOST_SEQ_OVERLAP")

.is_seq <- function(x) {
  !is.na(x) & nzchar(x) & x != "NA" & x != "-"
}

.gc_pct <- function(x) {
  gc <- nchar(gsub("[^GCgc]", "", x))
  at <- nchar(gsub("[^ACGTacgt]", "", x))
  ifelse(at > 0, round(100 * gc / at, 1), NA_real_)
}

.n_pct <- function(x) {
  n <- nchar(gsub("[^Nn]", "", x))
  l <- nchar(x)
  ifelse(l > 0, round(100 * n / l, 1), NA_real_)
}

add_seq_derivations <- function(d) {
  for (sc in intersect(.SEQ_COLS, colnames(d))) {
    v <- as.character(d[[sc]])
    ok <- .is_seq(v)
    d[[paste0(sc, "_LEN")]] <- ifelse(ok, nchar(v), NA_integer_)
    d[[paste0(sc, "_GC")]]  <- ifelse(ok, .gc_pct(v), NA_real_)
    d[[paste0(sc, "_N")]]   <- ifelse(ok, .n_pct(v),  NA_real_)
  }
  d
}

# ---- Read-level table (parseSAM.pl output) ---------------------------------
mdt <- NULL
if (!is.na(opt$master_tsv) && nzchar(opt$master_tsv) && file.exists(opt$master_tsv)) {
  mdt <- fread(opt$master_tsv, sep = "\t", header = TRUE)
  mdt <- add_seq_derivations(mdt)
  message("[3.Create_Sample_HTML.R] read-level table: ", nrow(mdt), " rows, ",
          ncol(mdt), " columns")
} else if (!is.na(opt$master_tsv) && nzchar(opt$master_tsv)) {
  message("[3.Create_Sample_HTML.R] --master_tsv not found: ", opt$master_tsv)
}
dt <- add_seq_derivations(dt)

# ---- Tiny HTML helpers ------------------------------------------------------
.esc <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x[is.na(x)] <- ""
  x
}

.fmt <- function(n) {
  # vectorised: count_card()/svg_barh() pass whole count vectors
  if (length(n) == 0) return("n/a")
  n   <- suppressWarnings(as.numeric(n))
  out <- rep("n/a", length(n))
  ok  <- !is.na(n)
  if (any(ok)) out[ok] <- formatC(n[ok], format = "d", big.mark = ",")
  out
}

.css <- '
:root { --fg:#222; --muted:#666; --bg:#fff; --accent:#1f6feb; --row:#f6f8fa;
        --pass:#1a7f37; --warn:#9a6700; --fail:#cf222e; --band:#e6edf7; }
* { box-sizing: border-box; }
body { font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;
       color:var(--fg); background:var(--bg); margin:0 auto; max-width:1100px;
       padding:16px 24px; line-height:1.45; }
h1 { font-size:22px; margin:4px 0 12px; }
h2 { font-size:16px; margin:18px 0 8px; color:var(--accent); }
h3 { font-size:13.5px; margin:10px 0 6px; color:var(--fg); }
.meta { color:var(--muted); font-size:13px; margin-bottom:16px; }
.grid { display:grid; grid-template-columns:1fr 1fr; gap:12px; }
.grid3 { display:grid; grid-template-columns:1fr 1fr 1fr; gap:12px; }
.card { border:1px solid #e1e4e8; border-radius:6px; padding:12px; background:#fff; }
.card.full { grid-column:1 / -1; }
.kpis { display:grid; grid-template-columns:repeat(auto-fit,minmax(150px,1fr));
        gap:10px; margin:8px 0 4px; }
.kpi { border:1px solid #e1e4e8; border-radius:6px; padding:10px 12px;
       background:var(--row); }
.kpi .v { font-size:21px; font-weight:600; line-height:1.15; }
.kpi .l { font-size:11.5px; color:var(--muted); text-transform:uppercase;
          letter-spacing:.03em; margin-top:2px; }
.kpi .s { font-size:11.5px; color:var(--muted); margin-top:3px; }
table { border-collapse:collapse; width:100%; font-size:12.5px; }
th, td { padding:5px 8px; text-align:left; border-bottom:1px solid #eee;
         white-space:nowrap; }
th { background:var(--row); position:sticky; top:0; }
tr:hover td { background:#fafbfc; }
table.compact th, table.compact td { white-space:normal; }
.nav { margin:12px 0 16px; font-size:13px; }
.nav a { color:var(--accent); text-decoration:none; margin-right:10px; }
img.circos { max-width:100%; height:auto; border:1px solid #e1e4e8;
             border-radius:4px; display:block; margin:0 auto; }
.hint { color:var(--muted); font-size:12px; margin-top:6px; }
.scroll { max-height:60vh; overflow:auto; border:1px solid #e1e4e8; border-radius:6px; }
.badge { display:inline-block; padding:1px 7px; border-radius:10px;
         font-size:11px; font-weight:600; color:#fff; }
.b-pass { background:var(--pass); } .b-warn { background:var(--warn); }
.b-fail { background:var(--fail); } .b-na { background:#8b949e; }
code.seq { font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:11px;
           word-break:break-all; color:var(--muted); }
svg { max-width:100%; height:auto; }
'

.js <- '
function filt(id){
  var q=document.getElementById(id+"-q").value.toLowerCase();
  var rows=document.querySelectorAll("#"+id+" tbody tr");
  rows.forEach(function(r){
    r.style.display=r.textContent.toLowerCase().indexOf(q)>=0?"":"none";
  });
}
'

write_html <- function(path, title, body_html) {
  doc <- paste0(
    '<!doctype html><html><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width,initial-scale=1"><title>',
    .esc(title), '</title><style>', .css, '</style><script>', .js,
    '</script></head><body>', body_html, '</body></html>'
  )
  writeLines(doc, path, useBytes = TRUE)
}

# ---- Inline SVG chart helpers (no plotting package) ------------------------
.PAL <- c("#1f6feb", "#2da44e", "#bf8700", "#cf222e", "#8250df",
          "#0969da", "#57606a", "#1a7f37", "#a40e26", "#6e40c9")

svg_barh <- function(labels, values, width = 520, bar_h = 20, gap = 6,
                     lab_w = 190, pal = .PAL) {
  ###########################################################################
  ## Horizontal bar chart for categorical counts. Emits a self-contained    ##
  ## <svg> string; nothing but base R string handling is involved.          ##
  ###########################################################################
  n <- length(values)
  if (n == 0) return('<div class="hint">No data.</div>')
  values <- as.numeric(values)
  vmax <- max(values, na.rm = TRUE)
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  plot_w <- width - lab_w - 60
  h <- n * (bar_h + gap) + gap
  parts <- character(0)
  for (i in seq_len(n)) {
    y  <- gap + (i - 1) * (bar_h + gap)
    bw <- max(1, round(plot_w * values[i] / vmax))
    col <- pal[((i - 1) %% length(pal)) + 1]
    parts <- c(parts, paste0(
      '<text x="', lab_w - 8, '" y="', y + bar_h - 5,
      '" text-anchor="end" font-size="11.5" fill="#222">',
      .esc(labels[i]), '</text>',
      '<rect x="', lab_w, '" y="', y, '" width="', bw, '" height="', bar_h,
      '" rx="2" fill="', col, '"/>',
      '<text x="', lab_w + bw + 6, '" y="', y + bar_h - 5,
      '" font-size="11.5" fill="#444">', .fmt(values[i]), '</text>'
    ))
  }
  paste0('<svg viewBox="0 0 ', width, ' ', h,
         '" width="100%" role="img" xmlns="http://www.w3.org/2000/svg">',
         paste(parts, collapse = ""), '</svg>')
}

svg_hist <- function(labels, values, width = 700, height = 220,
                     xlab = "", ylab = "", fill = "#1f6feb") {
  ###########################################################################
  ## Vertical histogram for distributions (read length, quality). x labels  ##
  ## are thinned so they stay legible however many bins FastQC produced.    ##
  ###########################################################################
  n <- length(values)
  if (n == 0) return('<div class="hint">No data.</div>')
  values <- as.numeric(values)
  vmax <- max(values, na.rm = TRUE)
  if (!is.finite(vmax) || vmax <= 0) vmax <- 1
  ml <- 58; mr <- 10; mt <- 10; mb <- 40
  pw <- width - ml - mr
  ph <- height - mt - mb
  bw <- pw / n
  every <- max(1, ceiling(n / 12))
  bars <- character(0)
  for (i in seq_len(n)) {
    bh <- ph * values[i] / vmax
    x  <- ml + (i - 1) * bw
    y  <- mt + ph - bh
    bars <- c(bars, paste0('<rect x="', round(x + bw * 0.08, 2),
                           '" y="', round(y, 2),
                           '" width="', round(bw * 0.84, 2),
                           '" height="', round(bh, 2),
                           '" fill="', fill, '"><title>',
                           .esc(labels[i]), ': ', .fmt(values[i]),
                           '</title></rect>'))
    if (((i - 1) %% every) == 0) {
      bars <- c(bars, paste0(
        '<text x="', round(x + bw / 2, 2), '" y="', mt + ph + 14,
        '" font-size="10" text-anchor="middle" fill="#666" ',
        'transform="rotate(-35 ', round(x + bw / 2, 2), ' ', mt + ph + 14, ')">',
        .esc(labels[i]), '</text>'))
    }
  }
  # y axis: 0 and max
  axes <- paste0(
    '<line x1="', ml, '" y1="', mt, '" x2="', ml, '" y2="', mt + ph,
    '" stroke="#d0d7de"/>',
    '<line x1="', ml, '" y1="', mt + ph, '" x2="', width - mr, '" y2="', mt + ph,
    '" stroke="#d0d7de"/>',
    '<text x="', ml - 6, '" y="', mt + 9,
    '" font-size="10" text-anchor="end" fill="#666">', .fmt(vmax), '</text>',
    '<text x="', ml - 6, '" y="', mt + ph,
    '" font-size="10" text-anchor="end" fill="#666">0</text>',
    if (nzchar(ylab)) paste0(
      '<text x="12" y="', mt + ph / 2,
      '" font-size="10.5" fill="#666" text-anchor="middle" transform="rotate(-90 12 ',
      mt + ph / 2, ')">', .esc(ylab), '</text>') else "",
    if (nzchar(xlab)) paste0(
      '<text x="', ml + pw / 2, '" y="', height - 2,
      '" font-size="10.5" fill="#666" text-anchor="middle">',
      .esc(xlab), '</text>') else ""
  )
  paste0('<svg viewBox="0 0 ', width, ' ', height,
         '" width="100%" role="img" xmlns="http://www.w3.org/2000/svg">',
         axes, paste(bars, collapse = ""), '</svg>')
}

.kpi <- function(value, label, sub = NULL) {
  paste0('<div class="kpi"><div class="v">', .esc(value), '</div>',
         '<div class="l">', .esc(label), '</div>',
         if (!is.null(sub)) paste0('<div class="s">', .esc(sub), '</div>') else "",
         '</div>')
}

# ---- Viral read accounting --------------------------------------------------
# "Reads that hit viral" is a READ-level quantity, so prefer the parseSAM.pl
# master table. Fall back to the site-level table when it carries read columns.
viral_read_stats <- function(mdt, dt, fq) {
  src <- if (!is.null(mdt)) mdt else dt
  from_master <- !is.null(mdt)
  s <- list(from_master = from_master)

  read_col <- intersect(c("READ", "read", "ccs_id", "READ_ID"), colnames(src))
  s$n_rows <- nrow(src)
  s$n_reads <- if (length(read_col))
                 length(unique(stats::na.omit(src[[read_col[1]]]))) else NA_integer_

  # A read "hit viral" if it carries a viral alignment: a usable HIV_SEQ, or a
  # positive HIV_LEN, or a resolved HIV_POS.
  hit <- rep(FALSE, nrow(src))
  if ("HIV_SEQ" %in% colnames(src)) hit <- hit | .is_seq(as.character(src$HIV_SEQ))
  if ("HIV_LEN" %in% colnames(src)) {
    hl <- suppressWarnings(as.numeric(gsub("///.*$", "", as.character(src$HIV_LEN))))
    hit <- hit | (!is.na(hl) & hl > 0)
  }
  if ("HIV_POS" %in% colnames(src)) {
    hp <- as.character(src$HIV_POS)
    hit <- hit | (!is.na(hp) & nzchar(hp) & hp != "UNK" & hp != "NA")
  }
  s$n_viral_rows <- sum(hit, na.rm = TRUE)
  s$n_viral_reads <- if (length(read_col))
                       length(unique(stats::na.omit(src[[read_col[1]]][hit]))) else
                       s$n_viral_rows

  if ("HIV_SEQ_LEN" %in% colnames(src)) {
    L <- suppressWarnings(as.numeric(src$HIV_SEQ_LEN))
    s$viral_bases  <- sum(L, na.rm = TRUE)
    s$median_len   <- if (any(!is.na(L))) round(stats::median(L, na.rm = TRUE)) else NA
    s$len_vec      <- L[!is.na(L)]
  } else {
    s$viral_bases <- NA; s$median_len <- NA; s$len_vec <- numeric(0)
  }

  cl <- intersect(c("clonal_id", "CLONE_ID"), colnames(src))
  s$n_clones <- if (length(cl))
                  length(unique(stats::na.omit(src[[cl[1]]]))) else NA_integer_
  pc <- intersect(c("PCR_ID", "pcr_id"), colnames(src))
  s$n_pcr <- if (length(pc))
               length(unique(stats::na.omit(src[[pc[1]]]))) else NA_integer_

  s$total_reads_fq <- if (!is.null(fq) && !is.null(fq$total_sequences))
                        fq$total_sequences else NA_real_
  s
}

viral_kpi_strip <- function(s) {
  pct <- if (!is.na(s$total_reads_fq) && s$total_reads_fq > 0 &&
             !is.na(s$n_viral_reads))
           sprintf("%.2f%% of %s sequenced",
                   100 * s$n_viral_reads / s$total_reads_fq,
                   .fmt(s$total_reads_fq)) else NULL
  paste0(
    '<div class="kpis">',
    .kpi(.fmt(s$n_viral_reads), "Reads hitting viral", pct),
    .kpi(.fmt(s$n_reads),  "Reads analysed",
         if (s$from_master) "from parseSAM master table" else "from site table"),
    .kpi(.fmt(s$n_clones), "Unique clones"),
    .kpi(.fmt(s$n_pcr),    "PCR-collapsed clones"),
    .kpi(if (is.na(s$viral_bases)) "n/a" else .fmt(s$viral_bases),
         "Viral bases recovered",
         if (!is.na(s$median_len)) paste0("median ", .fmt(s$median_len), " bp/read") else NULL),
    '</div>'
  )
}

viral_length_card <- function(s) {
  if (!length(s$len_vec)) return("")
  L <- s$len_vec
  brks <- pretty(range(L), n = 18)
  if (length(brks) < 2) brks <- c(min(L), max(L) + 1)
  cut_i <- cut(L, breaks = brks, include.lowest = TRUE, right = FALSE)
  tabv  <- table(cut_i)
  labs  <- sub(",", "\u2013", gsub("\\[|\\]|\\(|\\)", "", names(tabv)))
  paste0('<div class="card full"><h2>Recovered viral sequence length</h2>',
         svg_hist(labs, as.integer(tabv),
                  xlab = "HIV_SEQ length (bp)", ylab = "reads"),
         '<div class="hint">Derived from the HIV_SEQ column: length of the read ',
         'segment aligning to the viral reference, in reference orientation.</div>',
         '</div>')
}

# ---- Integrity / classification breakdowns ---------------------------------
# Data-driven: render a chart for whichever classification columns exist rather
# than hard-coding one schema.
.INTEGRITY_COLS <- c(
  FLANK_TYPE        = "Flanking structure",
  MATCH_TYPE        = "Match type",
  IPDA_INTACT       = "IPDA intact",
  IPDA_V2_INTACT    = "IPDA v2 intact",
  intactness        = "Intactness",
  provirus_class    = "Provirus class",
  viral_orientation = "Viral orientation"
)

count_card <- function(src, col, title, id) {
  if (is.null(src) || !(col %in% colnames(src))) return("")
  v <- as.character(src[[col]])
  v[is.na(v) | !nzchar(v)] <- "(unspecified)"
  tb <- sort(table(v), decreasing = TRUE)
  if (!length(tb)) return("")
  tot <- sum(tb)
  rows <- paste0("<tr><td>", .esc(names(tb)), "</td><td>", .fmt(as.integer(tb)),
                 "</td><td>", sprintf("%.1f%%", 100 * as.integer(tb) / tot),
                 "</td></tr>", collapse = "")
  paste0('<div class="card"><h2>', .esc(title), '</h2>',
         svg_barh(names(tb), as.integer(tb)),
         '<table class="compact"><thead><tr><th>', .esc(col),
         '</th><th>n</th><th>%</th></tr></thead><tbody>', rows,
         '</tbody></table>',
         '<div class="hint">', .fmt(tot), ' records.</div></div>')
}

integrity_section <- function(mdt, dt) {
  cards <- character(0)
  for (col in names(.INTEGRITY_COLS)) {
    # prefer the read-level table when the column exists in both
    src <- if (!is.null(mdt) && col %in% colnames(mdt)) mdt else dt
    cards <- c(cards, count_card(src, col, .INTEGRITY_COLS[[col]],
                                 paste0("cnt-", col)))
  }
  cards <- cards[nzchar(cards)]
  if (!length(cards)) return("")
  paste0('<h2>Integrity &amp; classification</h2><div class="grid">',
         paste(cards, collapse = ""), '</div>')
}

# ---- FastQC ----------------------------------------------------------------
resolve_fastqc <- function(explicit, project_dir, sample_id) {
  if (!is.na(explicit) && nzchar(explicit)) return(explicit)
  if (is.na(project_dir) || !nzchar(project_dir)) return(NA_character_)
  d <- file.path(project_dir, "output", "00_QualityControl", "fastqc")
  direct <- file.path(d, paste0(sample_id, "_fastqc.zip"))
  if (file.exists(direct)) return(direct)
  if (dir.exists(d)) {
    # tolerate lane/suffix decorations, e.g. 115-5_S3_L001_fastqc.zip
    cands <- list.files(d, pattern = "_fastqc\\.zip$", full.names = TRUE)
    hit <- cands[startsWith(basename(cands), sample_id)]
    if (length(hit)) return(hit[1])
  }
  # Return the expected location anyway so the report can say where it looked.
  direct
}

parse_fastqc_modules <- function(lines) {
  mods <- list(); cur <- NULL; buf <- character(0)
  for (ln in lines) {
    if (startsWith(ln, ">>END_MODULE")) {
      if (!is.null(cur)) mods[[cur$name]] <- list(status = cur$status, lines = buf)
      cur <- NULL; buf <- character(0)
    } else if (startsWith(ln, ">>")) {
      p <- strsplit(sub("^>>", "", ln), "\t", fixed = TRUE)[[1]]
      cur <- list(name = p[1], status = if (length(p) > 1) p[2] else NA_character_)
      buf <- character(0)
    } else if (!is.null(cur)) {
      buf <- c(buf, ln)
    }
  }
  mods
}

module_df <- function(mod) {
  if (is.null(mod)) return(NULL)
  ln   <- mod$lines
  hdr  <- ln[startsWith(ln, "#")]
  body <- ln[!startsWith(ln, "#") & nzchar(ln)]
  if (!length(body) || !length(hdr)) return(NULL)
  cols  <- strsplit(sub("^#", "", hdr[length(hdr)]), "\t", fixed = TRUE)[[1]]
  parts <- strsplit(body, "\t", fixed = TRUE)
  ncol_ <- length(cols)
  parts <- lapply(parts, function(p) { length(p) <- ncol_; p })
  m <- do.call(rbind, parts)
  df <- as.data.frame(m, stringsAsFactors = FALSE)
  names(df) <- cols
  df
}

read_fastqc_zip <- function(zip_path) {
  if (is.na(zip_path) || !nzchar(zip_path) || !file.exists(zip_path)) return(NULL)
  lst <- tryCatch(utils::unzip(zip_path, list = TRUE), error = function(e) NULL)
  if (is.null(lst) || !nrow(lst)) {
    message("[3.Create_Sample_HTML.R] could not list FastQC archive: ", zip_path)
    return(NULL)
  }
  out <- list(path = zip_path, modules = list(), summary = NULL,
              total_sequences = NA_real_)
  dat <- grep("fastqc_data\\.txt$", lst$Name, value = TRUE)
  if (length(dat)) {
    con <- unz(zip_path, dat[1])
    lines <- tryCatch(readLines(con, warn = FALSE), error = function(e) character(0))
    close(con)
    out$modules <- parse_fastqc_modules(lines)
  }
  smy <- grep("summary\\.txt$", lst$Name, value = TRUE)
  if (length(smy)) {
    con <- unz(zip_path, smy[1])
    s <- tryCatch(readLines(con, warn = FALSE), error = function(e) character(0))
    close(con)
    s <- s[nzchar(s)]
    if (length(s)) {
      p <- strsplit(s, "\t", fixed = TRUE)
      out$summary <- data.frame(
        status = vapply(p, function(x) x[1], character(1)),
        module = vapply(p, function(x) if (length(x) > 1) x[2] else NA_character_,
                        character(1)),
        stringsAsFactors = FALSE)
    }
  }
  bs <- module_df(out$modules[["Basic Statistics"]])
  if (!is.null(bs) && all(c("Measure", "Value") %in% names(bs))) {
    ts <- bs$Value[bs$Measure == "Total Sequences"]
    if (length(ts)) out$total_sequences <- suppressWarnings(as.numeric(ts[1]))
  }
  out
}

.status_badge <- function(st) {
  cls <- switch(toupper(as.character(st)),
                PASS = "b-pass", WARN = "b-warn", FAIL = "b-fail", "b-na")
  paste0('<span class="badge ', cls, '">', .esc(st), '</span>')
}

fastqc_section <- function(fq, sample_id, zip_guess) {
  if (is.null(fq)) {
    return(paste0(
      '<div class="card full"><h2>Read distribution (FastQC)</h2>',
      '<div class="hint">No FastQC archive loaded',
      if (!is.na(zip_guess)) paste0(' (looked for ', .esc(zip_guess), ')') else
        ' (pass --fastqc_zip or --project_dir)',
      '. Expected layout: ',
      '<code>&lt;project_dir&gt;/output/00_QualityControl/fastqc/&lt;sample&gt;_fastqc.zip</code>',
      '</div></div>'))
  }

  # Basic statistics
  bs <- module_df(fq$modules[["Basic Statistics"]])
  bs_html <- ""
  if (!is.null(bs) && all(c("Measure", "Value") %in% names(bs))) {
    keep <- bs[bs$Measure %in% c("Filename", "Total Sequences", "Sequences flagged as poor quality",
                                 "Sequence length", "%GC"), , drop = FALSE]
    if (nrow(keep)) {
      bs_html <- paste0(
        '<h3>Basic statistics</h3><table class="compact"><tbody>',
        paste0("<tr><th>", .esc(keep$Measure), "</th><td>",
               .esc(keep$Value), "</td></tr>", collapse = ""),
        '</tbody></table>')
    }
  }

  # Sequence length distribution -> the read distribution chart
  sld <- module_df(fq$modules[["Sequence Length Distribution"]])
  len_html <- '<div class="hint">Sequence Length Distribution module absent.</div>'
  if (!is.null(sld) && ncol(sld) >= 2) {
    labs <- as.character(sld[[1]])
    vals <- suppressWarnings(as.numeric(sld[[2]]))
    ok <- !is.na(vals)
    if (any(ok)) {
      len_html <- paste0(
        svg_hist(labs[ok], vals[ok], xlab = "read length (bp)", ylab = "reads"),
        '<div class="hint">FastQC module status: ',
        .status_badge(fq$modules[["Sequence Length Distribution"]]$status),
        '</div>')
    }
  }

  # Per sequence quality
  psq <- module_df(fq$modules[["Per sequence quality scores"]])
  q_html <- ""
  if (!is.null(psq) && ncol(psq) >= 2) {
    labs <- as.character(psq[[1]])
    vals <- suppressWarnings(as.numeric(psq[[2]]))
    ok <- !is.na(vals)
    if (any(ok)) {
      q_html <- paste0(
        '<h3>Per-sequence quality</h3>',
        svg_hist(labs[ok], vals[ok], height = 180,
                 xlab = "mean Phred score", ylab = "reads", fill = "#2da44e"))
    }
  }

  # Module status grid
  sm_html <- ""
  if (!is.null(fq$summary) && nrow(fq$summary)) {
    sm_html <- paste0(
      '<h3>Module status</h3><table class="compact"><tbody>',
      paste0("<tr><td>", .esc(fq$summary$module), "</td><td>",
             vapply(fq$summary$status, .status_badge, character(1)),
             "</td></tr>", collapse = ""),
      '</tbody></table>')
  }

  paste0(
    '<div class="card full"><h2>Read distribution (FastQC)</h2>',
    '<div class="hint">Source: ', .esc(basename(fq$path)), '</div>',
    len_html, q_html,
    '<div class="grid" style="margin-top:12px">',
    '<div>', bs_html, '</div><div>', sm_html, '</div></div>',
    '</div>')
}

# ---- Circos ----------------------------------------------------------------
encode_base64 <- function(path) {
  raw_bytes <- readBin(path, what = "raw", n = file.info(path)$size)
  if (requireNamespace("base64enc", quietly = TRUE))
    return(base64enc::base64encode(raw_bytes))
  if (requireNamespace("openssl", quietly = TRUE))
    return(openssl::base64_encode(raw_bytes))
  # Pure-R fallback (no package needed):
  alphabet <- c(LETTERS, letters, 0:9, "+", "/")
  n <- length(raw_bytes)
  pad_n <- (3L - (n %% 3L)) %% 3L
  raw_bytes <- c(raw_bytes, as.raw(rep(0L, pad_n)))
  ints <- as.integer(raw_bytes)
  triplets <- matrix(ints, nrow = 3L)
  big <- triplets[1L,] * 65536L + triplets[2L,] * 256L + triplets[3L,]
  idx <- cbind(bitwShiftR(big, 18) %% 64,
               bitwShiftR(big, 12) %% 64,
               bitwShiftR(big,  6) %% 64,
               big                  %% 64) + 1L
  chars <- matrix(alphabet[idx], ncol = 4, byrow = FALSE)
  if (pad_n > 0L) chars[nrow(chars), (4L - pad_n + 1L):4L] <- "="
  paste0(t(chars), collapse = "")
}

circos_card <- function(circos_png, src_override = NULL, embed = TRUE,
                        embed_max_kb = 2048, full = TRUE) {
  cls <- if (full) 'card full' else 'card'
  if (!is.null(src_override)) {
    return(paste0('<div class="', cls, '"><h2>Circos</h2>',
                  '<img class="circos" src="', .esc(src_override),
                  '" alt="Per-sample circos plot"/></div>'))
  }

  if (is.na(circos_png) || !file.exists(circos_png)) {
    cwd_pngs <- list.files(".", pattern = "\\.png$", full.names = FALSE)
    msg <- sprintf("No circos PNG found at %s. PNGs in cwd: %s",
                   ifelse(is.na(circos_png), "(unset)", circos_png),
                   if (length(cwd_pngs)) paste(cwd_pngs, collapse = ", ") else "(none)")
    message("[3.Create_Sample_HTML.R] ", msg)
    return(paste0(
      '<div class="', cls, '"><h2>Circos</h2>',
      '<div class="hint">', .esc(msg),
      '. Make sure make_circos.R ran before this report.</div></div>'
    ))
  }

  size_kb <- file.info(circos_png)$size / 1024

  # Default: embed as base64 data: URL so the single-page HTML is genuinely
  # self-contained -- no "is the PNG in the published dir" failure mode.
  if (embed && size_kb < embed_max_kb) {
    b64 <- encode_base64(circos_png)
    src <- paste0("data:image/png;base64,", b64)
    return(paste0(
      '<div class="', cls, '"><h2>Circos</h2>',
      '<img class="circos" src="', src,
      '" alt="Per-sample circos plot"/>',
      '<div class="hint">Inlaid as base64 (', round(size_kb, 1),
      ' KB) &mdash; this page is self-contained.</div></div>'
    ))
  }

  paste0('<div class="', cls, '"><h2>Circos</h2>',
         '<img class="circos" src="', .esc(basename(circos_png)),
         '" alt="Per-sample circos plot"/>',
         '<div class="hint">Referenced externally (', round(size_kb, 1),
         ' KB exceeds --embed_max_kb ', embed_max_kb, ').</div></div>')
}

# ---- Render a (potentially large) table as a string -----------------------
# Choose a focused set of columns for the inset table so the user actually
# sees the integration metadata, not 38 columns including the full viral
# nucleotide sequence.
.table_columns <- function(dt) {
  preferred <- c("READ", "sample_id", "clonal_id", "PCR_ID", "chrom", "position",
                 "FLANK_TYPE", "viral_orientation", "gene_name", "gene_id",
                 "read_support", "integration_site", "SHEAR_SITE",
                 "HIV_POS", "HIV_LEN",
                 "HIV_SEQ_LEN", "HIV_SEQ_GC", "HIV_SEQ_N",
                 "HOST_SEQ_OVERLAP_LEN",
                 "MATCH_TYPE", "IPDA_INTACT", "IPDA_V2_INTACT",
                 "in_repeat_t2t", "in_repeat_hg38",
                 "repeat_name_t2t", "repeat_name_hg38")
  cols <- intersect(preferred, colnames(dt))
  if (length(cols) == 0) {
    cols <- setdiff(colnames(dt), .SEQ_COLS)   # never dump raw sequence
  }
  cols
}

render_table <- function(dt_chunk, table_id) {
  cols <- .table_columns(dt_chunk)
  if (!nrow(dt_chunk) || !length(cols))
    return('<div class="hint">No rows to display.</div>')
  thead <- paste0("<thead><tr>",
                  paste0("<th>", .esc(cols), "</th>", collapse = ""),
                  "</tr></thead>")
  cell_mat <- vapply(cols, function(c) .esc(dt_chunk[[c]]),
                     character(nrow(dt_chunk)))
  if (is.null(dim(cell_mat))) cell_mat <- matrix(cell_mat, nrow = 1)
  rows <- apply(cell_mat, 1, function(r)
    paste0("<tr>",
           paste0("<td>", r, "</td>", collapse = ""),
           "</tr>"))
  tbody <- paste0("<tbody>", paste(rows, collapse = ""), "</tbody>")
  hidden <- intersect(.SEQ_COLS, colnames(dt_chunk))
  paste0(
    '<div><input id="', table_id,
    '-q" type="search" placeholder="filter&hellip;  (try: chr1, AluY, gene name)" ',
    'oninput="filt(\'', table_id, '\')" ',
    'style="width:100%;padding:6px 8px;margin-bottom:6px;"/></div>',
    '<div class="scroll"><table id="', table_id, '">', thead, tbody, '</table></div>',
    '<div class="hint">', nrow(dt_chunk),
    ' rows. Showing ', length(cols), '/', ncol(dt_chunk), ' columns.',
    if (length(hidden))
      paste0(' Raw sequence column(s) ', paste(hidden, collapse = ", "),
             ' are summarised as _LEN/_GC/_N and exported to the FASTA sidecar ',
             'rather than printed.') else "",
    ' Filter is client-side; full data is in the underlying CSV.</div>'
  )
}

summary_card <- function(dt) {
  tps    <- if ("timepoint" %in% colnames(dt))
              paste(sort(unique(stats::na.omit(dt$timepoint))), collapse = ", ") else "n/a"
  clones <- if ("clonal_id" %in% colnames(dt))
              length(unique(stats::na.omit(dt$clonal_id))) else NA_integer_
  in_rep <- if ("in_repeat_t2t" %in% colnames(dt))
              sum(as.logical(dt$in_repeat_t2t), na.rm = TRUE) else NA_integer_
  chr_n  <- if ("chrom" %in% colnames(dt))
              length(unique(stats::na.omit(dt$chrom))) else NA_integer_
  paste0(
    '<div class="card"><h2>Summary</h2>',
    '<div>Integrations: ',  nrow(dt), '</div>',
    '<div>Unique clones: ', if (is.na(clones)) "n/a" else clones, '</div>',
    '<div>Chromosomes hit: ', if (is.na(chr_n)) "n/a" else chr_n, '</div>',
    '<div>Timepoints: ', .esc(tps), '</div>',
    '<div>In RepeatMasker repeat (T2T): ',
       if (is.na(in_rep)) "n/a" else in_rep, '</div>',
    '</div>'
  )
}

# ---- FASTA sidecar for the recovered viral sequences -----------------------
write_seq_fasta <- function(src, outdir, sample_id) {
  if (is.null(src) || !("HIV_SEQ" %in% colnames(src))) return(NA_character_)
  v  <- as.character(src$HIV_SEQ)
  ok <- .is_seq(v)
  if (!any(ok)) return(NA_character_)
  idc <- intersect(c("READ", "clonal_id", "CLONE_ID", "PCR_ID"), colnames(src))
  ids <- if (length(idc)) {
    do.call(paste, c(lapply(idc, function(c) as.character(src[[c]])[ok]), sep = "|"))
  } else paste0(sample_id, "_", which(ok))
  path <- file.path(outdir, sprintf("%s.HIV_SEQ.fasta", sample_id))
  writeLines(as.vector(rbind(paste0(">", ids), v[ok])), path)
  message("[3.Create_Sample_HTML.R] wrote ", sum(ok), " viral sequences to ", path)
  path
}

# ---- Assemble shared blocks ------------------------------------------------
zip_guess <- resolve_fastqc(opt$fastqc_zip, opt$project_dir, opt$sample_id)
fq        <- read_fastqc_zip(zip_guess)
vstats    <- viral_read_stats(mdt, dt, fq)
seq_src   <- if (!is.null(mdt)) mdt else dt
fasta_out <- write_seq_fasta(seq_src, opt$outdir, opt$sample_id)

fasta_link <- if (!is.na(fasta_out))
  paste0('<div class="hint">Recovered viral sequences: <a href="',
         .esc(basename(fasta_out)), '">', .esc(basename(fasta_out)),
         '</a></div>') else ""

# ---- Single-page mode ------------------------------------------------------
write_single <- function() {
  # Stage circos image into outdir so a relative <img src=...> resolves when
  # the PNG is too big to inline.
  if (!is.na(opt$circos_png) && file.exists(opt$circos_png)) {
    dest <- file.path(opt$outdir, basename(opt$circos_png))
    if (normalizePath(opt$circos_png, mustWork = TRUE) !=
        suppressWarnings(normalizePath(dest, mustWork = FALSE))) {
      file.copy(opt$circos_png, dest, overwrite = TRUE)
    }
  }
  body <- paste0(
    '<h1>Viral integration report &mdash; ', .esc(opt$sample_id), '</h1>',
    '<div class="meta">Generated ',
       format(Sys.time(), "%Y-%m-%d %H:%M"), '</div>',
    '<h2>Viral read accounting</h2>',
    viral_kpi_strip(vstats), fasta_link,
    '<div class="grid">', summary_card(dt), '</div>',
    integrity_section(mdt, dt),
    '<div class="grid">',
      circos_card(opt$circos_png, embed_max_kb = opt$embed_max_kb),
      viral_length_card(vstats),
      fastqc_section(fq, opt$sample_id, zip_guess),
    '</div>',
    '<h2>Integration sites</h2>', render_table(dt, "tbl-main"),
    if (!is.null(mdt)) paste0('<h2>Read-level detail</h2>',
                              render_table(mdt, "tbl-reads")) else ""
  )
  out <- file.path(opt$outdir, sprintf("%s.report.html", opt$sample_id))
  write_html(out, sprintf("%s report", opt$sample_id), body)
  message("[3.Create_Sample_HTML.R] wrote single-page: ", out, " (",
          round(file.info(out)$size / 1024, 1), " KB)")
}

# ---- Multi-page mode -------------------------------------------------------
write_multi <- function() {
  multi_dir <- file.path(opt$outdir, sprintf("%s_report", opt$sample_id))
  dir.create(multi_dir, showWarnings = FALSE, recursive = TRUE)

  if (!is.na(opt$circos_png) && file.exists(opt$circos_png))
    file.copy(opt$circos_png, file.path(multi_dir, "circos.png"), overwrite = TRUE)
  if (!is.na(fasta_out) && file.exists(fasta_out))
    file.copy(fasta_out, file.path(multi_dir, basename(fasta_out)), overwrite = TRUE)

  # Filter out NA chromosomes from the per-chromosome split
  if ("chrom" %in% colnames(dt)) {
    chrs <- sort(unique(stats::na.omit(dt$chrom)))
    chrs <- chrs[chrs != "" & chrs != "NA"]
  } else {
    chrs <- character()
  }
  nav <- paste0(
    '<div class="nav">',
    '<a href="index.html">Index</a>',
    '<a href="circos.html">Circos</a>',
    '<a href="qc.html">Read QC</a>',
    '<a href="all.html">All sites</a>',
    if (!is.null(mdt)) '<a href="reads.html">Reads</a>' else "",
    paste0('<a href="chrom_', .esc(chrs), '.html">', .esc(chrs), '</a>',
           collapse = ""),
    '</div>'
  )

  body_idx <- paste0(
    '<h1>Viral integration report &mdash; ', .esc(opt$sample_id), '</h1>',
    '<div class="meta">Multi-page report. Generated ',
       format(Sys.time(), "%Y-%m-%d %H:%M"), '</div>', nav,
    '<h2>Viral read accounting</h2>',
    viral_kpi_strip(vstats), fasta_link,
    '<div class="grid">', summary_card(dt),
       circos_card(NA, src_override = "circos.png", full = FALSE), '</div>',
    integrity_section(mdt, dt),
    '<div class="grid">', viral_length_card(vstats), '</div>'
  )
  write_html(file.path(multi_dir, "index.html"),
             sprintf("%s index", opt$sample_id), body_idx)

  write_html(file.path(multi_dir, "circos.html"), "Circos",
             paste0('<h1>Circos plot</h1>', nav,
                    '<div class="card"><img class="circos" src="circos.png"/></div>'))

  write_html(file.path(multi_dir, "qc.html"), "Read QC",
             paste0('<h1>Read QC &mdash; ', .esc(opt$sample_id), '</h1>', nav,
                    '<div class="grid">', fastqc_section(fq, opt$sample_id, zip_guess),
                    viral_length_card(vstats), '</div>'))

  write_html(file.path(multi_dir, "all.html"), "All sites",
             paste0('<h1>All integration sites</h1>', nav,
                    '<div class="hint">Heaviest page in the report. ',
                    'Use per-chromosome pages for low-RAM machines.</div>',
                    render_table(dt, "tbl-all")))

  if (!is.null(mdt)) {
    write_html(file.path(multi_dir, "reads.html"), "Reads",
               paste0('<h1>Read-level detail</h1>', nav,
                      render_table(mdt, "tbl-reads")))
  }

  for (c in chrs) {
    sub <- dt[chrom == c]
    write_html(file.path(multi_dir, sprintf("chrom_%s.html", c)),
               sprintf("%s %s", opt$sample_id, c),
               paste0('<h1>', .esc(opt$sample_id), ' on ', .esc(c),
                      '</h1>', nav, render_table(sub, paste0("tbl-", c))))
  }

  message("[3.Create_Sample_HTML.R] wrote multi-page report into ", multi_dir,
          " (", length(chrs) + 4 + as.integer(!is.null(mdt)), " pages)")
}

mode <- tolower(opt$html_mode)
if (mode %in% c("single", "both")) { write_single() }
if (mode %in% c("multi",  "both")) { write_multi() }