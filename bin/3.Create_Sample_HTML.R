#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 3.Create_Sample_HTML.R   (base-R edition, no htmltools)
#
# Per-sample HTML report. Uses ONLY packages already present in the project
# SIFs (R_Genomics_v1_blast.def, Viral_Genomics_v5.def):
#       data.table, optparse  (and base R for HTML emission)
#
# Modes (--html_mode):
#   single  : one HTML file with external image references and lazy table
#   multi   : index + per-chromosome pages + circos + table pages (low-RAM)
#   both    : both (default; matches the original request)
#
# Inputs:
#   --integrations  CSV/TSV with at least one of each pair below:
#                       chrom        | chromosome
#                       position     | integration_site
#                       sample_id    | sample
#                       clonal_id    | clone_id
#                       read_support | reads_at_site
#                   plus optional: timepoint, viral_orientation, gene_name,
#                   in_repeat_t2t, in_repeat_hg38
#   --circos_png    PNG produced by make_circos.R (referenced externally,
#                   NOT base64-embedded). Optional — placeholder shown if
#                   missing.
#   --sample_id     Sample identifier (used in titles + file names)
#   --outdir        Where to write HTML(s)
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--integrations", type = "character"),
  make_option("--circos_png",   type = "character", default = NA),
  make_option("--sample_id",    type = "character", default = "sample"),
  make_option("--html_mode",    type = "character", default = "both"),
  make_option("--outdir",       type = "character", default = ".")
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
  clonal_id    = c("clonal_id",    "clone_id",         "Clone_ID",   "ClonalID"),
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

.css <- '
:root { --fg:#222; --muted:#666; --bg:#fff; --accent:#1f6feb; --row:#f6f8fa; }
* { box-sizing: border-box; }
body { font-family:-apple-system,Segoe UI,Helvetica,Arial,sans-serif;
       color:var(--fg); background:var(--bg); margin:0 auto; max-width:1100px;
       padding:16px 24px; line-height:1.45; }
h1 { font-size:22px; margin:4px 0 12px; }
h2 { font-size:16px; margin:18px 0 8px; color:var(--accent); }
.meta { color:var(--muted); font-size:13px; margin-bottom:16px; }
.grid { display:grid; grid-template-columns:1fr 1fr; gap:12px; }
.card { border:1px solid #e1e4e8; border-radius:6px; padding:12px; background:#fff; }
table { border-collapse:collapse; width:100%; font-size:12.5px; }
th, td { padding:5px 8px; text-align:left; border-bottom:1px solid #eee;
         white-space:nowrap; }
th { background:var(--row); position:sticky; top:0; }
tr:hover td { background:#fafbfc; }
.nav { margin:12px 0 16px; font-size:13px; }
.nav a { color:var(--accent); text-decoration:none; margin-right:10px; }
img.circos { max-width:100%; height:auto; border:1px solid #e1e4e8; }
.hint { color:var(--muted); font-size:12px; margin-top:6px; }
.scroll { max-height:60vh; overflow:auto; border:1px solid #e1e4e8; border-radius:6px; }
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
    '<!doctype html><html><head><meta charset="utf-8"><title>',
    .esc(title), '</title><style>', .css, '</style><script>', .js,
    '</script></head><body>', body_html, '</body></html>'
  )
  writeLines(doc, path, useBytes = TRUE)
}

# ---- Render a (potentially large) table as a string -----------------------
# Choose a focused set of columns for the inset table so the user actually
# sees the integration metadata, not 38 columns including the full viral
# nucleotide sequence. Falls back to all columns if the canonical ones aren't
# present.
.table_columns <- function(dt) {
  preferred <- c("sample_id", "clonal_id", "chrom", "position",
                 "viral_orientation", "gene_name", "gene_id",
                 "read_support", "integration_site",
                 "MATCH_TYPE", "IPDA_INTACT", "IPDA_V2_INTACT",
                 "in_repeat_t2t", "in_repeat_hg38",
                 "repeat_name_t2t", "repeat_name_hg38")
  cols <- intersect(preferred, colnames(dt))
  if (length(cols) == 0) cols <- colnames(dt)
  cols
}

render_table <- function(dt_chunk, table_id) {
  cols <- .table_columns(dt_chunk)
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
  paste0(
    '<div><input id="', table_id,
    '-q" type="search" placeholder="filter…  (try: chr1, AluY, gene name)" ',
    'oninput="filt(\'', table_id, '\')" ',
    'style="width:100%;padding:6px 8px;margin-bottom:6px;"/></div>',
    '<div class="scroll"><table id="', table_id, '">', thead, tbody, '</table></div>',
    '<div class="hint">', nrow(dt_chunk),
    ' rows. Showing ', length(cols), '/', ncol(dt_chunk),
    ' columns. Filter is client-side; full data is in the underlying CSV.</div>'
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

circos_card <- function(circos_png, src_override = NULL, embed = TRUE,
                        embed_max_kb = 500) {
  # If caller specified an explicit src (multi-page mode passes "circos.png"),
  # honour it as a relative file reference and skip base64 logic.
  if (!is.null(src_override)) {
    return(paste0('<div class="card"><h2>Circos</h2>',
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
      '<div class="card"><h2>Circos</h2>',
      '<div class="hint">', .esc(msg),
      '. Make sure make_circos.R ran before this report.</div></div>'
    ))
  }

  size_kb <- file.info(circos_png)$size / 1024

  # Default: embed as base64 data: URL when small. This makes the single-page
  # HTML self-contained — no "is the PNG in the published dir" failure mode.
  # Use whichever base64 implementation happens to be available; fall through
  # to a pure-base-R encoder so we have zero hard dependencies.
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
  if (embed && size_kb < embed_max_kb) {
    b64 <- encode_base64(circos_png)
    src <- paste0("data:image/png;base64,", b64)
    return(paste0(
      '<div class="card"><h2>Circos</h2>',
      '<img class="circos" src="', src,
      '" alt="Per-sample circos plot"/>',
      '<div class="hint">Embedded inline (', round(size_kb, 1),
      ' KB) — HTML is self-contained.</div></div>'
    ))
  }

  # Big PNG: fall back to external reference (and rely on the file_copy
  # in write_single() to put it next to the HTML).
  paste0('<div class="card"><h2>Circos</h2>',
         '<img class="circos" src="', .esc(basename(circos_png)),
         '" alt="Per-sample circos plot"/></div>')
}

# ---- Single-page mode ------------------------------------------------------
write_single <- function() {
  # Stage circos image into outdir so the relative <img src=...> resolves.
  if (!is.na(opt$circos_png) && file.exists(opt$circos_png)) {
    dest <- file.path(opt$outdir, basename(opt$circos_png))
    if (normalizePath(opt$circos_png, mustWork = TRUE) !=
        suppressWarnings(normalizePath(dest, mustWork = FALSE))) {
      file.copy(opt$circos_png, dest, overwrite = TRUE)
    }
  }
  body <- paste0(
    '<h1>Viral integration report — ', .esc(opt$sample_id), '</h1>',
    '<div class="meta">Generated ',
       format(Sys.time(), "%Y-%m-%d %H:%M"), '</div>',
    '<div class="grid">', summary_card(dt), circos_card(opt$circos_png), '</div>',
    '<h2>Integration sites</h2>', render_table(dt, "tbl-main")
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
    '<a href="all.html">All sites</a>',
    paste0('<a href="chrom_', .esc(chrs), '.html">', .esc(chrs), '</a>',
           collapse = ""),
    '</div>'
  )

  body_idx <- paste0(
    '<h1>Viral integration report — ', .esc(opt$sample_id), '</h1>',
    '<div class="meta">Multi-page report. Generated ',
       format(Sys.time(), "%Y-%m-%d %H:%M"), '</div>', nav,
    '<div class="grid">', summary_card(dt),
       circos_card(NA, src_override = "circos.png"), '</div>'
  )
  write_html(file.path(multi_dir, "index.html"),
             sprintf("%s index", opt$sample_id), body_idx)

  write_html(file.path(multi_dir, "circos.html"), "Circos",
             paste0('<h1>Circos plot</h1>', nav,
                    '<div class="card"><img class="circos" src="circos.png"/></div>'))


  write_html(file.path(multi_dir, "all.html"), "All sites",
             paste0('<h1>All integration sites</h1>', nav,
                    '<div class="hint">Heaviest page in the report. ',
                    'Use per-chromosome pages for low-RAM machines.</div>',
                    render_table(dt, "tbl-all")))

  for (c in chrs) {
    sub <- dt[chrom == c]
    write_html(file.path(multi_dir, sprintf("chrom_%s.html", c)),
               sprintf("%s %s", opt$sample_id, c),
               paste0('<h1>', .esc(opt$sample_id), ' on ', .esc(c),
                      '</h1>', nav, render_table(sub, paste0("tbl-", c))))
  }

  message("[3.Create_Sample_HTML.R] wrote multi-page report into ", multi_dir,
          " (", length(chrs) + 3, " pages)")
}

mode <- tolower(opt$html_mode)
if (mode %in% c("single", "both")) { write_single() }
if (mode %in% c("multi",  "both")) { write_multi() }
