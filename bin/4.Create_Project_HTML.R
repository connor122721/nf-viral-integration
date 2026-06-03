#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# 4.Create_Project_HTML.R   (base-R edition, no htmltools)
#
# Project-wide rollup over all samples. Uses ONLY packages already in the
# project SIFs:  data.table, optparse  (and base R for HTML emission).
#
# Inputs:
#   --integrations    Concatenated integrations TSV across samples (long,
#                     after assign_clonal_ids.R has added clonal_id).
#   --persistence     Wide clone × timepoint matrix (clonal_persistence_wide.tsv).
#   --summary         Per-clone summary (clonal_summary.tsv).
#   --circos_png      Project-wide circos PNG.
#   --sample_reports  Path to the directory holding per-sample HTML reports
#                     (used to make clickable links). Optional.
#   --html_mode       single | multi | both (default both)
#   --outdir          Output directory.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--integrations",   type = "character"),
  make_option("--persistence",    type = "character"),
  make_option("--summary",        type = "character"),
  make_option("--circos_png",     type = "character", default = NA),
  make_option("--sample_reports", type = "character", default = NA),
  make_option("--html_mode",      type = "character", default = "both"),
  make_option("--outdir",         type = "character", default = ".")
)))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dt   <- fread(opt$integrations)
pers <- if (!is.null(opt$persistence) && file.exists(opt$persistence))
          fread(opt$persistence) else NULL
summ <- if (!is.null(opt$summary)     && file.exists(opt$summary))
          fread(opt$summary)     else NULL

# ---- HTML helpers (same as sample script; kept inline to make this file
#      independent of the sample script if dropped in alone) -----------------
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
       color:var(--fg); background:var(--bg); margin:0 auto; max-width:1200px;
       padding:16px 24px; line-height:1.45; }
h1 { font-size:24px; margin:4px 0 12px; }
h2 { font-size:17px; margin:18px 0 8px; color:var(--accent); }
.meta { color:var(--muted); font-size:13px; margin-bottom:16px; }
.grid { display:grid; grid-template-columns:2fr 3fr; gap:14px; }
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

render_table <- function(dt_chunk, table_id) {
  cols <- colnames(dt_chunk)
  thead <- paste0("<thead><tr>",
                  paste0("<th>", .esc(cols), "</th>", collapse = ""),
                  "</tr></thead>")
  cell_mat <- vapply(cols, function(c) .esc(dt_chunk[[c]]),
                     character(nrow(dt_chunk)))
  if (is.null(dim(cell_mat))) cell_mat <- matrix(cell_mat, nrow = 1)
  rows <- apply(cell_mat, 1, function(r)
    paste0("<tr>", paste0("<td>", r, "</td>", collapse = ""), "</tr>"))
  tbody <- paste0("<tbody>", paste(rows, collapse = ""), "</tbody>")
  paste0(
    '<div><input id="', table_id,
    '-q" type="search" placeholder="filter…" oninput="filt(\'',
    table_id, '\')" style="width:100%;padding:6px 8px;margin-bottom:6px;"/></div>',
    '<div class="scroll"><table id="', table_id, '">', thead, tbody, '</table></div>',
    '<div class="hint">', nrow(dt_chunk), ' rows.</div>'
  )
}

# ---- Cards -----------------------------------------------------------------
project_summary <- function(dt) {
  paste0(
    '<div class="card"><h2>Project summary</h2>',
    '<div>Samples: ', length(unique(dt$sample_id)), '</div>',
    '<div>Total integrations: ', nrow(dt), '</div>',
    if ("clonal_id" %in% colnames(dt))
      paste0('<div>Unique clones: ', length(unique(dt$clonal_id)), '</div>') else "",
    if ("timepoint" %in% colnames(dt))
      paste0('<div>Timepoints: ',
             .esc(paste(sort(unique(stats::na.omit(dt$timepoint))),
                        collapse = ", ")),
             '</div>') else "",
    '</div>'
  )
}

sample_links <- function(dt, sample_reports_dir) {
  per_sample <- dt[, .(n_int   = .N,
                       n_clone = if ("clonal_id" %in% colnames(dt))
                                   uniqueN(clonal_id) else NA_integer_),
                   by = sample_id]
  setorder(per_sample, sample_id)
  rows <- vapply(seq_len(nrow(per_sample)), function(i) {
    sid  <- per_sample$sample_id[i]
    href <- if (!is.na(sample_reports_dir))
              file.path(basename(sample_reports_dir),
                        sprintf("%s.report.html", sid)) else "#"
    paste0('<tr><td><a href="', .esc(href), '">', .esc(sid),
           '</a></td><td>', per_sample$n_int[i], '</td><td>',
           if (is.na(per_sample$n_clone[i])) "n/a" else per_sample$n_clone[i],
           '</td></tr>')
  }, character(1))
  paste0(
    '<div class="card"><h2>Per-sample reports</h2>',
    '<table><thead><tr><th>sample_id</th><th>integrations</th>',
    '<th>clones</th></tr></thead><tbody>',
    paste(rows, collapse = ""),
    '</tbody></table></div>'
  )
}

circos_card <- function(circos_png, dest_dir, src_override = NULL) {
  if (is.null(src_override)) {
    if (is.na(circos_png) || !file.exists(circos_png)) {
      return('<div class="card"><h2>Project circos</h2><div class="hint">No circos PNG provided.</div></div>')
    }
    src <- "project_circos.png"
    file.copy(circos_png, file.path(dest_dir, src), overwrite = TRUE)
  } else {
    src <- src_override
  }
  paste0('<div class="card"><h2>Project circos</h2>',
         '<img class="circos" src="', .esc(src),
         '" alt="Project-wide circos"/></div>')
}

# ---- Single page -----------------------------------------------------------
write_single <- function() {
  body <- paste0(
    '<h1>HIV viral integration — project report</h1>',
    '<div class="meta">Generated ',
       format(Sys.time(), "%Y-%m-%d %H:%M"), '</div>',
    '<div class="grid">',
       project_summary(dt),
       circos_card(opt$circos_png, opt$outdir),
    '</div>',
    sample_links(dt, opt$sample_reports),
    if (!is.null(summ))
      paste0('<h2>Top persistent clones (by # timepoints, then samples)</h2>',
             render_table(head(summ, 100), "tbl-summary")) else "",
    if (!is.null(pers))
      paste0('<h2>Clone × timepoint persistence (full)</h2>',
             render_table(pers, "tbl-persistence")) else ""
  )
  out <- file.path(opt$outdir, "project.report.html")
  write_html(out, "Project report", body)
  message("[4.Create_Project_HTML.R] wrote single-page: ", out, " (",
          round(file.info(out)$size / 1024, 1), " KB)")
}

# ---- Multi page ------------------------------------------------------------
write_multi <- function() {
  multi_dir <- file.path(opt$outdir, "project_report")
  dir.create(multi_dir, showWarnings = FALSE, recursive = TRUE)
  if (!is.na(opt$circos_png) && file.exists(opt$circos_png))
    file.copy(opt$circos_png, file.path(multi_dir, "project_circos.png"),
              overwrite = TRUE)

  nav <- paste0('<div class="nav">',
                '<a href="index.html">Index</a>',
                '<a href="circos.html">Circos</a>',
                '<a href="samples.html">Samples</a>',
                '<a href="persistence.html">Persistence</a>',
                '<a href="summary.html">Top clones</a>',
                '</div>')

  write_html(file.path(multi_dir, "index.html"), "Project index",
             paste0('<h1>HIV integration — project</h1>', nav,
                    '<div class="grid">', project_summary(dt),
                    circos_card(NA, multi_dir, src_override = "project_circos.png"),
                    '</div>'))

  write_html(file.path(multi_dir, "circos.html"), "Circos",
             paste0('<h1>Project circos</h1>', nav,
                    '<div class="card"><img class="circos" src="project_circos.png"/></div>'))

  write_html(file.path(multi_dir, "samples.html"), "Samples",
             paste0('<h1>Samples</h1>', nav,
                    sample_links(dt, opt$sample_reports)))

  if (!is.null(summ))
    write_html(file.path(multi_dir, "summary.html"), "Top clones",
               paste0('<h1>Most persistent clones</h1>', nav,
                      render_table(summ, "tbl-summary")))

  if (!is.null(pers))
    write_html(file.path(multi_dir, "persistence.html"), "Persistence",
               paste0('<h1>Clone × timepoint persistence</h1>', nav,
                      render_table(pers, "tbl-pers")))

  message("[4.Create_Project_HTML.R] wrote multi-page report into ", multi_dir)
}

mode <- tolower(opt$html_mode)
if (mode %in% c("single", "both")) write_single()
if (mode %in% c("multi",  "both")) write_multi()
