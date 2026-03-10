#!/usr/bin/env Rscript
#
# Viral Integration HTML Report Generator
# Generates a self-contained, interactive HTML report for one or more samples.
#
# Usage:
#   Rscript create_report.R <results_dir> <output_prefix> [run_name]
#
#   results_dir    - Directory containing *_annotated.csv files (all samples)
#   output_prefix  - Prefix (and path) for the output HTML file
#   run_name       - Optional display name shown in the report header

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# ---------------------------------------------------------------------------
# 1. Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript create_report.R <results_dir> <output_prefix> [run_name]\n\n")
  cat("  results_dir    - Directory with *_annotated.csv files\n")
  cat("  output_prefix  - Prefix for HTML output (e.g. 'my_run')\n")
  cat("  run_name       - Display name for the report [optional]\n")
  quit(status = 1)
}

results_dir   <- args[1]
output_prefix <- args[2]
run_name      <- if (length(args) >= 3) args[3] else basename(results_dir)
report_date   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

cat("\n")
cat(strrep("=", 70), "\n")
cat("VIRAL INTEGRATION HTML REPORT GENERATOR\n")
cat(strrep("=", 70), "\n\n")

# ---------------------------------------------------------------------------
# 2. Discover and load all annotated CSV files
# ---------------------------------------------------------------------------
csv_files <- list.files(results_dir, pattern = "_annotated\\.csv$",
                        recursive = TRUE, full.names = TRUE)

if (length(csv_files) == 0) {
  cat("WARNING: No *_annotated.csv files found in", results_dir, "\n")
  cat("Searching for any *.csv files...\n")
  csv_files <- list.files(results_dir, pattern = "\\.csv$",
                          recursive = TRUE, full.names = TRUE)
}

if (length(csv_files) == 0) {
  stop("No CSV result files found in: ", results_dir)
}

cat(sprintf("Found %d result file(s):\n", length(csv_files)))
for (f in csv_files) cat("  -", f, "\n")

# Load and combine
all_data <- lapply(csv_files, function(f) {
  tryCatch({
    d <- read_csv(f, show_col_types = FALSE, guess_max = 5000)
    # Ensure a 'sample' column exists
    if (!"sample" %in% colnames(d)) {
      d$sample <- sub("_annotated\\.csv$", "", basename(f))
    }
    d
  }, error = function(e) {
    cat(sprintf("  WARNING: Could not read %s: %s\n", basename(f), conditionMessage(e)))
    NULL
  })
})
all_data <- bind_rows(Filter(Negate(is.null), all_data))

n_samples <- n_distinct(all_data$sample)
cat(sprintf("\nLoaded %d rows across %d sample(s)\n", nrow(all_data), n_samples))

# ---------------------------------------------------------------------------
# 3. Helper: encode a ggplot as an inline base64 PNG
# ---------------------------------------------------------------------------
plot_to_base64 <- function(gg_obj, width = 7, height = 5, res = 150) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  png(tmp, width = width, height = height, units = "in", res = res,
      bg = "transparent")
  print(gg_obj)
  dev.off()
  b64 <- base64enc::base64encode(tmp)
  sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;height:auto;">', b64)
}

# Fallback when base64enc is not available: write PNG to file and use file://
plot_to_img_tag <- function(gg_obj, fname, width = 7, height = 5, res = 150) {
  png(fname, width = width, height = height, units = "in", res = res,
      bg = "transparent")
  print(gg_obj)
  dev.off()
  sprintf('<img src="%s" style="max-width:100%%;height:auto;">', basename(fname))
}

# Try to load base64enc; fall back to file-based images
use_base64 <- requireNamespace("base64enc", quietly = TRUE)
img_dir <- paste0(output_prefix, "_report_files")
if (!use_base64) dir.create(img_dir, showWarnings = FALSE)

make_img <- function(gg_obj, label, width = 7, height = 5) {
  if (use_base64) {
    plot_to_base64(gg_obj, width = width, height = height)
  } else {
    fname <- file.path(img_dir, paste0(label, ".png"))
    plot_to_img_tag(gg_obj, fname, width = width, height = height)
  }
}

# ---------------------------------------------------------------------------
# 4. Build summary statistics
# ---------------------------------------------------------------------------

## Per-sample read counts
sample_counts <- all_data %>%
  count(sample, name = "n_reads") %>%
  arrange(desc(n_reads))

## Orientation summary
if ("viral_orientation" %in% colnames(all_data)) {
  ori_summary <- all_data %>%
    filter(!is.na(viral_orientation)) %>%
    count(sample, viral_orientation) %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()
} else {
  ori_summary <- tibble()
}

## Chromosome distribution
if ("chromosome" %in% colnames(all_data)) {
  chr_summary <- all_data %>%
    filter(!is.na(chromosome)) %>%
    count(chromosome, sort = TRUE) %>%
    head(20)
} else {
  chr_summary <- tibble()
}

## Insert length distribution
if ("INSERT_LEN" %in% colnames(all_data)) {
  len_stats <- all_data %>%
    filter(!is.na(INSERT_LEN), INSERT_LEN > 0) %>%
    group_by(sample) %>%
    summarise(
      n           = n(),
      mean_len    = round(mean(INSERT_LEN), 0),
      median_len  = round(median(INSERT_LEN), 0),
      min_len     = min(INSERT_LEN),
      max_len     = max(INSERT_LEN),
      .groups     = "drop"
    )
} else {
  len_stats <- tibble()
}

## Percent identity
if ("percent_identity" %in% colnames(all_data)) {
  pid_stats <- all_data %>%
    filter(!is.na(percent_identity)) %>%
    group_by(sample) %>%
    summarise(
      mean_pid   = round(mean(percent_identity), 1),
      median_pid = round(median(percent_identity), 1),
      min_pid    = round(min(percent_identity), 1),
      max_pid    = round(max(percent_identity), 1),
      .groups    = "drop"
    )
} else {
  pid_stats <- tibble()
}

## Gene annotation summary
if ("gene_name" %in% colnames(all_data)) {
  gene_summary <- all_data %>%
    filter(!is.na(gene_name)) %>%
    separate_rows(gene_name, sep = ";") %>%
    count(gene_name, sort = TRUE) %>%
    head(20)
} else {
  gene_summary <- tibble()
}

# ---------------------------------------------------------------------------
# 5. Build plots
# ---------------------------------------------------------------------------
theme_report <- function(dark = FALSE) {
  bg    <- if (dark) "#1e1e2e" else "#ffffff"
  fg    <- if (dark) "#cdd6f4" else "#333333"
  panel <- if (dark) "#313244" else "#f9f9f9"
  grid  <- if (dark) "#45475a" else "#e0e0e0"
  theme_minimal(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = bg,    color = NA),
      panel.background = element_rect(fill = panel, color = NA),
      text             = element_text(color = fg),
      axis.text        = element_text(color = fg),
      axis.title       = element_text(color = fg),
      plot.title       = element_text(color = fg, face = "bold"),
      plot.subtitle    = element_text(color = fg),
      legend.background = element_rect(fill = bg,   color = NA),
      legend.text      = element_text(color = fg),
      legend.title     = element_text(color = fg),
      panel.grid.major = element_line(color = grid),
      panel.grid.minor = element_line(color = grid, linetype = "dashed"),
      strip.text       = element_text(color = fg, face = "bold")
    )
}

pal_samples <- hue_pal()(max(n_samples, 1))

## Plot 1: reads per sample
p_reads <- ggplot(sample_counts, aes(x = reorder(sample, n_reads), y = n_reads,
                                      fill = sample)) +
  geom_col(show.legend = FALSE, width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = colorRampPalette(c("#74c0fc", "#4dabf7", "#339af0"))(nrow(sample_counts))) +
  scale_y_continuous(labels = label_comma()) +
  labs(title = "Reads per Sample", x = NULL, y = "Number of Reads") +
  theme_report()

## Plot 2: Orientation breakdown
if (nrow(ori_summary) > 0) {
  p_ori <- ggplot(ori_summary,
                  aes(x = sample, y = n, fill = viral_orientation)) +
    geom_col(position = "fill", width = 0.7) +
    scale_y_continuous(labels = label_percent()) +
    scale_fill_manual(values = c("5'" = "#51cf66", "3'" = "#ff6b6b"),
                      na.value = "#aaa") +
    labs(title = "Viral Orientation (Sense vs. Antisense)",
         x = NULL, y = "Proportion", fill = "Orientation") +
    theme_report() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p_ori <- NULL
}

## Plot 3: Chromosome distribution
if (nrow(chr_summary) > 0) {
  p_chr <- ggplot(chr_summary, aes(x = reorder(chromosome, n), y = n)) +
    geom_col(fill = "#74c0fc", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Top Integration Chromosomes",
         x = "Chromosome", y = "Number of Reads") +
    theme_report()
} else {
  p_chr <- NULL
}

## Plot 4: Insert length distribution
if ("INSERT_LEN" %in% colnames(all_data)) {
  len_data <- all_data %>% filter(!is.na(INSERT_LEN), INSERT_LEN > 0)
  if (nrow(len_data) > 0) {
    p_len <- ggplot(len_data, aes(x = INSERT_LEN, fill = sample)) +
      geom_histogram(bins = 60, alpha = 0.8, position = "identity") +
      scale_x_continuous(labels = label_comma()) +
      scale_y_continuous(labels = label_comma()) +
      scale_fill_manual(values = colorRampPalette(c("#74c0fc", "#f783ac", "#69db7c"))(n_samples)) +
      facet_wrap(~sample, scales = "free_y", ncol = min(n_samples, 3)) +
      labs(title = "Insert (Viral Sequence) Length Distribution",
           x = "Insert Length (bp)", y = "Count", fill = "Sample") +
      theme_report()
  } else {
    p_len <- NULL
  }
} else {
  p_len <- NULL
}

## Plot 5: Percent identity distribution
if ("percent_identity" %in% colnames(all_data)) {
  pid_data <- all_data %>% filter(!is.na(percent_identity))
  if (nrow(pid_data) > 0) {
    p_pid <- ggplot(pid_data, aes(x = percent_identity, fill = sample)) +
      geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
      scale_fill_manual(values = colorRampPalette(c("#74c0fc", "#f783ac", "#69db7c"))(n_samples)) +
      facet_wrap(~sample, scales = "free_y", ncol = min(n_samples, 3)) +
      labs(title = "Viral Sequence % Identity to Reference",
           x = "% Identity", y = "Count", fill = "Sample") +
      theme_report()
  } else {
    p_pid <- NULL
  }
} else {
  p_pid <- NULL
}

## Plot 6: Top genes
if (nrow(gene_summary) > 0) {
  p_genes <- ggplot(gene_summary, aes(x = reorder(gene_name, n), y = n)) +
    geom_col(fill = "#f783ac", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Top Integration-Bearing Genes",
         x = "Gene", y = "Reads with Integration") +
    theme_report()
} else {
  p_genes <- NULL
}

# ---------------------------------------------------------------------------
# 6. Helper: render a data frame as an HTML table
# ---------------------------------------------------------------------------
df_to_html_table <- function(df, id = NULL, caption = NULL) {
  if (nrow(df) == 0) return("<p><em>No data available.</em></p>")
  id_attr  <- if (!is.null(id))      sprintf(' id="%s"', id) else ""
  cap_html <- if (!is.null(caption))
                sprintf('<caption>%s</caption>', caption) else ""
  
  header <- paste(sprintf("<th>%s</th>", colnames(df)), collapse = "")
  rows   <- apply(df, 1, function(r) {
    cells <- paste(sprintf("<td>%s</td>", r), collapse = "")
    sprintf("<tr>%s</tr>", cells)
  })
  sprintf('<div class="table-wrap"><table%s>%s<thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
          id_attr, cap_html, header, paste(rows, collapse = "\n"))
}

# ---------------------------------------------------------------------------
# 7. Build section HTML
# ---------------------------------------------------------------------------

## Plots section
build_plot_section <- function(title, img_tag) {
  if (is.null(img_tag)) return("")
  sprintf('<div class="card">
  <h2>%s</h2>
  <div class="plot-container">%s</div>
</div>', title, img_tag)
}

reads_img <- make_img(p_reads, "reads_per_sample", width = 8, height = 4)
ori_img   <- if (!is.null(p_ori))   make_img(p_ori,   "orientation",    width = 8, height = 4) else NULL
chr_img   <- if (!is.null(p_chr))   make_img(p_chr,   "chromosome",     width = 8, height = 5) else NULL
len_img   <- if (!is.null(p_len))   make_img(p_len,   "insert_length",  width = 10, height = 4) else NULL
pid_img   <- if (!is.null(p_pid))   make_img(p_pid,   "percent_id",     width = 10, height = 4) else NULL
gene_img  <- if (!is.null(p_genes)) make_img(p_genes, "top_genes",      width = 8, height = 5) else NULL

## QC stats table
qc_table_df <- sample_counts
if (nrow(len_stats) > 0)
  qc_table_df <- left_join(qc_table_df, len_stats %>% dplyr::select(sample, mean_len, median_len),
                           by = "sample")
if (nrow(pid_stats) > 0)
  qc_table_df <- left_join(qc_table_df, pid_stats %>% dplyr::select(sample, mean_pid),
                           by = "sample")

## Integration sites table (first 200 rows for display)
integration_table_df <- all_data %>%
  dplyr::select(any_of(c("sample", "READ", "chromosome", "integration_site",
                          "viral_orientation", "viral_strand",
                          "INSERT_LEN", "percent_identity",
                          "gene_name", "gene_id"))) %>%
  head(200)

# ---------------------------------------------------------------------------
# 8. Assemble HTML
# ---------------------------------------------------------------------------
html_file <- paste0(output_prefix, "_report.html")

html_css <- '
<style>
/* ===== CSS Variables for theming ===== */
:root {
  --bg:       #ffffff;
  --bg2:      #f4f6fb;
  --card:     #ffffff;
  --border:   #d0d7de;
  --fg:       #1f2328;
  --fg2:      #57606a;
  --accent:   #0969da;
  --accent2:  #8250df;
  --success:  #1a7f37;
  --header:   #24292f;
  --header-fg:#ffffff;
  --shadow:   rgba(0,0,0,0.06);
  --badge-bg: #ddf4ff;
  --badge-fg: #0550ae;
  --th-bg:    #f6f8fa;
  --tr-hover: #f6f8fa;
}
[data-theme="dark"] {
  --bg:       #0d1117;
  --bg2:      #161b22;
  --card:     #1c2128;
  --border:   #30363d;
  --fg:       #e6edf3;
  --fg2:      #848d97;
  --accent:   #58a6ff;
  --accent2:  #d2a8ff;
  --success:  #3fb950;
  --header:   #161b22;
  --header-fg:#e6edf3;
  --shadow:   rgba(0,0,0,0.3);
  --badge-bg: #1f3a59;
  --badge-fg: #79c0ff;
  --th-bg:    #21262d;
  --tr-hover: #21262d;
}

/* ===== Base ===== */
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
  background: var(--bg2);
  color: var(--fg);
  line-height: 1.6;
  transition: background 0.3s, color 0.3s;
}
a { color: var(--accent); text-decoration: none; }
a:hover { text-decoration: underline; }

/* ===== Header ===== */
header {
  background: linear-gradient(135deg, #0969da 0%, #8250df 100%);
  color: #ffffff;
  padding: 2rem 2.5rem 1.5rem;
}
header h1 { font-size: 1.9rem; font-weight: 700; margin-bottom: 0.2rem; }
header .subtitle { opacity: 0.85; font-size: 0.95rem; }
header .meta { margin-top: 0.7rem; font-size: 0.8rem; opacity: 0.7; }

/* ===== Theme toggle ===== */
.theme-toggle {
  position: fixed; top: 1rem; right: 1rem; z-index: 999;
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 2rem;
  padding: 0.35rem 0.9rem;
  cursor: pointer;
  font-size: 0.85rem;
  color: var(--fg);
  box-shadow: 0 2px 8px var(--shadow);
  transition: background 0.3s, color 0.3s;
}
.theme-toggle:hover { background: var(--bg2); }

/* ===== Layout ===== */
.container { max-width: 1200px; margin: 0 auto; padding: 2rem 1.5rem; }

/* ===== Cards ===== */
.card {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 10px;
  padding: 1.6rem 2rem;
  margin-bottom: 1.8rem;
  box-shadow: 0 2px 8px var(--shadow);
  transition: background 0.3s, border-color 0.3s;
}
.card h2 {
  font-size: 1.15rem;
  font-weight: 700;
  margin-bottom: 1rem;
  color: var(--fg);
  border-bottom: 2px solid var(--border);
  padding-bottom: 0.5rem;
}

/* ===== Stats grid ===== */
.stats-grid {
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(180px, 1fr));
  gap: 1rem;
  margin-bottom: 0.5rem;
}
.stat-box {
  background: var(--bg2);
  border: 1px solid var(--border);
  border-radius: 8px;
  padding: 1rem 1.2rem;
  text-align: center;
}
.stat-box .stat-value {
  font-size: 1.7rem;
  font-weight: 700;
  color: var(--accent);
}
.stat-box .stat-label {
  font-size: 0.75rem;
  color: var(--fg2);
  margin-top: 0.2rem;
  text-transform: uppercase;
  letter-spacing: 0.05em;
}

/* ===== Tables ===== */
.table-wrap { overflow-x: auto; }
table {
  width: 100%;
  border-collapse: collapse;
  font-size: 0.82rem;
}
thead { background: var(--th-bg); }
th {
  padding: 0.55rem 0.75rem;
  text-align: left;
  font-weight: 600;
  border-bottom: 2px solid var(--border);
  color: var(--fg);
  white-space: nowrap;
}
td {
  padding: 0.45rem 0.75rem;
  border-bottom: 1px solid var(--border);
  color: var(--fg2);
  max-width: 260px;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}
tr:hover td { background: var(--tr-hover); }

/* ===== Badges ===== */
.badge {
  display: inline-block;
  background: var(--badge-bg);
  color: var(--badge-fg);
  border-radius: 1em;
  padding: 0.15em 0.6em;
  font-size: 0.75rem;
  font-weight: 600;
}

/* ===== Plot container ===== */
.plot-container { text-align: center; }
.plot-container img { max-width: 100%; border-radius: 6px; }

/* ===== Section nav ===== */
.section-nav {
  display: flex; flex-wrap: wrap; gap: 0.5rem;
  margin-bottom: 1.5rem;
}
.section-nav a {
  background: var(--card);
  border: 1px solid var(--border);
  border-radius: 1em;
  padding: 0.3em 0.8em;
  font-size: 0.82rem;
  color: var(--accent);
}
.section-nav a:hover { background: var(--bg2); }

/* ===== Footer ===== */
footer {
  text-align: center;
  font-size: 0.78rem;
  color: var(--fg2);
  padding: 2rem 0 1rem;
  border-top: 1px solid var(--border);
}
</style>
'

html_js <- '
<script>
(function() {
  var saved = localStorage.getItem("theme") || "light";
  document.documentElement.setAttribute("data-theme", saved);

  function toggleTheme() {
    var current = document.documentElement.getAttribute("data-theme");
    var next    = current === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    localStorage.setItem("theme", next);
    document.getElementById("toggle-btn").textContent =
      next === "dark" ? "☀️ Light mode" : "🌙 Dark mode";
  }

  document.addEventListener("DOMContentLoaded", function() {
    var btn = document.getElementById("toggle-btn");
    btn.textContent = saved === "dark" ? "☀️ Light mode" : "🌙 Dark mode";
    btn.addEventListener("click", toggleTheme);
  });
})();
</script>
'

# Summary stats boxes
total_reads <- sum(sample_counts$n_reads)
unique_sites <- if ("integration_site" %in% colnames(all_data))
  n_distinct(all_data$integration_site, na.rm = TRUE) else "—"
pct_oriented <- if ("viral_orientation" %in% colnames(all_data))
  sprintf("%.0f%%", 100 * mean(!is.na(all_data$viral_orientation))) else "—"
pct_annotated <- if ("gene_name" %in% colnames(all_data))
  sprintf("%.0f%%", 100 * mean(!is.na(all_data$gene_name))) else "—"

stats_html <- sprintf('
<div class="stats-grid">
  <div class="stat-box"><div class="stat-value">%s</div><div class="stat-label">Total Reads</div></div>
  <div class="stat-box"><div class="stat-value">%s</div><div class="stat-label">Samples</div></div>
  <div class="stat-box"><div class="stat-value">%s</div><div class="stat-label">Unique Integration Sites</div></div>
  <div class="stat-box"><div class="stat-value">%s</div><div class="stat-label">Orientation Called</div></div>
  <div class="stat-box"><div class="stat-value">%s</div><div class="stat-label">Reads in Genes</div></div>
</div>',
  format(total_reads, big.mark = ","),
  n_samples,
  if (is.numeric(unique_sites)) format(unique_sites, big.mark = ",") else unique_sites,
  pct_oriented,
  pct_annotated
)

# Build page sections
sections_html <- paste0(
  # Overview
  sprintf('<div class="card" id="overview">
<h2>📊 Run Overview</h2>
%s
%s
</div>', stats_html, df_to_html_table(qc_table_df, id = "qc-table", caption = "Per-sample QC Summary")),

  # Reads per sample plot
  build_plot_section("📁 Reads per Sample", reads_img),

  # Orientation
  if (!is.null(ori_img))
    build_plot_section("🔀 Viral Orientation (Sense / Antisense)", ori_img)
  else "",

  # Chromosome distribution
  if (!is.null(chr_img))
    build_plot_section("🧬 Integration Site Chromosomal Distribution", chr_img)
  else "",

  # Insert length
  if (!is.null(len_img))
    build_plot_section("📏 Viral Insert Length Distribution", len_img)
  else "",

  # % identity
  if (!is.null(pid_img))
    build_plot_section("🎯 Viral Sequence % Identity to Reference", pid_img)
  else "",

  # Top genes
  if (!is.null(gene_img))
    build_plot_section("🔬 Top Integration-Bearing Genes", gene_img)
  else "",

  # Integration sites table
  sprintf('<div class="card" id="integration-table">
<h2>📋 Integration Sites (first 200 rows)</h2>
<p style="font-size:0.8rem;color:var(--fg2);margin-bottom:0.8rem;">
  Showing up to 200 rows. Full data is available in the annotated CSV files.
</p>
%s
</div>', df_to_html_table(integration_table_df))
)

# Full HTML document
html_doc <- sprintf('<!DOCTYPE html>
<html lang="en" data-theme="light">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Viral Integration Report — %s</title>
  %s
  %s
</head>
<body>
<button class="theme-toggle" id="toggle-btn" aria-label="Toggle theme">🌙 Dark mode</button>

<header>
  <h1>🧬 Viral Integration Report</h1>
  <div class="subtitle">%s</div>
  <div class="meta">Generated: %s &nbsp;|&nbsp; Pipeline: nf-viral-integration</div>
</header>

<div class="container">

  <nav class="section-nav" aria-label="Jump to section">
    <a href="#overview">Overview</a>
    <a href="#integration-table">Integration Table</a>
  </nav>

  %s

</div>

<footer>
  Generated by <strong>nf-viral-integration</strong> &mdash; Connor S. Murray, PhD &mdash; University of Louisville SOM
</footer>

</body>
</html>',
  run_name,
  html_css,
  html_js,
  run_name,
  report_date,
  sections_html
)

writeLines(html_doc, html_file)
cat(sprintf("\n✅  HTML report written to: %s\n", html_file))
cat(strrep("=", 70), "\n\n")
