#!/usr/bin/env Rscript
#
# HIV Viral Integration – Per-Sample HTML Report
# Runs at the end of INTEGRATION_ANNOTATE, after the combined CSV is finalized.
#
# Usage:
#   Rscript create_sample_report.R <combined_csv> <sample_id> <output_html>

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(scales)
})

# -----------------------------------------------------------------------------
# 1. Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript create_sample_report.R <combined_csv> <sample_id> <output_html>\n")
  quit(status = 1)
}

combined_csv <- args[1]
sample_id <- args[2]
output_html <- args[3]
report_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

if (!file.exists(combined_csv)) {
  stop("Combined CSV not found: ", combined_csv)
}

cat("\n", strrep("=", 60), "\n")
cat(sprintf("PER-SAMPLE REPORT: %s\n", sample_id))
cat(strrep("=", 60), "\n\n")
cat("Input  :", combined_csv, "\n")
cat("Output :", output_html,  "\n\n")

# -----------------------------------------------------------------------------
# 2. Load combined CSV
# -----------------------------------------------------------------------------
dat <- tryCatch(
  read_csv(combined_csv, show_col_types = FALSE, guess_max = 10000),
  error = function(e) stop("Cannot read combined CSV: ", conditionMessage(e))
)

# Coerce columns that may be mis-typed as numeric by read_csv
char_cols <- c("GENE_MATCH_STRING", "IPDA_INTACT", "IPDA_V2_INTACT",
               "COMPLETE_5PRIME", "COMPLETE_3PRIME", "EPISOME_FLAG",
               "N_GAPS_5PRIME", "N_GAPS_3PRIME", "N_GAPS_TOTAL")
for (col in intersect(char_cols, colnames(dat))) dat[[col]] <- as.character(dat[[col]])

if (!"sample" %in% colnames(dat)) dat$sample <- sample_id
cat(sprintf("Loaded %d rows\n\n", nrow(dat)))

# -----------------------------------------------------------------------------
# 3. Derive analysis columns
# -----------------------------------------------------------------------------

# 3a. Intactness factor
if ("MATCH_TYPE" %in% colnames(dat)) {
  dat <- dat %>%
    mutate(intactness = factor(
      case_when(
        MATCH_TYPE == "INTACT" ~ "Intact",
        MATCH_TYPE == "PUTATIVELY INTACT" ~ "Putatively Intact",
        MATCH_TYPE == "INDETERMINATE" ~ "Indeterminate",
        MATCH_TYPE == "INTERNAL DELETION" ~ "Internal Deletion",
        MATCH_TYPE == "TRUNCATED" ~ "Truncated",
        MATCH_TYPE == "HEAVILY TRUNCATED" ~ "Heavily Truncated",
        TRUE ~ "Other / Unknown"),
      levels = c("Intact", "Putatively Intact", "Indeterminate",
                 "Internal Deletion", "Truncated", "Heavily Truncated",
                 "Other / Unknown")))
}

# 3b. Normalise viral_orientation
if ("viral_orientation" %in% colnames(dat)) {
  dat <- dat %>%
    mutate(viral_orientation = case_when(
      viral_orientation %in% c("5prime", "+", "plus",  "sense")     ~ "+",
      viral_orientation %in% c("3prime", "-", "minus", "antisense") ~ "-",
      TRUE ~ NA_character_))
}

# 3c. Perl strand normalisation
if ("STRAND" %in% colnames(dat)) {
  dat <- dat %>%
    mutate(perl_strand = case_when(
      STRAND %in% c("plus",  "+") ~ "+",
      STRAND %in% c("minus", "-") ~ "-",
      TRUE ~ NA_character_))
} else {
  dat <- dat %>% mutate(perl_strand = NA_character_)
}

# 3d. Strand concordance
if (all(c("viral_orientation", "perl_strand") %in% colnames(dat))) {
  dat <- dat %>%
    mutate(strand_concordance = case_when(
      is.na(viral_orientation) | is.na(perl_strand) ~ "one_call_missing",
      viral_orientation == perl_strand ~ "concordant",
      TRUE ~ "discordant"))
}

# 3e. Gene-match string → per-gene presence columns
parse_gene_coverage <- function(gms) {
  genes <- c("LTR5","GAG","POL","VIF","VPR","VPU","ENV","NEF","LTR3")
  if (is.na(gms) || nchar(trimws(gms)) < 9)
    return(setNames(rep(NA_integer_, 9), genes))
  chars <- strsplit(trimws(gms), "")[[1]]
  vals  <- suppressWarnings(as.integer(chars[seq_len(9)]))
  setNames(ifelse(is.na(vals), 0L, pmin(vals, 1L)), genes)
}
if ("GENE_MATCH_STRING" %in% colnames(dat)) {
  gcm  <- do.call(rbind, lapply(dat$GENE_MATCH_STRING, parse_gene_coverage))
  dat  <- bind_cols(dat, as_tibble(gcm))
}

# 3f. IPDA flags → logical
if ("IPDA_INTACT"    %in% colnames(dat))
  dat <- mutate(dat, ipda_intact    = as.logical(suppressWarnings(as.integer(IPDA_INTACT))))
if ("IPDA_V2_INTACT" %in% colnames(dat))
  dat <- mutate(dat, ipda_v2_intact = as.logical(suppressWarnings(as.integer(IPDA_V2_INTACT))))

# 3g. Chromosome factor  (chr1..22, chrX, chrY, chrM)
all_chrs <- unique(na.omit(as.character(dat$chromosome)))
num_part <- suppressWarnings(as.integer(sub("^chr", "", all_chrs)))
num_chrs_ord <- paste0("chr", sort(unique(na.omit(num_part[!is.na(num_part)]))))
sex_mt <- intersect(c("chrX","chrY","chrM"), all_chrs)
other_chrs <- sort(setdiff(all_chrs, c(num_chrs_ord, sex_mt)))
chr_order <- c(num_chrs_ord, sex_mt, other_chrs)
if ("chromosome" %in% colnames(dat))
  dat <- mutate(dat, chromosome = factor(as.character(chromosome), levels = chr_order))

# 3h. Episomal flag
episome_col <- "EPISOME_FLAG" %in% colnames(dat)
if ("integration_site" %in% colnames(dat) || episome_col) {
  dat <- dat %>%
    mutate(is_episomal = {
      site_na <- if ("integration_site" %in% colnames(dat))
        is.na(integration_site) | integration_site == "" else rep(FALSE, n())
      perl_ep <- if (episome_col)
        !is.na(EPISOME_FLAG) & as.character(EPISOME_FLAG) %in%
        c("1","TRUE","true","yes","YES") else rep(FALSE, n())
      site_na | perl_ep
    })
}

# -----------------------------------------------------------------------------
# 4. Summary statistics
# -----------------------------------------------------------------------------
n_reads <- nrow(dat)

intactness_counts <- if ("intactness" %in% colnames(dat)) {
  count(dat, intactness, .drop = FALSE)
} else tibble()

chr_summary <- if ("chromosome" %in% colnames(dat)) {
  dat %>% filter(!is.na(chromosome)) %>% count(chromosome)
} else tibble()

gene_cols <- intersect(c("LTR5","GAG","POL","VIF","VPR","VPU","ENV","NEF","LTR3"),
                       colnames(dat))
gene_pct <- if (length(gene_cols) > 0) {
  dat %>%
    summarise(across(all_of(gene_cols),
                     ~ round(100 * mean(.x, na.rm = TRUE), 1))) %>%
    pivot_longer(everything(), names_to = "gene", values_to = "pct") %>%
    mutate(gene = factor(gene, levels = gene_cols))
} else tibble()

pos_data <- if (all(c("chromosome","integration_site") %in% colnames(dat))) {
  dat %>%
    filter(!is.na(integration_site)) %>%
    mutate(pos = as.integer(str_extract(integration_site, "(?<=:)\\d+"))) %>%
    filter(!is.na(pos), !is.na(chromosome)) %>%
    mutate(chromosome = factor(as.character(chromosome), levels = rev(chr_order)))
} else tibble()

gene_summary <- if ("gene_name" %in% colnames(dat)) {
  dat %>%
    filter(!is.na(gene_name)) %>%
    separate_rows(gene_name, sep = ";") %>%
    filter(gene_name != "") %>%
    count(gene_name, sort = TRUE) %>%
    head(20)
} else tibble()

episomal_counts <- if ("is_episomal" %in% colnames(dat)) {
  tibble(type = c("Integrated","Episomal / No Host Site"),
         n    = c(sum(!dat$is_episomal, na.rm = TRUE),
                  sum( dat$is_episomal, na.rm = TRUE)))
} else tibble()

pid_rows <- if ("percent_identity" %in% colnames(dat))
  filter(dat, !is.na(percent_identity)) else tibble()

# -----------------------------------------------------------------------------
# 5. Shared styling
# -----------------------------------------------------------------------------
INTACT_COLS <- c(
  "Intact" = "#2f9e44",
  "Putatively Intact" = "#69db7c",
  "Indeterminate" = "#fcc419",
  "Internal Deletion" = "#ff922b",
  "Truncated" = "#f03e3e",
  "Heavily Truncated" = "#862e2e",
  "Other / Unknown" = "#adb5bd")
STRAND_COLS <- c("+" = "#51cf66", 
                 "-" = "#ff6b6b")

th <- function()
  theme_minimal(base_size = 12) +
  theme(plot.background  = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        panel.grid.major = element_line(color = "#e9ecef"),
        panel.grid.minor = element_line(color = "#f1f3f5", linetype = "dashed"),
        strip.text       = element_text(face = "bold"))

use_b64 <- requireNamespace("base64enc", quietly = TRUE)

plot_to_img <- function(gg, w = 8, h = 5) {
  if (is.null(gg)) return(NULL)
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  png(tmp, width = w, height = h, units = "in", res = 150, bg = "transparent")
  print(gg)
  dev.off()
  if (use_b64) {
    sprintf('<img src="data:image/png;base64,%s" style="max-width:100%%;height:auto;">',
            base64enc::base64encode(tmp))
  } else {
    out_png <- paste0(tools::file_path_sans_ext(output_html), "_",
                      gsub("[^a-zA-Z0-9]", "", basename(tmp)), ".png")
    file.copy(tmp, out_png)
    sprintf('<img src="%s" style="max-width:100%%;height:auto;">', basename(out_png))
  }
}

# -----------------------------------------------------------------------------
# 6. Build plots
# -----------------------------------------------------------------------------

## Intactness bar
p_intact <- if (nrow(intactness_counts) > 0) {
  ggplot(intactness_counts, aes(x = intactness, y = n, fill = intactness)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    scale_fill_manual(values = INTACT_COLS, drop = FALSE) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Proviral Intactness", x = NULL, y = "Reads") +
    th() + theme(axis.text.x = element_text(angle = 35, hjust = 1))
} else NULL

## Strand orientation
p_ori <- if ("viral_orientation" %in% colnames(dat)) {
  ori_d <- dat %>% filter(!is.na(viral_orientation)) %>% count(viral_orientation)
  if (nrow(ori_d) > 0)
    ggplot(ori_d, aes(x = viral_orientation, y = n, fill = viral_orientation)) +
    geom_col(width = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = STRAND_COLS) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Viral Strand Orientation",
         x = "Strand (+ sense / − antisense)", y = "Reads") +
    th()
  else NULL
} else NULL

## Chromosome distribution
p_chr <- if (nrow(chr_summary) > 0) {
  lvls <- intersect(chr_order, as.character(chr_summary$chromosome))
  chr_summary %>%
    mutate(chromosome = factor(as.character(chromosome), levels = rev(lvls))) %>%
    ggplot(aes(x = chromosome, y = n)) +
    geom_col(fill = "#4dabf7", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Integration Sites by Chromosome",
         x = "Chromosome", y = "Reads") +
    th()
} else NULL

## Genome-wide scatter
p_genome <- if (nrow(pos_data) > 0) {
  use_color <- "intactness" %in% colnames(pos_data)
  p <- ggplot(pos_data,
              aes(x = pos / 1e6, y = chromosome,
                  color = if (use_color) intactness else I("#4dabf7"))) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_x_continuous(labels = label_number(suffix = " Mb")) +
    labs(x = "Genomic Position (Mb)", y = "Chromosome",
         color = if (use_color) "Intactness" else NULL) +
    th() + theme(axis.text.y = element_text(size = 7))
  if (use_color) p + scale_color_manual(values = INTACT_COLS) else p
} else NULL

## HIV gene coverage
p_gene_cov <- if (nrow(gene_pct) > 0) {
  ggplot(gene_pct, aes(x = gene, y = pct, fill = pct)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.0f%%", pct)), vjust = -0.3, size = 3) +
    scale_fill_gradient2(low = "#f8d7da", mid = "#ffc107", high = "#198754",
                         midpoint = 50, limits = c(0, 100)) +
    scale_y_continuous(limits = c(0, 110),
                       labels = label_percent(scale = 1)) +
    labs(title = "HIV Gene Coverage (% of reads with gene detected)",
         x = "Viral Gene", y = "% Reads") +
    th()
} else NULL

## Top host genes
p_host_genes <- if (nrow(gene_summary) > 0) {
  ggplot(gene_summary, aes(x = reorder(gene_name, n), y = n)) +
    geom_col(fill = "#f783ac", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Top Integration-Bearing Host Genes",
         x = "Gene", y = "Reads with Integration") +
    th()
} else NULL

## Insert length
p_len <- if ("INSERT_LEN" %in% colnames(dat)) {
  ld <- dat %>% filter(!is.na(INSERT_LEN), INSERT_LEN > 0)
  if (nrow(ld) > 0)
    ggplot(ld, aes(x = INSERT_LEN)) +
    geom_histogram(bins = 55, fill = "#4dabf7", alpha = 0.85) +
    scale_x_continuous(labels = label_comma()) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Viral Insert Length Distribution",
         x = "Insert Length (bp)", y = "Count") + th()
  else NULL
} else NULL

## % Identity
p_pid <- if (nrow(pid_rows) > 0) {
  ggplot(pid_rows, aes(x = percent_identity)) +
    geom_histogram(bins = 40, fill = "#51cf66", alpha = 0.85) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Viral Sequence % Identity to Reference",
         x = "% Identity", y = "Count") + th()
} else NULL

## Episomal vs integrated
p_episomal <- if (nrow(episomal_counts) > 0) {
  ggplot(episomal_counts,
         aes(x = factor(type, levels = c("Integrated","Episomal / No Host Site")),
             y = n, fill = type)) +
    geom_col(width = 0.5, show.legend = FALSE) +
    scale_fill_manual(values = c("Integrated"             = "#4dabf7",
                                 "Episomal / No Host Site" = "#ff922b")) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Episomal vs Integrated Reads", x = NULL, y = "Reads") + th()
} else NULL

# -----------------------------------------------------------------------------
# 7. Render images
# -----------------------------------------------------------------------------
img <- list(
  intact = plot_to_img(p_intact, 8, 4),
  ori = plot_to_img(p_ori, 6, 4),
  chr = plot_to_img(p_chr, 8, 6),
  genome = plot_to_img(p_genome, 12, 7),
  gene_cov = plot_to_img(p_gene_cov, 9, 4),
  host_genes = plot_to_img(p_host_genes, 8, 5),
  len = plot_to_img(p_len, 8, 4),
  pid = plot_to_img(p_pid, 8, 4),
  episomal = plot_to_img(p_episomal, 6, 4))

# -----------------------------------------------------------------------------
# 8. HTML helpers 
# -----------------------------------------------------------------------------
card <- function(id, title, body_html) {
  if (is.null(body_html) || identical(body_html, "")) return("")
  sprintf('<div class="card" id="%s"><h2>%s</h2>%s</div>', id, title, body_html)
}
plot_card <- function(id, title, img_tag) {
  if (is.null(img_tag)) return("")
  card(id, title, sprintf('<div class="plot-container">%s</div>', img_tag))
}
stat_box <- function(value, label, color = "var(--accent)")
  sprintf('<div class="stat-box"><div class="stat-value" style="color:%s">%s</div>
           <div class="stat-label">%s</div></div>', color, value, label)

df_html <- function(df, id = NULL, cap = NULL) {
  if (is.null(df) || nrow(df) == 0) return("<p><em>No data.</em></p>")
  tbl_id <- if (!is.null(id)) id else paste0("dt_", sample.int(1e9, 1))
  hdr  <- paste(sprintf("<th>%s</th>", colnames(df)), collapse = "")
  rows <- lapply(seq_len(min(nrow(df), 500L)), function(i) {
    cells <- vapply(df[i,,drop=FALSE],
                    function(x) sprintf("<td>%s</td>",
                                        ifelse(is.na(x), "", as.character(x))),
                    character(1))
    sprintf("<tr>%s</tr>", paste(cells, collapse = ""))
  })
  note <- if (nrow(df) > 500L) {
    sprintf('<p style="font-size:.78rem;color:var(--fg2);">Showing 500 of %d rows. Full data in the combined CSV.</p>', nrow(df))
  } else ""
  sprintf('%s<div class="table-wrap"><table id="%s" class="dt-table display compact"
    style="width:100%%"><caption style="display:none">%s</caption>
    <thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
          note, tbl_id, if (!is.null(cap)) cap else "", hdr, paste(rows, collapse="\n"))
}

# -----------------------------------------------------------------------------
# 9. KPI row
# -----------------------------------------------------------------------------
pct_intact <- if ("intactness" %in% colnames(dat)) {
  n_i <- sum(dat$intactness %in% c("Intact","Putatively Intact"), na.rm = TRUE)
  sprintf("%.0f%%", 100 * n_i / max(n_reads, 1L))
} else "---"

sites_n <- if ("integration_site" %in% colnames(dat)) {
  format(n_distinct(dat$integration_site, na.rm = TRUE), big.mark = ",")
} else "---"

pct_genic <- if ("gene_name" %in% colnames(dat)) {
  sprintf("%.0f%%", 100 * mean(!is.na(dat$gene_name)))
} else "---"

n_episomal <- if ("is_episomal" %in% colnames(dat)) {
  format(sum(dat$is_episomal, na.rm = TRUE), big.mark = ",")
} else "---"

n_ipda <- if ("ipda_intact" %in% colnames(dat)) {
  format(sum(dat$ipda_intact, na.rm = TRUE), big.mark = ",")
} else "---"

kpi_html <- sprintf('<div class="stats-grid">%s</div>',
                    paste(stat_box(format(n_reads, big.mark=","), "Total Reads"),
                          stat_box(sites_n, "Unique Integration Sites"),
                          stat_box(pct_intact, "Intact / Put. Intact", "#2f9e44"),
                          stat_box(pct_genic, "Reads in Host Genes"),
                          stat_box(n_episomal, "Episomal Reads", "#ff922b"),
                          stat_box(n_ipda, "IPDA Intact", "#cc5de8"),
                          collapse = "\n"))

# Integration sites display table
display_cols <- intersect(
  c("READ","chromosome","integration_site","viral_orientation","perl_strand",
    "strand_concordance","INSERT_LEN","percent_identity","intactness",
    "GENE_MATCH_STRING","IPDA_INTACT","IPDA_V2_INTACT","EPISOME_FLAG","gene_name"),
  colnames(dat))
integration_tbl <- dplyr::select(dat, all_of(display_cols))

# -----------------------------------------------------------------------------
# 10. CSS
# -----------------------------------------------------------------------------
css <- '
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<style>
:root {
  --bg:#ffffff; --bg2:#f4f6fb; --card:#ffffff; --border:#d0d7de;
  --fg:#1f2328; --fg2:#57606a; --accent:#0969da;
  --header:#24292f; --header-fg:#ffffff; --shadow:rgba(0,0,0,0.06);
  --th-bg:#f6f8fa; --tr-hover:#f6f8fa; --nav-w:210px;
}
[data-theme="dark"] {
  --bg:#0d1117; --bg2:#161b22; --card:#1c2128; --border:#30363d;
  --fg:#e6edf3; --fg2:#848d97; --accent:#58a6ff;
  --header:#161b22; --header-fg:#e6edf3; --shadow:rgba(0,0,0,0.3);
  --th-bg:#21262d; --tr-hover:#21262d;
}
*,*::before,*::after{box-sizing:border-box;}
body{margin:0;font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;
  background:var(--bg);color:var(--fg);line-height:1.6;}
header{background:var(--header);color:var(--header-fg);
  padding:2rem 2.5rem 1.5rem;border-bottom:3px solid var(--accent);}
header h1{margin:0;font-size:1.8rem;}
header .subtitle{font-size:1rem;opacity:.85;margin-top:.3rem;}
header .meta{font-size:.8rem;opacity:.65;margin-top:.5rem;}
.layout{display:flex;min-height:calc(100vh - 120px);}
.sidebar{position:sticky;top:0;height:100vh;overflow-y:auto;
  width:var(--nav-w);min-width:var(--nav-w);
  background:var(--bg2);border-right:1px solid var(--border);
  padding:1.2rem .8rem;flex-shrink:0;}
.sidebar h3{font-size:.7rem;text-transform:uppercase;letter-spacing:.08em;
  color:var(--fg2);margin:0 0 .6rem .4rem;}
.sidebar a{display:block;padding:.4em .7em;border-radius:6px;font-size:.82rem;
  color:var(--fg);text-decoration:none;margin-bottom:.15rem;transition:background .15s;}
.sidebar a:hover,.sidebar a.active{background:var(--accent);color:#fff;}
.main-content{flex:1;min-width:0;padding:1.5rem 1.8rem 3rem;}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;
  padding:1.6rem 2rem;margin-bottom:1.8rem;box-shadow:0 2px 8px var(--shadow);}
.card h2{font-size:1.1rem;font-weight:700;margin-bottom:1rem;color:var(--fg);
  border-bottom:2px solid var(--border);padding-bottom:.5rem;}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(140px,1fr));
  gap:.9rem;margin-bottom:.5rem;}
.stat-box{background:var(--bg2);border:1px solid var(--border);
  border-radius:8px;padding:1rem 1.1rem;text-align:center;}
.stat-value{font-size:1.55rem;font-weight:700;}
.stat-label{font-size:.7rem;color:var(--fg2);margin-top:.2rem;
  text-transform:uppercase;letter-spacing:.05em;}
.table-wrap{overflow-x:auto;}
table.dt-table{width:100%;border-collapse:collapse;font-size:.8rem;}
table.dt-table thead{background:var(--th-bg);}
table.dt-table th{padding:.5rem .7rem;text-align:left;font-weight:600;
  border-bottom:2px solid var(--border);color:var(--fg);white-space:nowrap;}
table.dt-table td{padding:.4rem .7rem;border-bottom:1px solid var(--border);
  color:var(--fg2);max-width:260px;overflow:hidden;
  text-overflow:ellipsis;white-space:nowrap;}
table.dt-table tr:hover td{background:var(--tr-hover);}
.dataTables_wrapper .dataTables_filter input,
.dataTables_wrapper .dataTables_length select{
  background:var(--bg2);border:1px solid var(--border);color:var(--fg);
  border-radius:4px;padding:.2em .5em;}
.dataTables_wrapper .dataTables_info,
.dataTables_wrapper .dataTables_paginate{color:var(--fg2);font-size:.78rem;margin-top:.5rem;}
.dataTables_wrapper .paginate_button{
  background:var(--bg2);border:1px solid var(--border)!important;
  color:var(--fg)!important;border-radius:4px;padding:.2em .6em;
  margin:0 .1em;cursor:pointer;}
.dataTables_wrapper .paginate_button.current,
.dataTables_wrapper .paginate_button:hover{
  background:var(--accent)!important;color:#fff!important;
  border-color:var(--accent)!important;}
.plot-container{text-align:center;}
.plot-container img{max-width:100%;border-radius:6px;}
.theme-toggle{position:fixed;top:1rem;right:1rem;z-index:200;
  background:var(--card);border:1px solid var(--border);border-radius:2em;
  padding:.35em .9em;font-size:.8rem;cursor:pointer;color:var(--fg);}
footer{text-align:center;font-size:.78rem;color:var(--fg2);
  padding:2rem 0 1rem;border-top:1px solid var(--border);}
</style>'

# -----------------------------------------------------------------------------
# 11. JavaScript 
# -----------------------------------------------------------------------------
js <- '
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script>
(function () {
  var saved = localStorage.getItem("theme") || "light";
  document.documentElement.setAttribute("data-theme", saved);
  document.addEventListener("DOMContentLoaded", function () {
    var btn = document.getElementById("toggle-btn");
    btn.textContent = saved === "dark" ? "Light mode" : "Dark mode";
    btn.addEventListener("click", function () {
      var cur = document.documentElement.getAttribute("data-theme");
      var nxt = cur === "dark" ? "light" : "dark";
      document.documentElement.setAttribute("data-theme", nxt);
      localStorage.setItem("theme", nxt);
      btn.textContent = nxt === "dark" ? "Light mode" : "Dark mode";
    });
    var links    = document.querySelectorAll(".sidebar a");
    var sections = Array.from(links).map(function(a) {
      return document.getElementById(a.getAttribute("href").slice(1));
    }).filter(Boolean);
    function onScroll() {
      var scrollY = window.scrollY + 80;
      var active  = sections[0];
      sections.forEach(function(s) { if (s.offsetTop <= scrollY) active = s; });
      links.forEach(function(a) {
        a.classList.toggle("active", a.getAttribute("href") === "#" + active.id);
      });
    }
    window.addEventListener("scroll", onScroll, { passive: true });
    onScroll();
    if ($.fn.DataTable) {
      $("table.dt-table").each(function() {
        var tbl = $(this);
        var nRows = tbl.find("tbody tr").length;
        tbl.DataTable({
          pageLength: 25,
          lengthMenu: [[10, 25, 50, 100, -1], [10, 25, 50, 100, "All"]],
          scrollX: true,
          dom: nRows > 10
            ? \'<"dt-top"lf>rt<"dt-bot"ip>\'
            : \'<"dt-top"f>rt\',
          language: {
            search:       "Search / filter:",
            lengthMenu:   "Show _MENU_ rows",
            info:         "Showing _START_–_END_ of _TOTAL_ rows",
            infoFiltered: "(filtered from _MAX_)"
          }
        });
      });
    }
  });
}());
</script>'

# -----------------------------------------------------------------------------
# 12. Nav + sections
# -----------------------------------------------------------------------------
nav_items <- c(
  "overview" = "Overview",
  "intactness" = "Intactness",
  "orientation" = "Strand",
  "episomal" = "Episomal",
  "genome-view" = "Genome View",
  "chromosomes" = "Chromosomes",
  "gene-cov" = "Gene Coverage",
  "host-genes" = "Host Genes",
  "length" = "Insert Length",
  "pct-identity"= "% Identity",
  "sites-table" = "Integration Sites")

nav_html <- paste(sprintf('<a href="#%s">%s</a>', names(nav_items), nav_items),
                  collapse = "\n")

sections <- paste0(
  card("overview",    "Sample Overview", kpi_html),
  plot_card("intactness",  "Proviral Intactness",                  img$intact),
  plot_card("orientation", "Viral Strand Orientation (+ / −)",     img$ori),
  plot_card("episomal",    "Episomal vs Integrated Reads",          img$episomal),
  plot_card("genome-view", "Integration Sites – Genome-wide View", img$genome),
  plot_card("chromosomes", "Chromosomal Distribution",              img$chr),
  plot_card("gene-cov",    "HIV Gene Coverage",                     img$gene_cov),
  plot_card("host-genes",  "Top Integration-Bearing Host Genes",   img$host_genes),
  plot_card("length",      "Viral Insert Length Distribution",      img$len),
  plot_card("pct-identity","% Identity to Reference",               img$pid),
  card("sites-table", "Integration Sites",
       paste0('<p style="font-size:.8rem;color:var(--fg2);margin-bottom:.8rem;">',
              'viral_orientation and perl_strand both show + or −. ',
              'strand_concordance = agreement between the two calling methods. ',
              'Full data in the combined CSV.</p>',
              df_html(integration_tbl, id = "sites-dt",
                      cap = "Integration Sites"))))

# -----------------------------------------------------------------------------
# 13. Write HTML
# -----------------------------------------------------------------------------
html_doc <- paste0(
  '<!DOCTYPE html>\n',
  '<html lang="en" data-theme="light">\n',
  '<head>\n',
  '  <meta charset="UTF-8">\n',
  '  <meta name="viewport" content="width=device-width,initial-scale=1.0">\n',
  sprintf('  <title>HIV Integration Report – %s</title>\n', sample_id),
  css, "\n", js, "\n",
  '</head>\n<body>\n',
  '<button class="theme-toggle" id="toggle-btn">Dark mode</button>\n',
  '<header>\n',
  sprintf('  <h1>HIV Viral Integration Report</h1>\n'),
  sprintf('  <div class="subtitle">Sample: %s</div>\n', sample_id),
  sprintf('  <div class="meta">Generated: %s | nf-viral-integration v0.2 | Connor S. Murray PhD, University of Louisville</div>\n',
          report_date),
  '</header>\n',
  '<div class="layout">\n',
  sprintf('<nav class="sidebar" aria-label="Sections"><h3>Sections</h3>%s</nav>\n', nav_html),
  '<div class="main-content">\n',
  sections,
  '</div>\n</div>\n',
  '<footer>Generated by <strong>nf-viral-integration</strong> &mdash;',
  sprintf(' Sample: %s &mdash; Connor S. Murray, PhD &mdash; University of Louisville SOM</footer>\n',
          sample_id),
  '</body>\n</html>\n')

writeLines(html_doc, output_html)
cat(sprintf("\n  Per-sample HTML written to: %s\n", output_html))
cat(strrep("=", 60), "\n\n")
