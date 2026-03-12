#!/usr/bin/env Rscript
#
# HIV Viral Integration – HTML Report Generator
#
# Usage: Rscript create_report.R <results_dir> <output_prefix> [run_name]

library(tidyverse)
library(ggplot2)
library(scales)
  
# -----------------------------------------------------------------------------
# 1. Arguments
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: Rscript create_report.R <results_dir> <output_prefix> [run_name]\n")
  quit(status = 1)
}

results_dir <- args[1]
output_prefix <- args[2]
run_name <- if (length(args) >= 3) args[3] else basename(results_dir)
report_date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

cat("\n", strrep("=", 70), "\n")
cat("HIV VIRAL INTEGRATION HTML REPORT (v4)\n")
cat(strrep("=", 70), "\n\n")

# -----------------------------------------------------------------------------
# 2. Load data
# -----------------------------------------------------------------------------
csv_files <- list.files(results_dir, pattern = "combined\\.csv$",
                        recursive = TRUE, full.names = TRUE)
if (length(csv_files) == 0) {
  cat("NOTE: No *_combined.csv found - falling back to *annotated.csv\n")
  csv_files <- list.files(results_dir, pattern = "_annotated\\.csv$",
                          recursive = TRUE, full.names = TRUE)
}
if (length(csv_files) == 0) stop("No result CSV files found in: ", results_dir)

cat(sprintf("Found %d file(s):\n", length(csv_files)))
for (f in csv_files) cat("  -", f, "\n")

COMBINED_CSV_HEADER <- paste(
  "READ,RTF_NUM,HUMAN_GROUP,INSERT,INSERT_LEN,LEFT_FLANK,RIGHT_FLANK,HUMAN_CHECK",
  "HUMAN_ALTS,HIV_DIR_ERR,FLANK_DIR_ERR,HUMAN_MAP_ERR,OVERLAP_ERR,UNMAPPED",
  "viral_sequence,viral_seq_length,viral_orientation,viral_strand,alignment_score",
  "ref_start,ref_end,percent_identity,integration_site,viral_region,chromosome,sample",
  "gene_name,gene_id,STRAND,GENE_MATCH_STRING,MATCH_TYPE,IPDA_INTACT,IPDA_V2_INTACT",
  "COMPLETE_5PRIME,N_GAPS_5PRIME,N_GAPS_3PRIME,N_GAPS_TOTAL,COMPLETE_3PRIME,EPISOME_FLAG",
  sep = ",")

repair_combined_csv_header <- function(f) {
  # Read first line to check for header
  first <- readLines(f, n = 1)
  # If the first field looks like a read ID (contains '/'), the header is missing
  if (grepl("^[^,]+/[^,]+/", first)) {
    cat(sprintf("  NOTE: %s appears to be missing its header – prepending now\n",
                basename(f)))
    lines <- readLines(f)
    tmp   <- tempfile(fileext = ".csv")
    writeLines(c(COMBINED_CSV_HEADER, lines), tmp)
    return(tmp)
  }
  f  # already has a header
}

all_data <- bind_rows(lapply(csv_files, function(f) {
  tryCatch({
    f_use <- if (grepl("combined\\.csv$", basename(f))) repair_combined_csv_header(f) else f
    d <- read_csv(f_use, show_col_types = FALSE, guess_max = 10000)
    if (!"sample" %in% colnames(d)) {
      d$sample <- sub("_?combined\\.csv$|_?annotated\\.csv$", "", basename(f))
    }
    d
  }, error = function(e) {
    cat(sprintf("  WARNING: Could not read %s: %s\n",
                basename(f), conditionMessage(e)))
    NULL
  })
}))

n_samples <- n_distinct(all_data$sample)
cat(sprintf("\nLoaded %d rows across %d sample(s)\n", nrow(all_data), n_samples))

# -----------------------------------------------------------------------------
# 3. Normalise / derive key columns
# -----------------------------------------------------------------------------

# 3a. Intactness from MATCH_TYPE (Perl output)
if ("MATCH_TYPE" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(intactness = factor(
        case_when(MATCH_TYPE == "INTACT" ~ "Intact",
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

# 3b. Normalise viral_orientation to + / - only (no "ambiguous", no "5prime")
#     R annotation script may write "5prime" / "3prime"; normalise here.
if ("viral_orientation" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(viral_orientation = case_when(
      viral_orientation %in% c("5prime", "+", "plus",  "sense") ~ "+",
      viral_orientation %in% c("3prime", "-", "minus", "antisense") ~ "-",
      TRUE ~ NA_character_))
}

# 3c. Perl strand (STRAND column from findViralGenes.pl)
if ("STRAND" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(
      perl_strand = case_when(
        STRAND %in% c("plus", "+") ~ "+",
        STRAND %in% c("minus", "-") ~ "-",
        TRUE ~ NA_character_))
} else {
  all_data <- all_data %>% mutate(perl_strand = NA_character_)
}

# 3d. Strand concordance between R (SAM/k-mer) and Perl (BLAST) calls
if ("viral_orientation" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(
      strand_concordance = case_when(
        is.na(viral_orientation) | is.na(perl_strand) ~ "one_call_missing",
        viral_orientation == perl_strand ~ "concordant",
        TRUE ~ "discordant"))
}

# 3e. Parse 9-character gene-match string into per-gene presence columns
#     Position order: LTR5 GAG POL VIF VPR VPU ENV NEF LTR3
parse_gene_coverage <- function(gms) {
  genes <- c("LTR5", "GAG", "POL", "VIF", "VPR", "VPU", "ENV", "NEF", "LTR3")
  if (is.na(gms) || nchar(trimws(gms)) < 9) {
    return(setNames(rep(NA_integer_, 9), genes))
  }
  chars <- strsplit(trimws(gms), "")[[1]]
  vals  <- suppressWarnings(as.integer(chars[seq_len(9)]))
  setNames(ifelse(is.na(vals), 0L, pmin(vals, 1L)), genes)
}

if ("GENE_MATCH_STRING" %in% colnames(all_data)) {
  gene_cov_mat <- do.call(rbind,
    lapply(all_data$GENE_MATCH_STRING, parse_gene_coverage))
  all_data <- bind_cols(all_data, as_tibble(gene_cov_mat))
}

# 3f. IPDA flags to logical
if ("IPDA_INTACT" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(ipda_intact =
             as.logical(suppressWarnings(as.integer(IPDA_INTACT))))
}
if ("IPDA_V2_INTACT" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(ipda_v2_intact =
             as.logical(suppressWarnings(as.integer(IPDA_V2_INTACT))))
}

# 3g. Chromosome factor – strictly ordered chr1..chr22, chrX, chrY, chrM
all_chrs <- unique(na.omit(as.character(all_data$chromosome)))
num_part <- suppressWarnings(as.integer(sub("^chr", "", all_chrs)))
num_chrs_ord <- paste0("chr", sort(unique(na.omit(num_part[!is.na(num_part)]))))
sex_mt <- intersect(c("chrX", "chrY", "chrM"), all_chrs)
other_chrs <- sort(setdiff(all_chrs, c(num_chrs_ord, sex_mt)))
chr_order <- c(num_chrs_ord, sex_mt, other_chrs)

if ("chromosome" %in% colnames(all_data)) {
  all_data <- all_data %>%
    mutate(chromosome = factor(as.character(chromosome),
                               levels = chr_order))
}

# 3h. Episomal flag: reads with no valid host integration site.
episome_from_col <- "EPISOME_FLAG" %in% colnames(all_data)
if ("integration_site" %in% colnames(all_data) || episome_from_col) {
  all_data <- all_data %>%
    mutate(
      is_episomal = {
        site_na <- if ("integration_site" %in% colnames(all_data)) {
                     is.na(integration_site) | integration_site == ""
                   } else { rep(FALSE, n()) }
        perl_ep  <- if (episome_from_col) {
                     !is.na(EPISOME_FLAG) & as.character(EPISOME_FLAG) %in%
                       c("1", "TRUE", "true", "yes", "YES")
                   } else { rep(FALSE, n()) }
        site_na | perl_ep
      })
}

# 3i. Reference-selection data
ref_comp_files    <- list.files(results_dir, pattern = "_mapping_comparison\\.txt$",
                                recursive = TRUE, full.names = TRUE)
ref_detail_files  <- list.files(results_dir, pattern = "_detailed_metrics\\.txt$",
                                recursive = TRUE, full.names = TRUE)

read_ref_file <- function(files, label) {
  if (length(files) == 0) return(tibble())
  bind_rows(lapply(files, function(f) {
    tryCatch({
      # Detect format: detailed_metrics.txt is a free-text report, not TSV
      first_line <- readLines(f, n = 1)
      if (grepl("^===|^Reference:", first_line)) {
        # Parse the free-text detailed metrics format
        lines   <- readLines(f)
        records <- list()
        cur     <- list()
        for (ln in lines) {
          if (grepl("^Reference:", ln)) {
            if (length(cur) > 0) records <- c(records, list(cur))
            cur <- list(reference = trimws(sub("^Reference:", "", ln)))
          } else if (grepl("^\\s+.+:", ln)) {
            kv    <- strsplit(trimws(ln), ":\\s*")[[1]]
            key   <- trimws(kv[1])
            val   <- trimws(paste(kv[-1], collapse = ":"))
            # strip trailing parenthetical e.g. "0 (0.0%)"
            val   <- trimws(sub("\\s*\\(.*\\)$", "", val))
            key   <- gsub("\\s+", "_", tolower(key))
            cur[[key]] <- suppressWarnings(as.numeric(val))
          }
        }
        if (length(cur) > 0) records <- c(records, list(cur))
        d <- bind_rows(lapply(records, as_tibble))
      } else {
        d <- read_tsv(f, show_col_types = FALSE, guess_max = 5000)
      }
      # Inject sample name from filename if not already a column
      if (!"sample" %in% tolower(colnames(d))) {
        d$sample <- sub("_mapping_comparison\\.txt$", "",
                    sub("_detailed_metrics\\.txt$",   "", basename(f)))
      }
      # Move sample to first column for readability
      d <- dplyr::select(d, sample, everything())
    }, error = function(e) {
      cat(sprintf("  WARNING: Could not read %s (%s): %s\n",
                  label, basename(f), conditionMessage(e)))
      NULL
    })
  }))
}

ref_comp_data   <- read_ref_file(ref_comp_files,   "mapping_comparison")
ref_detail_data <- read_ref_file(ref_detail_files, "detailed_metrics")

# Normalise column names to lower_snake_case for downstream use
if (nrow(ref_comp_data) > 0) {
  colnames(ref_comp_data) <- tolower(gsub("[ .]", "_", colnames(ref_comp_data)))
}
if (nrow(ref_detail_data) > 0) {
  colnames(ref_detail_data) <- tolower(gsub("[ .]", "_", colnames(ref_detail_data)))
}

cat(sprintf("Reference-selection: %d comparison rows, %d detailed rows\n",
            nrow(ref_comp_data), nrow(ref_detail_data)))

# -----------------------------------------------------------------------------
# 4. Summary statistics
# -----------------------------------------------------------------------------
sample_counts <- all_data %>%
  count(sample, name = "n_reads") %>%
  arrange(desc(n_reads))

# Detect samples with too few reads to have been fully analyzed
LOW_READ_THRESHOLD <- 2L
low_read_samples <- sample_counts %>% filter(n_reads < LOW_READ_THRESHOLD)
if (nrow(low_read_samples) > 0) {
  cat(sprintf(
    "\nWARNING: %d sample(s) have fewer than %d read(s) and were likely not analyzed:\n",
    nrow(low_read_samples), LOW_READ_THRESHOLD))
  for (i in seq_len(nrow(low_read_samples))) {
    cat(sprintf("  - %s (%d read(s))\n",
                low_read_samples$sample[i], low_read_samples$n_reads[i]))
  }
}

if ("chromosome" %in% colnames(all_data)) {
  chr_summary <- all_data %>%
    filter(!is.na(chromosome)) %>%
    count(chromosome, sort = FALSE) %>%   # preserve factor order, not count order
    mutate(chromosome = factor(chromosome, levels = chr_order))
} else {
  chr_summary <- tibble()
}

if ("INSERT_LEN" %in% colnames(all_data)) {
  len_stats <- all_data %>%
    filter(!is.na(INSERT_LEN), INSERT_LEN > 0) %>%
    group_by(sample) %>%
    summarise(n = n(),
              mean_len = round(mean(INSERT_LEN)),
              median_len = round(median(INSERT_LEN)),
              min_len = min(INSERT_LEN),
              max_len = max(INSERT_LEN),
              .groups = "drop")} else {
  len_stats <- tibble()
}

if ("viral_orientation" %in% colnames(all_data)) {
  ori_summary <- all_data %>%
    filter(!is.na(viral_orientation)) %>%
    count(sample, viral_orientation) %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()} else {
  ori_summary <- tibble()
}

if ("strand_concordance" %in% colnames(all_data)) {
  concordance_summary <- all_data %>%
    filter(!is.na(strand_concordance)) %>%
    count(sample, strand_concordance) %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()} else {
  concordance_summary <- tibble()
}

if ("intactness" %in% colnames(all_data)) {
  intactness_summary <- all_data %>%
    count(sample, intactness) %>%
    group_by(sample) %>%
    mutate(pct = round(100 * n / sum(n), 1)) %>%
    ungroup()} else {
  intactness_summary <- tibble()
}

if ("gene_name" %in% colnames(all_data)) {
  gene_summary <- all_data %>%
    filter(!is.na(gene_name)) %>%
    separate_rows(gene_name, sep = ";") %>%
    filter(gene_name != "") %>%
    count(gene_name, sort = TRUE) %>%
    head(20)} else {
  gene_summary <- tibble()
}

if ("ipda_intact" %in% colnames(all_data)) {
  ipda_summary <- all_data %>%
    group_by(sample) %>%
    summarise(
      n_ipda_intact = sum(ipda_intact, na.rm = TRUE),
      n_ipda_v2_intact = if ("ipda_v2_intact" %in% colnames(all_data))
                           sum(ipda_v2_intact, na.rm = TRUE) else 0L,
      .groups = "drop")} else {
  ipda_summary <- tibble()
}

# Diversity stats – computed only from rows with real percent_identity values.
if ("percent_identity" %in% colnames(all_data)) {
  pid_rows <- all_data %>% filter(!is.na(percent_identity))
  if (nrow(pid_rows) > 0) {
    diversity_stats <- pid_rows %>%
      group_by(sample) %>%
      summarise(n_seqs = n(),
                mean_pid = round(mean(percent_identity), 2),
                median_pid = round(median(percent_identity), 2),
                sd_pid = round(sd(percent_identity), 2),
                iqr_pid = round(IQR(percent_identity), 2),
                .groups = "drop")
  } else {
    diversity_stats <- tibble()
  }
} else {
  pid_rows <- tibble()
  diversity_stats <- tibble()
}

gene_cols <- intersect(
  c("LTR5","GAG","POL","VIF","VPR","VPU","ENV","NEF","LTR3"),
  colnames(all_data))

if (length(gene_cols) > 0) {
  genome_cov_summary <- all_data %>%
    group_by(sample) %>%
    summarise(across(all_of(gene_cols),
                     ~ round(100 * mean(.x, na.rm = TRUE), 1)),
              .groups = "drop")} else {
  genome_cov_summary <- tibble()
}

# Episomal per-sample counts
if ("is_episomal" %in% colnames(all_data)) {
  episomal_summary <- all_data %>%
    group_by(sample) %>%
    summarise(n_episomal = sum( is_episomal, na.rm = TRUE),
              n_integrated = sum(!is_episomal, na.rm = TRUE),
              pct_episomal = round(100 * mean(is_episomal, na.rm = TRUE), 1),
              .groups = "drop")
} else {
  episomal_summary <- tibble()
}

# -----------------------------------------------------------------------------
# 5. Plot helpers
# -----------------------------------------------------------------------------
use_base64 <- requireNamespace("base64enc", quietly = TRUE)
img_dir <- paste0(output_prefix, "_report_files")
if (!use_base64) dir.create(img_dir, showWarnings = FALSE)

plot_to_base64 <- function(gg, w = 8, h = 5, res = 150) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  png(tmp, width = w, height = h, units = "in", res = res, bg = "transparent")
  print(gg)
  dev.off()
  sprintf(
    '<img src="data:image/png;base64,%s" style="max-width:100%%;height:auto;">',
    base64enc::base64encode(tmp))
}

plot_to_file <- function(gg, label, w = 8, h = 5) {
  fn <- file.path(img_dir, paste0(label, ".png"))
  png(fn, width = w, height = h, units = "in", res = 150, bg = "transparent")
  print(gg)
  dev.off()
  sprintf('<img src="%s" style="max-width:100%%;height:auto;">', basename(fn))
}

make_img <- function(gg, label, w = 8, h = 5) {
  if (is.null(gg)) return(NULL)
  if (use_base64) plot_to_base64(gg, w, h) else plot_to_file(gg, label, w, h)
}

# Colour palettes
PALETTE <- c("#4dabf7","#51cf66","#ff6b6b","#fcc419","#cc5de8",
             "#20c997","#f783ac","#a9e34b","#74c0fc","#ff922b")

INTACT_COLS <- c(
  "Intact" = "#2f9e44",
  "Putatively Intact" = "#69db7c",
  "Indeterminate" = "#fcc419",
  "Internal Deletion" = "#ff922b",
  "Truncated" = "#f03e3e",
  "Heavily Truncated" = "#862e2e",
  "Other / Unknown" = "#adb5bd")

STRAND_COLS  <- c("+" = "#51cf66", "-" = "#ff6b6b")

CONCORD_COLS <- c("concordant" = "#2f9e44",
                  "discordant" = "#e03131",
                  "one_call_missing" = "#adb5bd")

th <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      panel.grid.major = element_line(color = "#e9ecef"),
      panel.grid.minor = element_line(color = "#f1f3f5", linetype = "dashed"),
      strip.text = element_text(face = "bold"))}

# -----------------------------------------------------------------------------
# 6. Build all plots
# -----------------------------------------------------------------------------

## 6a. Reads per sample
p_reads <- {
  ggplot(sample_counts,
         aes(x = reorder(sample, n_reads), y = n_reads, fill = sample)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    coord_flip() +
    scale_fill_manual(
      values = colorRampPalette(PALETTE)(nrow(sample_counts))) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Reads per Sample", x = NULL, y = "Count") +
    th()
}

## 6b. Intactness stacked proportion
if (nrow(intactness_summary) > 0) {
  p_intact <- ggplot(intactness_summary,
                     aes(x = sample, y = n, fill = intactness)) +
    geom_col(position = "fill", width = 0.7) +
    scale_y_continuous(labels = label_percent()) +
    scale_fill_manual(values = INTACT_COLS, na.value = "#adb5bd",
                      drop = FALSE) +
    labs(title = "Proviral Intactness per Sample",
         x = NULL, y = "Proportion", fill = "Category") +
    th() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p_intact <- NULL
}

## 6c. Intactness absolute counts
if (nrow(intactness_summary) > 0) {
  p_intact_abs <- ggplot(intactness_summary,
                         aes(x = intactness, y = n, fill = intactness)) +
    geom_col(show.legend = FALSE, width = 0.7) +
    facet_wrap(~sample, scales = "free_y", ncol = min(n_samples, 3)) +
    scale_fill_manual(values = INTACT_COLS, drop = FALSE) +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Proviral Intactness - Absolute Counts",
         x = NULL, y = "Reads") +
    th() +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8))
} else {
  p_intact_abs <- NULL
}

## 6d. Viral strand orientation (+/-)
if (nrow(ori_summary) > 0) {
  p_ori <- ggplot(ori_summary,
                  aes(x = sample, y = n, fill = viral_orientation)) +
    geom_col(position = "fill", width = 0.7) +
    scale_y_continuous(labels = label_percent()) +
    scale_fill_manual(values = STRAND_COLS, na.value = "#adb5bd",
                      labels = c("+" = "+ (sense / forward)",
                                 "-" = "- (antisense / RC)")) +
    labs(title = "Viral Strand Orientation (+ / -)",
         x = NULL, y = "Proportion", fill = "Strand") +
    th() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p_ori <- NULL
}

## 6f. Chromosome distribution – ordered chr1..chr22, chrX, chrY, chrM
if (nrow(chr_summary) > 0) {
  chr_levels_present <- intersect(chr_order, as.character(chr_summary$chromosome))
  p_chr <- chr_summary %>%
    mutate(chromosome = factor(as.character(chromosome),
                               levels = rev(chr_levels_present))) %>%
    ggplot(aes(x = chromosome, y = n)) +
    geom_col(fill = "#4dabf7", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Integration Sites - Chromosomal Distribution",
         x = "Chromosome", y = "Number of Reads") +
    th()
} else {
  p_chr <- NULL
}

## 6g. Genome-wide Manhattan-style plot – y-axis chr1 (top) to chrM (bottom)
if (all(c("chromosome", "integration_site") %in% colnames(all_data))) {
  pos_data <- all_data %>%
    filter(!is.na(integration_site)) %>%
    mutate(pos = as.integer(str_extract(integration_site, "(?<=:)\\d+"))) %>%
    filter(!is.na(pos), !is.na(chromosome))

  if (nrow(pos_data) > 0) {
    pos_data <- pos_data %>%
      mutate(chromosome = factor(as.character(chromosome),
                                 levels = rev(chr_order)))

    if ("intactness" %in% colnames(pos_data)) {
      color_aes   <- aes(x = pos / 1e6, y = chromosome, color = intactness)
      color_scale <- scale_color_manual(values = INTACT_COLS, name = "Intactness")
    } else {
      color_aes   <- aes(x = pos / 1e6, y = chromosome, color = sample)
      color_scale <- scale_color_manual(
        values = colorRampPalette(PALETTE)(n_samples), name = "Sample")
    }
    p_genome <- ggplot(pos_data, color_aes) +
      geom_point(alpha = 0.7, size = 1.8) +
      scale_x_continuous(labels = label_number(suffix = " Mb")) +
      color_scale +
      labs(title = "Integration Sites - Genome-wide View",
           x = "Position (Mb)", y = "Chromosome") +
      th() +
      theme(axis.text.y = element_text(size = 7))
  } else {
    p_genome <- NULL
  }
} else {
  p_genome <- NULL
}

## 6h. Insert length distribution
if ("INSERT_LEN" %in% colnames(all_data)) {
  ld <- all_data %>% filter(!is.na(INSERT_LEN), INSERT_LEN > 0)
  if (nrow(ld) > 0) {
    p_len <- ggplot(ld, aes(x = INSERT_LEN, fill = sample)) +
      geom_histogram(bins = 60, alpha = 0.85, position = "identity") +
      scale_x_continuous(labels = label_comma()) +
      scale_y_continuous(labels = label_comma()) +
      scale_fill_manual(
        values = colorRampPalette(PALETTE)(n_samples)) +
      facet_wrap(~sample, scales = "free_y", ncol = min(n_samples, 3)) +
      labs(title = "Viral Insert Length Distribution",
           x = "Insert Length (bp)", y = "Count") +
      th()
  } else {
    p_len <- NULL
  }
} else {
  p_len <- NULL
}

## 6i. Percent identity histogram
if ("percent_identity" %in% colnames(all_data)) {
  pd <- all_data %>% filter(!is.na(percent_identity))
  if (nrow(pd) > 0) {
    p_pid <- ggplot(pd, aes(x = percent_identity, fill = sample)) +
      geom_histogram(bins = 40, alpha = 0.85, position = "identity") +
      scale_fill_manual(
        values = colorRampPalette(PALETTE)(n_samples)) +
      facet_wrap(~sample, scales = "free_y", ncol = min(n_samples, 3)) +
      labs(title = "Viral Sequence % Identity to Reference",
           x = "% Identity", y = "Count") +
      th()
  } else {
    p_pid <- NULL
  }
} else {
  p_pid <- NULL
}

## 6j. HIV gene coverage heatmap
if (nrow(genome_cov_summary) > 0 && length(gene_cols) >= 3) {
  gc_long <- genome_cov_summary %>%
    pivot_longer(all_of(gene_cols),
                 names_to = "gene", values_to = "pct") %>%
    mutate(gene = factor(gene, levels = gene_cols))

  p_gene_cov <- ggplot(gc_long, aes(x = gene, y = sample, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.0f%%", pct)),
              size = 3, color = "black") +
    scale_fill_gradient2(low = "#f8d7da", mid = "#ffc107",
                         high = "#198754", midpoint = 50,
                         limits = c(0, 100),
                         name = "% Reads\nwith gene") +
    labs(title = "HIV Gene Coverage Across Samples",
         x = "Viral Gene", y = "Sample") +
    th()
} else {
  p_gene_cov <- NULL
}

## 6k. Top integration-bearing host genes
if (nrow(gene_summary) > 0) {
  p_host_genes <- ggplot(gene_summary,
                         aes(x = reorder(gene_name, n), y = n)) +
    geom_col(fill = "#f783ac", width = 0.7) +
    coord_flip() +
    scale_y_continuous(labels = label_comma()) +
    labs(title = "Top Integration-Bearing Host Genes",
         x = "Gene", y = "Reads with Integration") +
    th()
} else {
  p_host_genes <- NULL
}

## 6l. Diversity violin + boxplot
if (nrow(pid_rows) > 0) {
  p_diversity <- ggplot(pid_rows,
                        aes(x = sample, y = percent_identity, fill = sample)) +
    geom_violin(alpha = 0.5, show.legend = FALSE) +
    geom_boxplot(width = 0.18, alpha = 0.9, outlier.size = 0.6,
                 show.legend = FALSE) +
    scale_fill_manual(
      values = colorRampPalette(PALETTE)(n_samples)) +
    scale_y_continuous(limits = c(NA, 100)) +
    labs(title = "Viral Sequence Diversity (% Identity to Reference)",
         subtitle = "Lower = more divergent from reference",
         x = NULL, y = "% Identity") +
    th() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p_diversity <- NULL
}

## 6m. Episomal vs integrated stacked bar
if (nrow(episomal_summary) > 0) {
  ep_long <- episomal_summary %>%
    pivot_longer(cols = c(n_integrated, n_episomal),
                 names_to = "type", values_to = "count") %>%
    mutate(type = factor(recode(type,
                                n_integrated = "Integrated",
                                n_episomal   = "Episomal / No Host Site"),
                         levels = c("Integrated", "Episomal / No Host Site")))

  p_episomal <- ggplot(ep_long, aes(x = sample, y = count, fill = type)) +
    geom_col(position = "fill", width = 0.7) +
    scale_y_continuous(labels = label_percent()) +
    scale_fill_manual(values = c("Integrated" = "#4dabf7",
                                 "Episomal / No Host Site" = "#ff922b")) +
    labs(title = "Episomal vs Integrated Reads per Sample",
         x = NULL, y = "Proportion", fill = NULL) +
    th() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
} else {
  p_episomal <- NULL
}

## 6n. Reference-selection: mapped-reads bar per sample x reference
if (nrow(ref_comp_data) > 0) {
  cnt_col <- intersect(c("mapped_reads","total_mapped","n_mapped",
                          "reads_mapped","mapped","count"),
                        colnames(ref_comp_data))[1]
  ref_col <- intersect(c("reference","viral_reference","ref","ref_name",
                          "reference_name","genome"),
                        colnames(ref_comp_data))[1]
  smp_col <- intersect(c("sample","sample_id"), colnames(ref_comp_data))[1]

  if (!is.na(cnt_col) && !is.na(ref_col) && !is.na(smp_col)) {
    p_ref_reads <- ref_comp_data %>%
      rename(ref_reads = !!cnt_col,
             reference = !!ref_col,
             smp       = !!smp_col) %>%
      mutate(ref_reads = suppressWarnings(as.numeric(ref_reads)),
             reference = str_trunc(as.character(reference), 40)) %>%
      filter(!is.na(ref_reads)) %>%
      ggplot(aes(x = reorder(reference, ref_reads), y = ref_reads,
                 fill = smp)) +
      geom_col(position = "dodge", width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = colorRampPalette(PALETTE)(n_samples),
                        name = "Sample") +
      scale_y_continuous(labels = label_comma()) +
      labs(title = "Reference-Selection: Mapped Reads per Viral Reference",
           subtitle = "Taller bar = better-matching reference for that sample",
           x = "Viral Reference", y = "Mapped Reads") +
      th()
  } else {
    p_ref_reads <- NULL
  }
} else {
  p_ref_reads <- NULL
}

## 6o. Reference-selection: heatmap of a quality metric across all candidates
if (nrow(ref_comp_data) > 0) {
  heat_col <- intersect(c("coverage","breadth_coverage","pct_coverage",
                           "percent_coverage","mean_depth","depth",
                           "avg_identity","mean_identity","identity",
                           "score","avg_mapq","avg_edit_dist"),
                         colnames(ref_comp_data))[1]
  ref_col2 <- intersect(c("reference","viral_reference","ref","ref_name",
                           "reference_name","genome"),
                         colnames(ref_comp_data))[1]
  smp_col2 <- intersect(c("sample","sample_id"), colnames(ref_comp_data))[1]

  if (!is.na(heat_col) && !is.na(ref_col2) && !is.na(smp_col2)) {
    heat_label <- tools::toTitleCase(gsub("_", " ", heat_col))
    p_ref_heat <- ref_comp_data %>%
      rename(metric    = !!heat_col,
             reference = !!ref_col2,
             smp       = !!smp_col2) %>%
      mutate(metric    = suppressWarnings(as.numeric(metric)),
             reference = str_trunc(as.character(reference), 40)) %>%
      filter(!is.na(metric)) %>%
      ggplot(aes(x = smp, y = reference, fill = metric)) +
      geom_tile(color = "white", linewidth = 0.3) +
      scale_fill_gradient2(low = "#f8d7da", mid = "#ffc107",
                           high = "#198754", midpoint = median(
                             suppressWarnings(as.numeric(
                               ref_comp_data[[heat_col]]))[!is.na(
                               suppressWarnings(as.numeric(
                                 ref_comp_data[[heat_col]])))]),
                           name = heat_label) +
      labs(title = sprintf("Reference-Selection: %s per Candidate", heat_label),
           x = "Sample", y = "Viral Reference") +
      th() +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            axis.text.y = element_text(size = 7))
  } else {
    p_ref_heat <- NULL
  }
} else {
  p_ref_heat <- NULL
}

# -----------------------------------------------------------------------------
# 7. Render images
# -----------------------------------------------------------------------------
img <- list(
  reads = make_img(p_reads, "reads", 8, 4),
  intact = make_img(p_intact, "intactness", 9, 5),
  intact_abs = make_img(p_intact_abs, "intact_abs", 10,
                        max(4, n_samples * 1.5)),
  ori = make_img(p_ori, "orientation", 8, 4),
  chr = make_img(p_chr, "chromosome", 8, 6),
  genome = make_img(p_genome, "genome", 12, 7),
  len = make_img(p_len, "insert_len", 10, 4),
  pid = make_img(p_pid, "pct_id", 10, 4),
  gene_cov = make_img(p_gene_cov, "gene_cov", 9,
                        max(4, n_samples * 0.8 + 2)),
  host_genes = make_img(p_host_genes, "host_genes", 8, 5),
  diversity = make_img(p_diversity, "diversity", 8, 5),
  episomal = make_img(p_episomal, "episomal", 8, 4),
  ref_reads = make_img(p_ref_reads, "ref_reads", 10,
                        max(5, nrow(ref_comp_data) * 0.35 + 2)),
  ref_heat = make_img(p_ref_heat, "ref_heat", 10,
                        max(5, nrow(ref_comp_data) * 0.35 + 2)))

# -----------------------------------------------------------------------------
# 8. HTML helpers
# -----------------------------------------------------------------------------
card <- function(id, title, body_html) {
  if (is.null(body_html) || identical(body_html, "")) return("")
  sprintf('<div class="card" id="%s"><h2>%s</h2>%s</div>',
          id, title, body_html)
}

plot_card <- function(id, title, img_tag) {
  if (is.null(img_tag)) return("")
  card(id, title,
       sprintf('<div class="plot-container">%s</div>', img_tag))
}

df_html <- function(df, id = NULL, cap = NULL) {
  if (is.null(df) || nrow(df) == 0) return("<p><em>No data.</em></p>")
  tbl_id <- if (!is.null(id)) id else paste0("dt_", sample.int(1e9, 1))
  cap_h  <- if (!is.null(cap)) sprintf('<caption>%s</caption>', cap) else ""
  hdr    <- paste(sprintf("<th>%s</th>", colnames(df)), collapse = "")
  rows   <- lapply(seq_len(nrow(df)), function(i) {
    cells <- vapply(df[i, , drop = FALSE], function(x)
      sprintf("<td>%s</td>", ifelse(is.na(x), "", as.character(x))),
      character(1))
    sprintf("<tr>%s</tr>", paste(cells, collapse = ""))
  })
  sprintf(
    '<div class="table-wrap"><table id="%s" class="dt-table display compact" style="width:100%%"><caption style="display:none">%s</caption><thead><tr>%s</tr></thead><tbody>%s</tbody></table></div>',
    tbl_id, if (!is.null(cap)) cap else "", hdr, paste(rows, collapse = "\n"))
}

stat_box <- function(value, label, color = "var(--accent)") {
  sprintf(
    '<div class="stat-box"><div class="stat-value" style="color:%s">%s</div><div class="stat-label">%s</div></div>',
    color, value, label)
}

# -----------------------------------------------------------------------------
# 9. KPI row
# -----------------------------------------------------------------------------
total_reads <- sum(sample_counts$n_reads)

if ("intactness" %in% colnames(all_data)) {
  n_intact   <- sum(all_data$intactness %in%
                      c("Intact","Putatively Intact"), na.rm = TRUE)
  pct_intact <- sprintf("%.0f%%", 100 * n_intact / total_reads)
} else {
  pct_intact <- "---"
}

if ("integration_site" %in% colnames(all_data)) {
  sites_label <- format(n_distinct(all_data$integration_site, na.rm = TRUE),
                        big.mark = ",")
} else {
  sites_label <- "---"
}

if ("gene_name" %in% colnames(all_data)) {
  pct_genic <- sprintf("%.0f%%", 100 * mean(!is.na(all_data$gene_name)))
} else {
  pct_genic <- "---"
}

if ("is_episomal" %in% colnames(all_data)) {
  episomal_label <- format(sum(all_data$is_episomal, na.rm = TRUE),
                           big.mark = ",")
} else {
  episomal_label <- "---"
}

n_chromosomes <- n_distinct(all_data$chromosome, na.rm = TRUE)

kpi_html <- sprintf('<div class="stats-grid">%s</div>',
  paste(
    stat_box(format(total_reads, big.mark = ","), "Total Reads"),
    stat_box(n_samples, "Samples"),
    stat_box(sites_label, "Unique Integration Sites"),
    stat_box(pct_intact, "Intact / Put. Intact", "#2f9e44"),
    stat_box(pct_genic, "Reads in Host Genes"),
    stat_box(episomal_label, "Episomal Reads", "#ff922b"),
    stat_box(n_chromosomes, "Chromosomes Affected"),
    collapse = "\n"))

# QC table
qc_df <- sample_counts
if (nrow(len_stats) > 0) {
  qc_df <- left_join(qc_df,
                     len_stats %>% dplyr::select(sample, mean_len, median_len),
                     by = "sample")
}
if (nrow(ipda_summary) > 0) {
  qc_df <- left_join(qc_df, ipda_summary, by = "sample")
}
if (nrow(diversity_stats) > 0) {
  qc_df <- left_join(qc_df,
                     diversity_stats %>% dplyr::select(sample, mean_pid, sd_pid),
                     by = "sample")
}
if (nrow(episomal_summary) > 0) {
  qc_df <- left_join(qc_df,
                     episomal_summary %>% dplyr::select(sample, n_episomal, pct_episomal),
                     by = "sample")
}

# Integration sites table
display_cols <- intersect(
  c("sample", "READ", "chromosome", "integration_site",
    "viral_orientation", "perl_strand", "strand_concordance",
    "INSERT_LEN", "percent_identity",
    "intactness", "GENE_MATCH_STRING",
    "IPDA_INTACT", "IPDA_V2_INTACT", "EPISOME_FLAG",
    "gene_name"),
  colnames(all_data))
integration_tbl <- all_data %>%
  dplyr::select(all_of(display_cols))

# -----------------------------------------------------------------------------
# 10. CSS – left-sidebar navigation layout
# -----------------------------------------------------------------------------
css <- '
<link rel="stylesheet" href="https://cdn.datatables.net/1.13.7/css/jquery.dataTables.min.css">
<style>
:root {
  --bg:#ffffff; --bg2:#f4f6fb; --card:#ffffff; --border:#d0d7de;
  --fg:#1f2328; --fg2:#57606a; --accent:#0969da;
  --header:#24292f; --header-fg:#ffffff; --shadow:rgba(0,0,0,0.06);
  --th-bg:#f6f8fa; --tr-hover:#f6f8fa;
  --nav-w:210px;
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
.sidebar{
  position:sticky;top:0;height:100vh;overflow-y:auto;
  width:var(--nav-w);min-width:var(--nav-w);
  background:var(--bg2);border-right:1px solid var(--border);
  padding:1.2rem .8rem;flex-shrink:0;}
.sidebar h3{font-size:.7rem;text-transform:uppercase;letter-spacing:.08em;
  color:var(--fg2);margin:0 0 .6rem .4rem;}
.sidebar a{
  display:block;padding:.4em .7em;border-radius:6px;
  font-size:.82rem;color:var(--fg);text-decoration:none;
  margin-bottom:.15rem;transition:background .15s;}
.sidebar a:hover,.sidebar a.active{background:var(--accent);color:#fff;}
.main-content{flex:1;min-width:0;padding:1.5rem 1.8rem 3rem;}
.card{background:var(--card);border:1px solid var(--border);border-radius:10px;
  padding:1.6rem 2rem;margin-bottom:1.8rem;box-shadow:0 2px 8px var(--shadow);}
.card h2{font-size:1.1rem;font-weight:700;margin-bottom:1rem;color:var(--fg);
  border-bottom:2px solid var(--border);padding-bottom:.5rem;}
.stats-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));
  gap:1rem;margin-bottom:.5rem;}
.stat-box{background:var(--bg2);border:1px solid var(--border);
  border-radius:8px;padding:1rem 1.2rem;text-align:center;}
.stat-value{font-size:1.6rem;font-weight:700;}
.stat-label{font-size:.72rem;color:var(--fg2);margin-top:.2rem;
  text-transform:uppercase;letter-spacing:.05em;}
.table-wrap{overflow-x:auto;}
table.dt-table{width:100%;border-collapse:collapse;font-size:.8rem;}
table.dt-table thead{background:var(--th-bg);}
table.dt-table th{padding:.5rem .7rem;text-align:left;font-weight:600;
  border-bottom:2px solid var(--border);color:var(--fg);white-space:nowrap;}
table.dt-table td{padding:.4rem .7rem;border-bottom:1px solid var(--border);color:var(--fg2);
  max-width:260px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;}
table.dt-table tr:hover td{background:var(--tr-hover);}
/* DataTables widget overrides */
.dataTables_wrapper .dataTables_filter input,
.dataTables_wrapper .dataTables_length select{
  background:var(--bg2);border:1px solid var(--border);color:var(--fg);
  border-radius:4px;padding:.2em .5em;}
.dataTables_wrapper .dataTables_info,
.dataTables_wrapper .dataTables_paginate{color:var(--fg2);font-size:.78rem;margin-top:.5rem;}
.dataTables_wrapper .paginate_button{
  background:var(--bg2);border:1px solid var(--border)!important;
  color:var(--fg)!important;border-radius:4px;padding:.2em .6em;margin:0 .1em;cursor:pointer;}
.dataTables_wrapper .paginate_button.current,
.dataTables_wrapper .paginate_button:hover{
  background:var(--accent)!important;color:#fff!important;border-color:var(--accent)!important;}
/* Low-read warning banner */
.low-read-warning{
  background:#fff3cd;border:1px solid #ffc107;border-radius:8px;
  padding:.8rem 1.2rem;margin-bottom:1rem;font-size:.85rem;color:#856404;}
[data-theme="dark"] .low-read-warning{
  background:#2d2200;border-color:#856404;color:#ffc107;}
.low-read-warning strong{display:block;margin-bottom:.3rem;}
.low-read-warning ul{margin:.2rem 0 0 1.2rem;padding:0;}
.plot-container{text-align:center;}
.plot-container img{max-width:100%;border-radius:6px;}
.theme-toggle{position:fixed;top:1rem;right:1rem;z-index:200;
  background:var(--card);border:1px solid var(--border);border-radius:2em;
  padding:.35em .9em;font-size:.8rem;cursor:pointer;color:var(--fg);}
footer{text-align:center;font-size:.78rem;color:var(--fg2);
  padding:2rem 0 1rem;border-top:1px solid var(--border);}
</style>
'

# -----------------------------------------------------------------------------
# 11. JavaScript (theme toggle + sidebar active-link tracking)
# -----------------------------------------------------------------------------
js <- '
<script src="https://code.jquery.com/jquery-3.7.0.min.js"></script>
<script src="https://cdn.datatables.net/1.13.7/js/jquery.dataTables.min.js"></script>
<script>
(function () {
  /* ── Theme toggle ─────────────────────────────────────────────── */
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

    /* ── Sidebar active-link tracking ────────────────────────────── */
    var links    = document.querySelectorAll(".sidebar a");
    var sections = Array.from(links).map(function(a) {
      return document.getElementById(a.getAttribute("href").slice(1));
    }).filter(Boolean);
    function onScroll() {
      var scrollY = window.scrollY + 80;
      var active  = sections[0];
      sections.forEach(function(s) { if (s.offsetTop <= scrollY) active = s; });
      links.forEach(function(a) {
        a.classList.toggle("active",
          a.getAttribute("href") === "#" + active.id);
      });
    }
    window.addEventListener("scroll", onScroll, { passive: true });
    onScroll();

    /* ── DataTables init ──────────────────────────────────────────── */
    /* Integration Sites table: full-featured with many rows */
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
</script>
'

# -----------------------------------------------------------------------------
# 12. Nav items + sections
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
  "diversity" = "Diversity",
  "ref-select" = "Ref Selection",
  "sites-table" = "Sites Table")

nav_html <- paste(
  sprintf('<a href="#%s">%s</a>', names(nav_items), nav_items),
  collapse = "\n")

# Low-read warning banner (injected into Overview card)
low_read_html <- if (nrow(low_read_samples) > 0) {
  items <- paste(
    sprintf("<li><strong>%s</strong> — %d read(s) (below threshold of %d)</li>",
            low_read_samples$sample, low_read_samples$n_reads,
            LOW_READ_THRESHOLD),
    collapse = "\n")
  sprintf(
    '<div class="low-read-warning"><strong>⚠ Low-Read Samples — Analysis May Be Incomplete</strong>
The following sample(s) had fewer than %d read(s) in the combined CSV and were likely
skipped or partially processed during the pipeline run:
<ul>%s</ul>
Check upstream INTEGRATION_ANNOTATE / SELECT_BEST_REFERENCE outputs for these samples.
</div>', LOW_READ_THRESHOLD, items)
} else {
  ""
}

sec_overview  <- card("overview", "Run Overview",
  paste0(low_read_html,
         kpi_html,
         df_html(qc_df, id = "qc-table",
                 cap = "Per-sample QC Summary")))
sec_reads     <- plot_card("reads-per-sample",
                           "Reads per Sample", img$reads)
sec_intact    <- plot_card("intactness",
                           "Proviral Intactness (Proportion)", img$intact)
sec_intact_a  <- plot_card("intactness-abs",
                           "Proviral Intactness (Absolute Counts)",
                           img$intact_abs)
sec_ori       <- plot_card("orientation",
                           "Viral Strand Orientation (+ / -)", img$ori)
sec_episomal  <- plot_card("episomal",
                           "Episomal vs Integrated Reads", img$episomal)
sec_genome    <- plot_card("genome-view",
                           "Integration Sites - Genome-wide", img$genome)
sec_chr       <- plot_card("chromosomes",
                           "Chromosomal Distribution", img$chr)
sec_genecov   <- plot_card("gene-cov",
                           "HIV Gene Coverage per Sample", img$gene_cov)
sec_hostgene  <- plot_card("host-genes",
                           "Top Integration-Bearing Host Genes",
                           img$host_genes)
sec_len       <- plot_card("length",
                           "Viral Insert Length Distribution", img$len)
sec_div       <- plot_card("diversity",
                           "Sequence Diversity (violin + box)", img$diversity)
sec_pid       <- plot_card("pct-identity",
                           "% Identity to Reference", img$pid)

# Reference-selection section: table + plots
ref_sel_body <- {
  parts <- character(0)
  if (nrow(ref_comp_data) > 0) {
    parts <- c(parts,
      '<p style="font-size:.8rem;color:var(--fg2);margin-bottom:.8rem;">',
      'Mapping statistics for every viral reference candidate evaluated by ',
      'SELECT_BEST_REFERENCE.  The selected reference for each sample is the ',
      'one with the most mapped reads after duplicate removal.</p>',
      df_html(ref_comp_data, id = "ref-comp-table",
              cap = "Mapping Comparison – All Reference Candidates"))
  }
  if (!is.null(img$ref_reads)) {
    parts <- c(parts,
      sprintf('<div class="plot-container" style="margin-top:1.2rem;">%s</div>',
              img$ref_reads))
  }
  if (!is.null(img$ref_heat)) {
    parts <- c(parts,
      sprintf('<div class="plot-container" style="margin-top:1.2rem;">%s</div>',
              img$ref_heat))
  }
  if (nrow(ref_detail_data) > 0) {
    parts <- c(parts,
      '<h3 style="font-size:.9rem;margin-top:1.4rem;">Detailed Metrics</h3>',
      df_html(ref_detail_data, id = "ref-detail-table",
              cap = "Detailed Reference Metrics"))
  }
  if (length(parts) == 0) {
    '<p><em>No reference-selection files found in results directory.  '
    'Ensure SELECT_BEST_REFERENCE output (*_mapping_comparison.txt, '
    '*_detailed_metrics.txt) was copied into the report directory.</em></p>'
  } else {
    paste(parts, collapse = "\n")
  }
}
sec_ref_select <- card("ref-select", "Viral Reference Selection", ref_sel_body)

sec_table     <- card("sites-table",
                      "Integration Sites",
  paste0(
    '<p style="font-size:.8rem;color:var(--fg2);margin-bottom:.8rem;">',
    'viral_orientation and perl_strand both show + or -. ',
    'strand_concordance = agreement between the two calling methods. ',
    'Full data in the combined CSV files.</p>',
    df_html(integration_tbl)))

sections <- paste0(
  sec_overview, sec_reads,
  sec_intact,   sec_intact_a,
  sec_ori,      sec_episomal,
  sec_genome,   sec_chr,
  sec_genecov,  sec_hostgene,
  sec_len,      sec_div, sec_pid,
  sec_ref_select,
  sec_table)

# -----------------------------------------------------------------------------
# 13. Write HTML
# -----------------------------------------------------------------------------
html_doc <- paste0(
  '<!DOCTYPE html>\n',
  '<html lang="en" data-theme="light">\n',
  '<head>\n',
  '  <meta charset="UTF-8">\n',
  '  <meta name="viewport" content="width=device-width,initial-scale=1.0">\n',
  sprintf('  <title>HIV Viral Integration Report - %s</title>\n', run_name),
  css, "\n",
  js,  "\n",
  '</head>\n',
  '<body>\n',
  '<button class="theme-toggle" id="toggle-btn">Dark mode</button>\n',
  '<header>\n',
  '  <h1>HIV Viral Integration Report</h1>\n',
  sprintf('  <div class="subtitle">%s</div>\n', run_name),
  sprintf('  <div class="meta">Generated: %s | nf-viral-integration v0.2 | Connor S. Murray PhD, University of Louisville</div>\n',
          report_date),
  '</header>\n',
  '<div class="layout">\n',
  sprintf('<nav class="sidebar" aria-label="Sections"><h3>Sections</h3>%s</nav>\n',
          nav_html),
  '<div class="main-content">\n',
  sections,
  '</div>\n',
  '</div>\n',
  '<footer>Generated by <strong>nf-viral-integration</strong> &mdash; Connor S. Murray, PhD &mdash; University of Louisville SOM</footer>\n',
  '</body>\n',
  '</html>\n')

html_file <- paste0(output_prefix, "_report.html")
writeLines(html_doc, html_file)
cat(sprintf("\n  HTML report written to: %s\n", html_file))

# -----------------------------------------------------------------------------
# 14. FASTA output for MAFFT
#     Header: >sample|READ_ID|intactness|strand|chromosome|integration_site
# -----------------------------------------------------------------------------
cat("\nWriting FASTA for MAFFT alignment...\n")

if ("viral_sequence" %in% colnames(all_data)) {
  has_read <- "READ"              %in% colnames(all_data)
  has_int  <- "intactness"        %in% colnames(all_data)
  has_ori  <- "viral_orientation" %in% colnames(all_data)
  has_chr  <- "chromosome"        %in% colnames(all_data)
  has_site <- "integration_site"  %in% colnames(all_data)

  fasta_data <- all_data %>%
    filter(!is.na(viral_sequence),
           nchar(trimws(as.character(viral_sequence))) > 0) %>%
    mutate(
      read_id_col = if (has_read) as.character(READ)                               else as.character(row_number()),
      intact_col  = if (has_int)  gsub(" ", "_", as.character(intactness))         else "NA",
      strand_col  = if (has_ori)  replace_na(as.character(viral_orientation), "NA") else "NA",
      chr_col     = if (has_chr)  replace_na(as.character(chromosome),        "NA") else "NA",
      site_col    = if (has_site) replace_na(as.character(integration_site),  "NA") else "NA",
      fasta_id    = paste(sample, read_id_col, intact_col,
                          strand_col, chr_col, site_col, sep = "|"),
      fasta_id    = gsub("[[:space:]]", "_", fasta_id),
      seq_clean   = gsub("N", "", as.character(viral_sequence))
    ) %>%
    filter(nchar(seq_clean) >= 100)

  if (nrow(fasta_data) > 0) {
    fasta_lines <- unlist(lapply(seq_len(nrow(fasta_data)), function(i) {
      c(paste0(">", fasta_data$fasta_id[i]),
        fasta_data$seq_clean[i])
    }))
    fasta_out <- paste0(output_prefix, "_for_mafft.fasta")
    writeLines(fasta_lines, fasta_out)
    cat(sprintf("  %d sequences written to: %s\n",
                nrow(fasta_data), fasta_out))
    cat("  Suggested MAFFT command:\n")
    cat(sprintf("    mafft --auto --thread -1 %s > %s_aligned.fasta\n",
                fasta_out, output_prefix))
  } else {
    cat("  No sequences >= 100 bp available for FASTA output.\n")
  }
} else {
  cat("  NOTE: viral_sequence column not found; skipping FASTA output.\n")
}

cat(strrep("=", 70), "\n\n")
