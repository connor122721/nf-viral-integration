#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# assign_clonal_ids.R
#
# Assigns persistent clonal IDs to viral integration sites so the same
# integration can be tracked across time-course samples.
#
# A "clone" is defined as a cluster of integration calls whose breakpoints sit
# within --window bp of one another on the same chromosome and (optionally) the
# same strand. The canonical clone ID is built from the chromosome and the
# integer median breakpoint of the cluster:
#       chrN_<position>            (strand-agnostic, default)
#       chrN_<position>_+/-        (with --use_strand)
#
# Column aliasing (NEW):
#   The combined.csv produced by the upstream merge uses upstream-friendly
#   names; we transparently map them to the canonical names this script needs:
#       sample          -> sample_id
#       chromosome      -> chrom
#       integration_pos -> position    (preferred — already an integer)
#       canonical_pos   -> position    (fallback)
#
# Samplesheet merge (NEW):
#   The combined.csv has only `sample`; it does NOT carry patient_id or
#   timepoint. When --samplesheet is provided we left-join it by the `sample`
#   column (or `sample_id` if that's the header) to pull in patient_id +
#   timepoint, which assign_clonal_ids needs for its longitudinal output.
#
# Inputs / outputs unchanged from the previous version.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(GenomicRanges)
})

option_list <- list(
  make_option("--integrations", type = "character",
              help = "Path to integrations CSV/TSV (required)"),
  make_option("--samplesheet",  type = "character", default = NA,
              help = "Optional sample sheet CSV/TSV with sample, timepoint, patient_id"),
  make_option("--window",       type = "integer",   default = 10L,
              help = "Clustering window in bp [default %default]"),
  make_option("--use_strand",   action = "store_true", default = FALSE,
              help = "Require same strand within cluster"),
  make_option("--out_long",     type = "character", default = "integrations_with_clonal_id.tsv"),
  make_option("--out_wide",     type = "character", default = "clonal_persistence_wide.tsv"),
  make_option("--out_summary",  type = "character", default = "clonal_summary.tsv")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(!is.null(opt$integrations))

# ---- Helper: rename a column to a canonical name if any alias is present ----
apply_alias <- function(d, canon, candidates) {
  if (canon %in% colnames(d)) return(d)
  hit <- intersect(candidates, colnames(d))
  if (length(hit) > 0) setnames(d, hit[1], canon)
  d
}

# ---- Load integrations + alias columns --------------------------------------
intg <- fread(opt$integrations)

intg <- apply_alias(intg, "sample_id", c("sample_id", "sample", "Sample", "SampleID"))
intg <- apply_alias(intg, "chrom",     c("chrom", "chromosome", "Chrom", "Chromosome"))
intg <- apply_alias(intg, "position",  c("position", "integration_pos", "canonical_pos",
                                         "Position"))
intg <- apply_alias(intg, "read_support", c("read_support", "reads_at_site",
                                            "ReadSupport", "Reads_At_Site"))

required_cols <- c("sample_id", "chrom", "position")
missing <- setdiff(required_cols, colnames(intg))
if (length(missing) > 0) {
  stop("integrations file is missing required column(s): ",
       paste(missing, collapse = ", "),
       "  (after aliasing; available: ",
       paste(colnames(intg), collapse = ", "), ")")
}

if (!"strand" %in% colnames(intg)) intg[, strand := "*"]
suppressWarnings(intg[, position := as.integer(position)])

# ---- Drop any pre-existing columns we are about to re-derive ---------------
# The upstream combined.csv already carries `canonical_pos` (and sometimes
# `clonal_id` / `clone_id`). If we left them in place, the merge below would
# silently rename them to canonical_pos.x / canonical_pos.y and the later
# setorder() would fail with "columns are not in the data.table".
for (col in c("canonical_pos", "clonal_id", "cluster_idx")) {
  if (col %in% colnames(intg)) {
    intg[, (col) := NULL]
  }
}

# Drop rows with NA chrom/position (no host integration to track)
n_before <- nrow(intg)
intg <- intg[!is.na(chrom) & chrom != "" & chrom != "NA"
             & !is.na(position) & position > 0]
if (nrow(intg) < n_before) {
  message(sprintf("[assign_clonal_ids.R] dropped %d/%d rows with NA chrom/position",
                  n_before - nrow(intg), n_before))
}

if (nrow(intg) == 0) {
  message("[assign_clonal_ids.R] no rows with mappable host coordinates; ",
          "writing empty outputs and exiting cleanly.")
  fwrite(data.table(clonal_id = character(0)),                opt$out_long,    sep = "\t")
  fwrite(data.table(clonal_id = character(0)),                opt$out_wide,    sep = "\t")
  fwrite(data.table(clonal_id = character(0), n_samples = 0L), opt$out_summary, sep = "\t")
  quit(status = 0)
}

# ---- Optional samplesheet merge for patient_id + timepoint ------------------
if (!is.na(opt$samplesheet) && file.exists(opt$samplesheet)) {
  ss <- fread(opt$samplesheet)
  ss <- apply_alias(ss, "sample_id", c("sample_id", "sample", "Sample", "SampleID"))
  ss_keep <- intersect(c("sample_id", "patient_id", "timepoint"), colnames(ss))
  if (!"sample_id" %in% ss_keep) {
    stop("samplesheet has no 'sample' or 'sample_id' column;  available: ",
         paste(colnames(ss), collapse = ", "))
  }
  ss <- unique(ss[, ..ss_keep])
  message(sprintf("[assign_clonal_ids.R] merging samplesheet (%d unique samples, cols: %s)",
                  nrow(ss), paste(ss_keep, collapse = ", ")))
  intg <- merge(intg, ss, by = "sample_id", all.x = TRUE,
                suffixes = c("", ".ss"))
  # If integrations already had timepoint/patient_id (mostly empty) prefer the
  # samplesheet value when ours is NA.
  for (col in c("timepoint", "patient_id")) {
    ss_col <- paste0(col, ".ss")
    if (ss_col %in% colnames(intg)) {
      intg[is.na(get(col)) | get(col) == "", (col) := get(ss_col)]
      intg[, (ss_col) := NULL]
    }
  }
}

if (!"timepoint" %in% colnames(intg))   intg[, timepoint  := NA_character_]
if (!"patient_id" %in% colnames(intg))  intg[, patient_id := NA_character_]

# ---- Cluster integrations within window using GRanges::reduce ---------------
gr <- GRanges(
  seqnames = intg$chrom,
  ranges   = IRanges(start = pmax(1L, intg$position - opt$window),
                     end   = intg$position + opt$window),
  strand   = if (opt$use_strand) intg$strand else "*"
)

clusters <- reduce(gr, ignore.strand = !opt$use_strand)
hits     <- findOverlaps(gr, clusters, ignore.strand = !opt$use_strand)
intg[, cluster_idx := subjectHits(hits)[match(seq_len(.N), queryHits(hits))]]

clone_pos <- intg[, .(canonical_pos = as.integer(round(median(position)))),
                  by = .(cluster_idx, chrom)]
if (opt$use_strand) {
  clone_pos[, strand := intg[match(cluster_idx, intg$cluster_idx), strand]]
  clone_pos[, clonal_id := sprintf("%s_%d_%s", chrom, canonical_pos, strand)]
} else {
  clone_pos[, clonal_id := sprintf("%s_%d", chrom, canonical_pos)]
}

intg <- merge(intg, clone_pos[, .(cluster_idx, clonal_id, canonical_pos)],
              by = "cluster_idx", all.x = TRUE)

# Reorder columns
front <- c("clonal_id", "patient_id", "sample_id", "timepoint",
           "chrom", "position", "canonical_pos", "strand")
front <- intersect(front, colnames(intg))
rest  <- setdiff(colnames(intg), c(front, "cluster_idx"))
intg  <- intg[, c(front, rest), with = FALSE]

setorder(intg, chrom, canonical_pos, sample_id)

fwrite(intg, opt$out_long, sep = "\t")

# ---- Wide persistence matrix -----------------------------------------------
val_col <- if ("read_support" %in% colnames(intg)) "read_support" else NULL
if (!is.null(val_col)) {
  wide <- dcast(intg, clonal_id + chrom + canonical_pos ~ timepoint,
                value.var = val_col, fun.aggregate = sum, fill = 0L)
} else {
  wide <- dcast(intg, clonal_id + chrom + canonical_pos ~ timepoint,
                value.var = "sample_id", fun.aggregate = length, fill = 0L)
}
fwrite(wide, opt$out_wide, sep = "\t")

# ---- Per-clone summary ------------------------------------------------------
summary_dt <- intg[, .(
  n_samples       = uniqueN(sample_id),
  n_timepoints    = uniqueN(timepoint),
  first_timepoint = if (all(is.na(timepoint))) NA_character_ else min(timepoint, na.rm = TRUE),
  last_timepoint  = if (all(is.na(timepoint))) NA_character_ else max(timepoint, na.rm = TRUE),
  total_reads     = if (!is.null(val_col)) sum(get(val_col), na.rm = TRUE) else NA_integer_,
  patient_id      = paste(unique(stats::na.omit(patient_id)), collapse = ",")
), by = .(clonal_id, chrom, canonical_pos)]
setorder(summary_dt, -n_timepoints, -n_samples)
fwrite(summary_dt, opt$out_summary, sep = "\t")

message("[assign_clonal_ids.R] ", nrow(intg), " integrations collapsed into ",
        nrow(summary_dt), " clones (window = ", opt$window, " bp).")
