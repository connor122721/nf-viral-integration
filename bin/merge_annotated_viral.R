#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# merge_annotated_viral.R
#
# This script does a proper outer merge using data.table, preserving the
# column names and the canonical column order. Both files are expected to be
# CSVs WITHOUT a header (they are at the point of the join in the .nf script);
# the column names are passed via --annot_cols and --viral_cols so the merge
# can produce a correctly-labeled output.
#
# Usage:
#   Rscript merge_annotated_viral.R \
#       --annotated sampleX_annotated.csv \
#       --viral sampleX.viral_tmp.csv \
#       --annot_key_col 7 \
#       --viral_key_col 1 \
#       --out sampleX.combined.csv
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})
opt <- parse_args(OptionParser(option_list = list(
  make_option("--annotated", type = "character"),
  make_option("--viral", type = "character"),
  make_option("--annot_key_col", type = "integer", default = 7L,
              help = "1-based key column index in the annotated file"),
  make_option("--viral_key_col", type = "integer", default = 1L,
              help = "1-based key column index in the viral file"),
  make_option("--out", type = "character", default = "combined.csv")
)))
annot <- fread(opt$annotated,  fill = TRUE)
viral <- fread(opt$viral, sep = ",", fill = TRUE)
out <- merge(annot, viral, by = "READ")
fwrite(out, opt$out, sep = ",", quote = TRUE, na = "")
message(sprintf("[merge_annotated_viral.R] %d annotated + %d viral -> %d combined (%d cols)",
                nrow(annot), nrow(viral), nrow(out), ncol(out)))
