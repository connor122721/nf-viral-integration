#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# merge_repeatmasker_annotation.R
#
# Joins the bedtools-closest output produced by REPEATMASKER_OVERLAP onto the
# main integrations table so each integration row gains:
#
#       repeat_name        e.g. AluY
#       repeat_class       parsed from RepeatMasker name (LINE/SINE/LTR/...)
#       repeat_family      parsed sub-class
#       distance_to_repeat 0 = integration sits *inside* a repeat
#       in_repeat          TRUE if distance == 0
#
# Usage:
#   Rscript merge_repeatmasker_annotation.R \
#       --integrations integrations_with_clonal_id.tsv \
#       --rmsk_overlap sampleX.t2t.integrations.repeatmasker.tsv \
#       --genome_label t2t \
#       --out integrations.t2t_repeats.tsv
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--integrations",  type = "character"),
  make_option("--rmsk_overlap",  type = "character"),
  make_option("--genome_label",  type = "character", default = "t2t"),
  make_option("--out",           type = "character",
              default = "integrations.with_repeats.tsv")
)))

intg <- fread(opt$integrations)
ovl  <- fread(opt$rmsk_overlap)

# UCSC RepeatMasker bigBed packs class/family into the `name` field as
# "name#class/family" (e.g. AluY#SINE/Alu). Some tracks use only "name".
ovl[, c("repeat_name_raw", "repeat_class_family") :=
      tstrsplit(rmsk_name, "#", fixed = TRUE, fill = NA_character_)]
ovl[, c("repeat_class", "repeat_family") :=
      tstrsplit(repeat_class_family, "/", fixed = TRUE, fill = NA_character_)]
ovl[is.na(repeat_name_raw), repeat_name_raw := rmsk_name]

ovl_min <- ovl[, .(
  chrom              = int_chrom,
  position           = as.integer(int_start),       # bed half-open
  repeat_name        = repeat_name_raw,
  repeat_class,
  repeat_family,
  distance_to_repeat = as.integer(distance_to_repeat),
  in_repeat          = distance_to_repeat == 0L
)]

# Suffix columns by genome so the same integration row can carry both T2T and
# HG38 repeat context simultaneously (they almost always disagree).
sfx <- paste0("_", opt$genome_label)
setnames(ovl_min,
         old = c("repeat_name", "repeat_class", "repeat_family",
                 "distance_to_repeat", "in_repeat"),
         new = paste0(c("repeat_name", "repeat_class", "repeat_family",
                        "distance_to_repeat", "in_repeat"), sfx))

intg[, position := as.integer(position)]
out <- merge(intg, ovl_min, by = c("chrom", "position"),
             all.x = TRUE, sort = FALSE)
fwrite(out, opt$out, sep = "\t")

message("[merge_repeatmasker_annotation.R] joined ",
        sum(!is.na(out[[paste0("repeat_name", sfx)]])),
        " / ", nrow(out), " integrations to ", opt$genome_label, " repeats.")
