#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# make_circos.R   (no-circlize edition)
#
# Renders a circos-style plot of viral integration sites using ggplot2 in
# polar coordinates. This avoids a dependency on the `circlize` package, which
# is NOT installed in either of the project SIFs (R_Genomics_v1_blast.def,
# Viral_Genomics_v5.def). It uses only packages that ARE in both:
#       ggplot2, data.table, optparse, dplyr/tidyverse
#
# Two modes (--mode):
#   sample   : one set of points, colored by viral_orientation if available.
#   project  : points colored by sample_id, with concentric persistence rings
#              underneath that show clone presence per timepoint.
#
# Two genome modes (--genome):
#   hg38     : chromosome lengths come from a built-in vector (GRCh38 sizes).
#   t2t      : chromosome lengths from a built-in vector (CHM13v2.0 sizes)
#              OR override with --chrom_sizes (UCSC chrom.sizes file:
#              two columns: chrom <TAB> length).
#
# Output: PNG file (default 1200×1200 at 150 dpi) suitable for embedding via
# <img src="..."> in the HTML reports.
#
# Usage:
#   Rscript make_circos.R --integrations sample.tsv --mode sample \
#       --genome t2t --label "Sample 1" --out sample.circos.png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--integrations", type = "character"),
  make_option("--mode",         type = "character", default = "sample"),
  make_option("--genome",       type = "character", default = "hg38"),
  make_option("--chrom_sizes",  type = "character", default = NA,
              help = "Optional UCSC chrom.sizes override"),
  make_option("--label",        type = "character", default = ""),
  make_option("--out",          type = "character", default = "circos.png"),
  make_option("--width",        type = "integer",   default = 1200),
  make_option("--height",       type = "integer",   default = 1200),
  make_option("--res",          type = "integer",   default = 150)
)))

# ---- Built-in chromosome sizes ---------------------------------------------
# GRCh38 primary assembly (chr1-22, X, Y, M)
hg38_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468,  chrX=156040895,  chrY=57227415,
  chrM=16569
)
# T2T-CHM13v2.0
t2t_sizes <- c(
  chr1=248387328, chr2=242696752, chr3=201105948, chr4=193574945,
  chr5=182045439, chr6=172126628, chr7=160567428, chr8=146259331,
  chr9=150617247, chr10=134758134, chr11=135127769, chr12=133324548,
  chr13=113566686, chr14=101161492, chr15=99753195, chr16=96330374,
  chr17=84276897,  chr18=80542538,  chr19=61707364,  chr20=66210255,
  chr21=45090682,  chr22=51324926,  chrX=154259566,  chrY=62460029,
  chrM=16569
)

load_sizes <- function(genome, chrom_sizes_file) {
  if (!is.na(chrom_sizes_file) && file.exists(chrom_sizes_file)) {
    df <- read.table(chrom_sizes_file, sep = "\t", header = FALSE,
                     col.names = c("chrom", "length"),
                     stringsAsFactors = FALSE)
    sz <- setNames(as.numeric(df$length), df$chrom); return(sz)
  }
  switch(genome, hg38 = hg38_sizes, t2t = t2t_sizes,
         stop("Unknown genome: ", genome,
              " (provide --chrom_sizes to override)"))
}

sizes <- load_sizes(opt$genome, opt$chrom_sizes)

# ---- Load integrations and align with chromosome ordering -------------------
dt <- fread(opt$integrations)

# Column aliasing: accept upstream-friendly names from the combined.csv.
# Prefer the integer column (integration_pos / canonical_pos) over
# "integration_site" which is a string like "chr1:40058490".
.alias <- list(
  chrom        = c("chrom",        "chromosome",       "Chrom",      "Chromosome"),
  position     = c("position",     "integration_pos",  "canonical_pos",
                   "Position",     "integration_site"),
  sample_id    = c("sample_id",    "sample",           "Sample",     "SampleID"),
  read_support = c("read_support", "reads_at_site",    "ReadSupport","Reads_At_Site")
)
for (canon in names(.alias)) {
  if (canon %in% colnames(dt)) next
  cands <- intersect(.alias[[canon]], colnames(dt))
  if (length(cands) > 0) setnames(dt, cands[1], canon)
}

required <- c("chrom", "position", "sample_id")
missing  <- setdiff(required, colnames(dt))
if (length(missing) > 0) stop("integrations missing column(s): ",
                              paste(missing, collapse = ", "),
                              "  (had: ", paste(colnames(dt), collapse = ", "), ")")

# Coerce + drop rows where the integration site couldn't be placed on the host
suppressWarnings(dt[, position := as.integer(position)])
dt[, chrom := as.character(chrom)]
n_before <- nrow(dt)
dt <- dt[!is.na(chrom) & chrom != "" & chrom != "NA"
         & !is.na(position) & position > 0]
if (nrow(dt) < n_before) {
  message(sprintf("[make_circos.R] dropped %d/%d rows with NA chrom/position (no host integration).",
                  n_before - nrow(dt), n_before))
}

if (nrow(dt) == 0) {
  # Render a polite "no data" PNG instead of failing the whole pipeline.
  png(opt$out, width = opt$width, height = opt$height, res = opt$res, type = "cairo")
  par(mar = c(0,0,0,0)); plot.new()
  text(0.5, 0.55, opt$label, cex = 1.1)
  text(0.5, 0.45, "No host-integration sites to plot.", cex = 0.9, col = "#666")
  dev.off()
  message("[make_circos.R] no plottable rows; wrote placeholder ", opt$out)
  quit(status = 0)
}

# Restrict to chromosomes we have lengths for; warn if we drop anything
unknown <- setdiff(unique(dt$chrom), names(sizes))
if (length(unknown) > 0) {
  message("[make_circos.R] dropping integrations on unknown chromosomes: ",
          paste(unknown, collapse = ", "))
  dt <- dt[chrom %in% names(sizes)]
}

if (nrow(dt) == 0) {
  png(opt$out, width = opt$width, height = opt$height, res = opt$res, type = "cairo")
  par(mar = c(0,0,0,0)); plot.new()
  text(0.5, 0.55, opt$label, cex = 1.1)
  text(0.5, 0.45,
       sprintf("All chromosomes unrecognised for genome '%s'.", opt$genome),
       cex = 0.9, col = "#666")
  dev.off()
  message("[make_circos.R] all chromosomes unrecognised; wrote placeholder ", opt$out)
  quit(status = 0)
}

# Build a global angular coordinate per integration:
#   theta = (cumulative_offset_of_chrom + position) / total_genome_length
#         expressed as a fraction of the circle.
chrom_order <- names(sizes)
chrom_dt <- data.table(
  chrom  = chrom_order,
  length = as.numeric(sizes[chrom_order])
)
chrom_dt[, start_offset := cumsum(c(0, head(length, -1)))]
chrom_dt[, end_offset   := start_offset + length]
chrom_dt[, mid_offset   := (start_offset + end_offset) / 2]
total_len <- sum(chrom_dt$length)

# Inter-chromosome gaps: drop a small fraction of the circle as gap so chroms
# are visually separated. We bake this in by inflating each "slot" length.
gap_frac     <- 0.005                        # 0.5% per gap
n_chrom      <- nrow(chrom_dt)
gap_total    <- gap_frac * n_chrom * total_len
slot_total   <- total_len + gap_total
chrom_dt[, slot_start := start_offset + (seq_len(.N) - 1) * gap_frac * total_len]
chrom_dt[, slot_end   := slot_start + length]
chrom_dt[, slot_mid   := (slot_start + slot_end) / 2]
to_theta <- function(off) (off / slot_total) * 2 * pi

dt <- merge(dt, chrom_dt[, .(chrom, slot_start)], by = "chrom", all.x = TRUE)
dt[, theta := to_theta(slot_start + position)]

# ---- Helper: chromosome ring geometry ---------------------------------------
make_ring <- function(chrom_dt, r0, r1) {
  rbindlist(lapply(seq_len(nrow(chrom_dt)), function(i) {
    n_seg <- 50
    th <- seq(to_theta(chrom_dt$slot_start[i]),
              to_theta(chrom_dt$slot_end[i]), length.out = n_seg)
    data.table(
      chrom = chrom_dt$chrom[i],
      x_in  = r0 * cos(th), y_in  = r0 * sin(th),
      x_out = r1 * cos(th), y_out = r1 * sin(th),
      group = i, theta = th
    )
  }))
}

# Ring polygons for chromosome ideogram
ring <- make_ring(chrom_dt, r0 = 0.92, r1 = 1.00)
ring_polys <- ring[, .(
  x = c(x_out, rev(x_in)), y = c(y_out, rev(y_in)),
  fill = chrom[1]
), by = group]

# Chromosome labels just outside the ring
chrom_dt[, lab_theta := to_theta(slot_mid)]
chrom_dt[, lab_x := 1.06 * cos(lab_theta)]
chrom_dt[, lab_y := 1.06 * sin(lab_theta)]
# Strip "chr" prefix for label brevity
chrom_dt[, label := sub("^chr", "", chrom)]

# Alternating ideogram fill (light/dark)
chrom_fill <- rep(c("#cccccc", "#888888"), length.out = nrow(chrom_dt))
names(chrom_fill) <- chrom_dt$chrom

# ---- Build base plot --------------------------------------------------------
base_p <- ggplot() +
  # Chromosome ring
  geom_polygon(data = ring_polys,
               aes(x = x, y = y, group = group, fill = fill),
               color = NA) +
  scale_fill_manual(values = chrom_fill, guide = "none") +
  # Chromosome labels
  geom_text(data = chrom_dt,
            aes(x = lab_x, y = lab_y, label = label),
            size = 3.0, color = "#333333") +
  coord_fixed() +
  theme_void() +
  theme(plot.title    = element_text(hjust = 0.5, size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 10,
                                     color = "#666666")) +
  ggtitle(opt$label,
          subtitle = sprintf("genome: %s   |   integrations: %d   |   chroms: %d",
                             opt$genome, nrow(dt), nrow(chrom_dt)))

# ---- Mode-specific point and ring layers -----------------------------------
add_points <- function(p, dt) {
  if (nrow(dt) == 0) return(p)

  # Radius for points: just inside the chromosome ring.
  pt_r <- 0.86
  if (opt$mode == "sample") {
    color_col <- if ("viral_orientation" %in% colnames(dt))
                   dt$viral_orientation else rep("integration", nrow(dt))
    pt_dt <- data.table(
      x      = pt_r * cos(dt$theta),
      y      = pt_r * sin(dt$theta),
      color  = color_col,
      size   = if ("read_support" %in% colnames(dt))
                 pmin(4, 1 + log10(dt$read_support + 1)) else 1.5
    )
    p <- p + geom_point(data = pt_dt,
                        aes(x = x, y = y, color = color, size = size),
                        alpha = 0.85) +
      scale_color_manual(values = c("+" = "#1f77b4", "-" = "#d62728",
                                    "*" = "#444444",
                                    "integration" = "#1f77b4"),
                         na.value = "#444444",
                         name = "viral orientation") +
      scale_size_identity()

  } else if (opt$mode == "project") {
    pt_dt <- data.table(
      x         = pt_r * cos(dt$theta),
      y         = pt_r * sin(dt$theta),
      sample_id = dt$sample_id,
      size      = if ("read_support" %in% colnames(dt))
                    pmin(4, 1 + log10(dt$read_support + 1)) else 1.5
    )
    p <- p + geom_point(data = pt_dt,
                        aes(x = x, y = y, color = sample_id, size = size),
                        alpha = 0.85) +
      scale_size_identity() +
      guides(color = guide_legend(title = "sample", ncol = 1,
                                   override.aes = list(size = 3)))

    # Persistence rings: one inner ring per timepoint
    if ("timepoint" %in% colnames(dt)) {
      tps <- sort(unique(stats::na.omit(dt$timepoint)))
      n   <- length(tps)
      if (n > 0) {
        # Reserve r in [0.55, 0.82] for persistence rings
        r_outer <- 0.82
        r_inner <- 0.55
        ring_h  <- (r_outer - r_inner) / max(n, 1)
        for (i in seq_along(tps)) {
          tp  <- tps[i]
          sub <- dt[timepoint == tp]
          if (nrow(sub) == 0) next
          rr0 <- r_outer - i * ring_h + 0.005
          rr1 <- r_outer - (i - 1) * ring_h - 0.005
          # tiny rectangles in polar
          tick_th_w <- 2 * pi * 0.0008  # tick width
          tick_dt <- data.table(
            theta = sub$theta,
            r0 = rr0, r1 = rr1, tw = tick_th_w
          )
          tick_polys <- rbindlist(lapply(seq_len(nrow(tick_dt)), function(j) {
            th0 <- tick_dt$theta[j] - tick_dt$tw[j] / 2
            th1 <- tick_dt$theta[j] + tick_dt$tw[j] / 2
            data.table(
              x = c(rr0 * cos(th0), rr1 * cos(th0),
                    rr1 * cos(th1), rr0 * cos(th1)),
              y = c(rr0 * sin(th0), rr1 * sin(th0),
                    rr1 * sin(th1), rr0 * sin(th1)),
              group = paste0("t", i, "_", j)
            )
          }))
          p <- p + geom_polygon(data = tick_polys,
                                aes(x = x, y = y, group = group),
                                fill = "#333333", color = NA)
          # Timepoint label at 9 o'clock
          p <- p + annotate("text",
                            x = -1.10, y = (rr0 + rr1) / 2,
                            label = as.character(tp),
                            size = 2.6, hjust = 1, color = "#555555")
        }
      }
    }
  }
  p
}

p <- add_points(base_p, dt)

# ---- Save -------------------------------------------------------------------
ggsave(opt$out, plot = p,
       width  = opt$width  / opt$res,
       height = opt$height / opt$res,
       dpi    = opt$res, units = "in", bg = "white")

message("[make_circos.R] wrote ", opt$out, " (",
        round(file.info(opt$out)$size / 1024, 1), " KB)")
