#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# make_circos.R  (basic, ggplot2-only)
#
# Circos-style genome plot of HIV integration sites from the per-sample
# *.annotated.csv masterframe (one row per READ).
#
# Rows are collapsed to unique sites (chromosome + integration_pos + STRAND):
#   n_reads     = number of reads at the site
#   n_molecules = unique PCR_ID  (PCR-deduplicated support)
#   n_clones    = unique CLONE_ID
# Rows with chromosome == "HIV" / NA (no host flank) are dropped.
#
# Tracks (outside -> in):
#   chromosome ideogram | lollipops (height = log10 molecules, colour = strand
#   or sample) | optional gene labels for the top N sites
#
# Usage:
#   # one sample
#   Rscript make_circos.R --integrations 167-5.annotated.csv --genome t2t \
#       --out 167-5.circos.png
#   # several samples on one plot (comma-separated), coloured by sample
#   Rscript make_circos.R --integrations a.annotated.csv,b.annotated.csv \
#       --genome t2t --label "Project" --out project.circos.png
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--integrations", type = "character",
              help = "annotated.csv file(s), comma-separated"),
  make_option("--genome",       type = "character", default = "t2t",
              help = "t2t or hg38 [default %default]"),
  make_option("--chrom_sizes",  type = "character", default = NA,
              help = "Optional chrom.sizes (chrom<TAB>length) override"),
  make_option("--label",        type = "character", default = NA,
              help = "Plot title [default: sample name]"),
  make_option("--top_n",        type = "integer",   default = 10,
              help = "Label top N sites by molecule count with gene name [default %default]"),
  make_option("--out",          type = "character", default = "circos.png"),
  make_option("--width",        type = "double",    default = 8),
  make_option("--height",       type = "double",    default = 8),
  make_option("--dpi",          type = "integer",   default = 300)
)))
if (is.null(opt$integrations)) stop("--integrations is required")

# ---- Chromosome sizes -------------------------------------------------------
hg38_sizes <- c(
  chr1=248956422, chr2=242193529, chr3=198295559, chr4=190214555,
  chr5=181538259, chr6=170805979, chr7=159345973, chr8=145138636,
  chr9=138394717, chr10=133797422, chr11=135086622, chr12=133275309,
  chr13=114364328, chr14=107043718, chr15=101991189, chr16=90338345,
  chr17=83257441,  chr18=80373285,  chr19=58617616,  chr20=64444167,
  chr21=46709983,  chr22=50818468,  chrX=156040895,  chrY=57227415)
t2t_sizes <- c(
  chr1=248387328, chr2=242696752, chr3=201105948, chr4=193574945,
  chr5=182045439, chr6=172126628, chr7=160567428, chr8=146259331,
  chr9=150617247, chr10=134758134, chr11=135127769, chr12=133324548,
  chr13=113566686, chr14=101161492, chr15=99753195, chr16=96330374,
  chr17=84276897,  chr18=80542538,  chr19=61707364,  chr20=66210255,
  chr21=45090682,  chr22=51324926,  chrX=154259566,  chrY=62460029)

if (!is.na(opt$chrom_sizes)) {
  cs    <- fread(opt$chrom_sizes, header = FALSE, col.names = c("chrom", "length"),
                 select = 1:2)                       # auto-detects tab/comma
  sizes <- setNames(as.numeric(cs$length), cs$chrom)
} else {
  sizes <- switch(opt$genome, hg38 = hg38_sizes, t2t = t2t_sizes,
                  stop("Unknown --genome '", opt$genome, "'; use --chrom_sizes"))
}

# ---- Read masterframe(s) ----------------------------------------------------
files <- trimws(strsplit(opt$integrations, ",")[[1]])
keep  <- c("READ", "chromosome", "integration_pos", "STRAND",
           "CLONE_ID", "PCR_ID", "gene_name", "region")

dt <- rbindlist(lapply(files, function(f) {
  hdr <- names(fread(f, nrows = 0))
  x   <- fread(f, select = intersect(keep, hdr), na.strings = c("NA", ""))
  need <- c("chromosome", "integration_pos")
  if (!all(need %in% names(x)))
    stop(f, " is missing: ", paste(setdiff(need, names(x)), collapse = ", "))
  x[, sample_id := sub("\\.annotated\\.csv$", "", basename(f))]
  x
}), fill = TRUE)

for (col in setdiff(keep, names(dt))) dt[, (col) := NA]
dt[, integration_pos := suppressWarnings(as.numeric(integration_pos))]

n_reads_in <- nrow(dt)
dt <- dt[!is.na(chromosome) & chromosome != "HIV" & !is.na(integration_pos)]
off <- setdiff(unique(dt$chromosome), names(sizes))
if (length(off)) message("[circos] dropping sites on: ", paste(off, collapse = ", "))
dt <- dt[chromosome %in% names(sizes)]
message(sprintf("[circos] %d / %d reads have a host integration site",
                nrow(dt), n_reads_in))

# ---- Collapse reads -> unique sites -----------------------------------------
sites <- dt[, .(
  n_reads     = .N,
  n_molecules = if (all(is.na(PCR_ID)))   .N else uniqueN(PCR_ID,   na.rm = TRUE),
  n_clones    = if (all(is.na(CLONE_ID))) NA_integer_ else uniqueN(CLONE_ID, na.rm = TRUE),
  gene_name   = gene_name[!is.na(gene_name)][1],
  region      = region[!is.na(region)][1]
), by = .(sample_id, chromosome, integration_pos, STRAND)]
sites[, STRAND := fifelse(is.na(STRAND), "unknown", as.character(STRAND))]

samples  <- unique(sites$sample_id)
multi    <- length(samples) > 1
title    <- if (!is.na(opt$label)) opt$label else paste(samples, collapse = ", ")

if (nrow(sites) == 0) {
  p <- ggplot() + theme_void() +
    annotate("text", 0, 0, label = paste0(title, "\nNo host integration sites")) +
    xlim(-1, 1) + ylim(-1, 1)
  ggsave(opt$out, p, width = opt$width, height = opt$height, dpi = opt$dpi, bg = "white")
  message("[circos] no sites; wrote placeholder ", opt$out); quit(status = 0)
}

# ---- Genome -> angle --------------------------------------------------------
chr <- data.table(chrom = names(sizes), length = as.numeric(sizes))
gap <- 0.004 * sum(chr$length)                     # gap between chromosomes
chr[, start := cumsum(c(0, head(length + gap, -1)))]
chr[, end   := start + length]
total <- max(chr$end) + gap
theta <- function(g) pi / 2 - 2 * pi * g / total   # clockwise from 12 o'clock

sites[chr, on = .(chromosome = chrom), g := i.start + integration_pos]
sites[, th := theta(g)]

# ---- Geometry ---------------------------------------------------------------
R_ideo_in <- 1.00; R_ideo_out <- 1.07
R_base    <- 0.55; R_top      <- 0.95             # lollipop track

arc_poly <- function(a0, a1, r0, r1, n = 60) {
  a <- seq(a0, a1, length.out = n)
  data.table(x = c(r1 * cos(a), r0 * cos(rev(a))),
             y = c(r1 * sin(a), r0 * sin(rev(a))))
}
ideo <- chr[, arc_poly(theta(start), theta(end), R_ideo_in, R_ideo_out),
            by = .(chrom, idx = seq_len(nrow(chr)))]
ideo[, fill := ifelse(idx %% 2 == 0, "#8c8c8c", "#c8c8c8")]

chr[, `:=`(lab_th = theta((start + end) / 2), lab = sub("^chr", "", chrom))]

# track background + guide circles at 1, 10, 100, 1000 molecules
max_mol <- max(sites$n_molecules)
ymax    <- max(1, ceiling(log10(max_mol + 1)))
rad     <- function(n) R_base + (R_top - R_base) * log10(n + 1) / ymax
circle  <- function(r, n = 360) { a <- seq(0, 2 * pi, length.out = n)
                                  data.table(x = r * cos(a), y = r * sin(a)) }
guides_r <- rbindlist(lapply(10^(0:ymax), function(n)
  if (rad(n) <= R_top + 1e-6) circle(rad(n))[, n := n]))
guide_labs <- data.table(n = 10^(0:ymax))[, r := rad(n)][r <= R_top + 1e-6]

sites[, r := rad(n_molecules)]
sites[, colour := if (multi) sample_id else STRAND]
setorder(sites, n_molecules)                        # biggest drawn on top

strand_cols <- c(plus = "#1f77b4", minus = "#d62728", unknown = "#7f7f7f")

# ---- Plot -------------------------------------------------------------------
p <- ggplot() +
  geom_path(data = guides_r, aes(x, y, group = n), colour = "grey88", linewidth = 0.3) +
  geom_text(data = guide_labs, aes(0, r, label = n), size = 2.2,
            colour = "grey55", vjust = -0.3) +
  geom_polygon(data = ideo, aes(x, y, group = chrom, fill = I(fill))) +
  geom_text(data = chr, aes(1.13 * cos(lab_th), 1.13 * sin(lab_th), label = lab),
            size = 3, colour = "grey25") +
  geom_segment(data = sites,
               aes(x = R_base * cos(th), y = R_base * sin(th),
                   xend = r * cos(th), yend = r * sin(th), colour = colour),
               linewidth = 0.3, alpha = 0.6) +
  geom_point(data = sites,
             aes(x = r * cos(th), y = r * sin(th), colour = colour),
             size = 1.2, alpha = 0.9) +
  coord_fixed(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25)) +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, colour = "grey40"),
        legend.position = "bottom") +
  labs(title = title,
       subtitle = sprintf("%s  |  %d unique sites  |  %d reads  |  height = log10(molecules + 1)",
                          if (is.na(opt$chrom_sizes)) opt$genome else basename(opt$chrom_sizes),
                          nrow(sites), sum(sites$n_reads)),
       colour = if (multi) "Sample" else "Strand")

if (!multi) p <- p + scale_colour_manual(values = strand_cols)

# Label top-N sites with gene name (or chr:pos when intergenic)
if (opt$top_n > 0) {
  # greedy: take sites by molecule count, skip ones too close to a kept label
  cand <- sites[order(-n_molecules)]
  kept <- integer(0)
  for (i in seq_len(nrow(cand))) {
    if (length(kept) >= opt$top_n) break
    d <- abs(cand$th[i] - cand$th[kept]); d <- pmin(d, 2 * pi - d)
    if (all(d > 0.06)) kept <- c(kept, i)
  }
  top <- cand[kept]
  top[, lab := fifelse(is.na(gene_name),
                       sprintf("%s:%s", chromosome, format(integration_pos, big.mark = ",")),
                       gene_name)]
  top[, `:=`(lx = (R_base - 0.04) * cos(th), ly = (R_base - 0.04) * sin(th),
             hj = ifelse(cos(th) >= 0, 1, 0))]
  p <- p + geom_text(data = top, aes(lx, ly, label = lab, hjust = hj),
                     size = 2.4, fontface = "italic", colour = "grey20")
}

ggsave(opt$out, p, width = opt$width, height = opt$height, dpi = opt$dpi, bg = "white")
message("[circos] wrote ", opt$out)
