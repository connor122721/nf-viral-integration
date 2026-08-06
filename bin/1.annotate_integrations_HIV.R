#!/usr/bin/env Rscript
# Annotate HIV integration sites from a *_combined.csv against a host GTF
# Rscript annotate_integrations_HIV.R <combined.csv> <host.gtf> [rmsk.bed|NA] [out.tsv]

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges) 
  library(rtracklayer)
})

#setwd("X:BMG/smithlab/smrtcap_v4_07_2026/work/0c/4049e3c4f9af93390d7be1976e41aa")

## ---- args -----------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1 && nzchar(args[1])) args[1] else "151_15.combined.csv"
gtf_file <- if (length(args) >= 2 && nzchar(args[2])) args[2] else "host.gtf"
repeatmasker_bed <- if (length(args) >= 3 && !args[3] %in% c("", "NA")) args[3] else "repeatmasker/t2t.repeatmasker.bed.gz"
out_csv <- if (length(args) >= 4 && nzchar(args[4])) args[4] else
  sub("\\.csv$", "_annotated.tsv", basename(input_csv))
pos_col <- "POS"

first_non_na <- function(x) { x <- x[!is.na(x)]; if (length(x)) x[1] else NA }
pick_col     <- function(d, cands) { h <- cands[cands %in% names(d)]; if (length(h)) h[1] else NA_character_ }

## ---- 1. read input & derive site coordinates ------------------------------
df <- fread(input_csv, fill = TRUE)
df[, integration_pos := suppressWarnings(as.integer(get(pos_col)))]
df[, chromosome := sub("_[0-9]+$", "", CLONE_ID)]     # chr3_1 -> chr3

has_site <- grepl("^chr", df$chromosome) & !is.na(df$integration_pos)
cat(sprintf("  - %d / %d rows are host integration sites\n", sum(has_site), nrow(df)))

sites <- unique(df[has_site, .(chromosome, integration_pos)])
if (!nrow(sites)) {
  message("No host integration sites in ", basename(input_csv), " — writing pass-through.")
  write.csv(x = df, file = out_csv, quote = F, row.names = F); quit(save = "no", status = 0)
}
site_gr <- GRanges(sites$chromosome, IRanges(sites$integration_pos, width = 1L))
cat(sprintf("  - %d unique host sites\n", length(site_gr)))

## ---- 2. load GTF and DE-DUPLICATE -----------------------------------------
# Restrict the load to the site regions; fall back to a full import if the
# indexed 'which' path is unsupported for this file.
gtf <- tryCatch(import(gtf_file, which = site_gr),
                error = function(e) subsetByOverlaps(import(gtf_file), site_gr))
cat(sprintf("  - Loaded %d GTF features overlapping sites\n", length(gtf)))

gtf_df <- as.data.frame(gtf)

# Resolve attribute names that vary between GTF flavours (RefSeq vs Ensembl vs
# gffread/AGAT output) so we don't silently annotate with all-NA columns.
gene_name_col <- pick_col(gtf_df, c("gene_name", "gene", "Name", "gene_symbol"))
biotype_col <- pick_col(gtf_df, c("gene_biotype", "biotype", "gene_type"))
gtf_df$.gene_name <- if (!is.na(gene_name_col)) as.character(gtf_df[[gene_name_col]]) else NA_character_
gtf_df$.biotype <- if (!is.na(biotype_col)) as.character(gtf_df[[biotype_col]]) else NA_character_

# (2a) Drop exact-duplicate feature rows (same locus/strand/type/gene/transcript)
key_cols <- intersect(c("seqnames", "start", "end", "strand", "type",
                        "gene_id", "transcript_id"), names(gtf_df))
n_before <- nrow(gtf_df)
gtf_df <- gtf_df[!duplicated(gtf_df[, key_cols, drop = FALSE]), ]
cat(sprintf("  - GTF de-dup: removed %d duplicate feature rows (%d -> %d)\n",
            n_before - nrow(gtf_df), n_before, nrow(gtf_df)))

# Rebuild the (deduped) feature GRanges used for region typing
gtf <- makeGRangesFromDataFrame(gtf_df, keep.extra.columns = TRUE)

# (2b) Collapse to ONE gene-body range per gene locus (gene_id x seqnames x strand).
#      This guarantees the gene-overlap step can never emit a gene twice.
gene_bodies <- gtf_df %>%
  filter(!is.na(gene_id)) %>%
  group_by(gene_id, seqnames, strand) %>%
  summarise(gene_name = first_non_na(.gene_name),
            biotype = first_non_na(.biotype),
            start = min(start),
            end = max(end),
            .groups = "drop") %>%
  distinct(gene_id, seqnames, strand, .keep_all = TRUE)   # (2c) belt-and-braces
cat(sprintf("  - %d unique gene loci after collapse\n", nrow(gene_bodies)))

genes <- GRanges(gene_bodies$seqnames,
                 IRanges(gene_bodies$start, gene_bodies$end),
                 strand = gene_bodies$strand,
                 gene_id = gene_bodies$gene_id,
                 gene_name = gene_bodies$gene_name,
                 biotype = gene_bodies$biotype)

## ---- 3. per-site region typing (most-specific wins) -----------------------
ft <- as.character(gtf$type)
region <- rep("intergenic", length(site_gr))
region[overlapsAny(site_gr, gtf[ft %in% c("gene", "transcript", "intron")])] <- "intron"
region[overlapsAny(site_gr, gtf[ft == "exon"])] <- "exon"
region[overlapsAny(site_gr, gtf[ft %in% c("five_prime_UTR", "5UTR")])] <- "5UTR"
region[overlapsAny(site_gr, gtf[ft %in% c("three_prime_UTR", "3UTR")])] <- "3UTR"
region[overlapsAny(site_gr, gtf[ft == "CDS"])] <- "CDS"
region[overlapsAny(site_gr, gtf[ft == "start_codon"])] <- "start_codon"
region[overlapsAny(site_gr, gtf[ft == "stop_codon"])] <- "stop_codon"

site_ann <- data.table(chromosome = sites$chromosome,
                       integration_pos = sites$integration_pos,
                       region = region,
                       gene_name = NA_character_, gene_id = NA_character_,
                       biotype = NA_character_)

## ---- 4. gene overlap (unique() is redundant given (2b) but kept as a guard)-
ov <- findOverlaps(site_gr, genes)
if (length(ov)) {
  a <- data.table(i = queryHits(ov),
                  gene_name = genes$gene_name[subjectHits(ov)],
                  gene_id = genes$gene_id[subjectHits(ov)],
                  biotype = genes$biotype[subjectHits(ov)])
  a <- a[, .(gene_name = paste(unique(na.omit(gene_name)), collapse = ";"),
             gene_id = paste(unique(na.omit(gene_id)), collapse = ";"),
             biotype = paste(unique(na.omit(biotype)), collapse = ";")), by = i]
  site_ann[a$i, `:=`(gene_name = a$gene_name, gene_id = a$gene_id, biotype = a$biotype)]
  for (col in c("gene_name", "gene_id", "biotype"))
    site_ann[get(col) == "", (col) := NA_character_]
}
cat(sprintf("  - %d / %d sites annotated with a gene\n",
            sum(!is.na(site_ann$gene_name)), nrow(site_ann)))

## ---- 5. RepeatMasker overlap (optional) -----------------------------------
site_ann[, `:=`(in_repeat = 0L, repeat_name = NA_character_,
                repeat_class = NA_character_, repeat_family = NA_character_)]
if (!is.null(repeatmasker_bed) && file.exists(repeatmasker_bed)) {
  cat(sprintf("  - Loading RepeatMasker BED: %s\n", basename(repeatmasker_bed)))
  rmsk <- fread(repeatmasker_bed, header = FALSE, select = c(1:3, 10:12),
                col.names = c("chrom", "start", "end",
                              "rpt_name", "rpt_class", "rpt_family"))
  rmsk[rpt_family == "", rpt_family := NA_character_]
  rmsk_gr <- GRanges(rmsk$chrom, IRanges(rmsk$start + 1L, rmsk$end),   # BED 0-based -> 1-based
                     rpt_name = rmsk$rpt_name, rpt_class = rmsk$rpt_class,
                     rpt_family = rmsk$rpt_family)
  ro <- findOverlaps(site_gr, rmsk_gr)
  if (length(ro)) {
    r <- data.table(i = queryHits(ro),
                    rpt_name = rmsk_gr$rpt_name[subjectHits(ro)],
                    rpt_class = rmsk_gr$rpt_class[subjectHits(ro)],
                    rpt_family = rmsk_gr$rpt_family[subjectHits(ro)])
    r <- r[, .(repeat_name = paste(unique(na.omit(rpt_name)), collapse = ";"),
               repeat_class = paste(unique(na.omit(rpt_class)), collapse = ";"),
               repeat_family = paste(unique(na.omit(rpt_family)), collapse = ";")), by = i]
    site_ann[r$i, `:=`(in_repeat = 1L, repeat_name = r$repeat_name,
                       repeat_class = r$repeat_class, repeat_family = r$repeat_family)]
  }
  cat(sprintf("  - %d / %d sites overlap a repeat\n",
              sum(site_ann$in_repeat), nrow(site_ann)))
} else {
  cat("  - RepeatMasker BED not provided, skipping\n")
}

## ---- 6. join site annotation back onto every read row ---------------------
out <- data.table(merge(df, site_ann, by = c("chromosome", "integration_pos"),
             all.x = TRUE, sort = FALSE))
out <- out[!duplicated(out)]
setcolorder(out, c("READ", "chromosome", "integration_pos"))
write.csv(x = out, file = out_csv, quote = F, row.names = F)
cat(sprintf("wrote %s  (%d read rows, %d unique sites)\n",
            out_csv, nrow(out), nrow(site_ann)))
