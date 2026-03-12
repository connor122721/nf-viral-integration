#!/usr/bin/env Rscript
#
# Viral Integration Orientation and Sequence Analysis 
#
# Usage:
#   Rscript simple_annotate_bam_v2.R \
#       <input.csv> <unmasked.fasta> <hiv_ref.fasta> <gtf> \
#       <output_prefix> [iteration1.sam]
#
# The optional 6th argument is the path to the iteration-1 SAM file
# (e.g. sampleID.dedup.1.sam).  When supplied orientation is read directly
# from FLAG bits.  When absent the k-mer fallback handles all sequences.
#
# Output
# ------
#   <prefix>_annotated.csv          – per-read annotated table
#   <prefix>_viral_sequences.fasta  – extracted viral sequences (FASTA)
#   <prefix>_similarity_matrix.csv  – pairwise % identity matrix
#   <prefix>_report.txt             – plain-text summary

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
  library(ape)
})

# ============================================================================
# 0. Arguments
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript simple_annotate_bam_v2.R",
      "<input.csv> <unmasked.fasta> <hiv_ref.fasta> <gtf>",
      "<output_prefix> [iteration1.sam]\n")
  quit(status = 1)
}

input_csv <- args[1]
input_fasta <- args[2]
hiv_ref_fasta <- args[3]
gtf_file <- args[4]
output_prefix <- args[5]
iter1_sam <- if (length(args) >= 6) args[6] else NULL

# Test
# input_csv="X:/nf-viral-integration_t2t/test_results/04_final_results/92UG_029.subset.csv"; input_fasta="X:/nf-viral-integration_t2t/test_results/unmasked_sequences/92UG_029.subset.unmasked.fa"; hiv_ref_fasta="X:/nf-viral-integration_t2t/data/hiv_genome_panal/A1.UG..UG031.AB098330.fasta"; gtf_file="X:/nf-viral-integration_t2t/test_results/genome_files/host.gtf"; output_prefix="test1"; iter1_sam="X:/nf-viral-integration_t2t/test_results/02_iterative_masking/92UG_029.subset/92UG_029.subset.dedup.1.sam"

cat("\n", strrep("=", 70), "\n")
cat("VIRAL ORIENTATION AND GENOMIC ANNOTATION ANALYSIS (v3 – SAM/k-mer)\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# STEP 1: Read CSV data
# ============================================================================
cat("Step 1: Reading integration data...\n")
df <- read_csv(input_csv, show_col_types = FALSE) %>%
  filter(INSERT_LEN >= 200)
cat(sprintf("  - Total reads with INSERT_LEN >= 200: %d\n", nrow(df)))

# ============================================================================
# STEP 2: Read FASTA sequences
# ============================================================================
cat("\nStep 2: Reading FASTA sequences...\n")
fasta <- readDNAStringSet(input_fasta)
cat(sprintf("  - Loaded %d sequences\n", length(fasta)))

extract_viral <- function(seq_string) gsub("N", "", seq_string)

seq_lookup <- tibble(
  read_id = names(fasta),
  full_sequence = as.character(fasta)) %>%
  mutate(viral_sequence = map_chr(full_sequence, extract_viral),
         viral_seq_length = nchar(viral_sequence),
         read_id_base = str_remove(read_id, "/\\d+$"))

# ============================================================================
# STEP 3: Read HIV reference genome
# ============================================================================
cat("\nStep 3: Reading HIV reference genome...\n")
hiv_ref <- readDNAStringSet(hiv_ref_fasta)
cat(sprintf("  - Reference: %s\n", names(hiv_ref)[1]))
cat(sprintf("  - Length   : %d bp\n", width(hiv_ref)[1]))

# ============================================================================
# STEP 4: Determine viral orientation (SAM-flag approach + k-mer fallback)
# ============================================================================
cat("\nStep 4: Determining viral orientation...\n")

seq_lookup <- seq_lookup %>%
  mutate(orientation = NA_character_,
         strand = NA_character_,
         alignment_score = NA_real_,
         ref_start = NA_integer_,
         ref_end = NA_integer_,
         percent_identity = NA_real_)

# ---------------------------------------------------------------------------
# 4a. Parse SAM FLAG bits when iteration-1 SAM is available
# ---------------------------------------------------------------------------
sam_orient <- tibble()

if (!is.null(iter1_sam) && file.exists(iter1_sam)) {
  cat(sprintf("  - Parsing SAM flags from: %s\n", basename(iter1_sam)))

  sam_lines <- readLines(iter1_sam, warn = FALSE)
  sam_lines <- sam_lines[!startsWith(sam_lines, "@")] # drop header

  if (length(sam_lines) > 0) {
    # Read all columns so that optional fields (NM tag) are available.
    # Columns 1-11 are mandatory SAM fields; optional tags follow from col 12 onward.
    sam_dt <- fread(text = paste(sam_lines, collapse = "\n"),
                    sep = "\t", header = FALSE, fill = TRUE)
    setnames(sam_dt,
             seq_len(min(ncol(sam_dt), 11L)),
             c("qname","flag","rname","pos","mapq","cigar",
               "rnext","pnext","tlen","seq","qual")[seq_len(min(ncol(sam_dt), 11L))])

    # ---------------------------------------------------------------------------
    # CIGAR helpers
    # cigar_ref_len : number of reference bases consumed (M/D/N/X/= ops)
    # Used to compute percent identity against the viral reference section
    # that was actually aligned.
    # ---------------------------------------------------------------------------
    cigar_ref_len <- function(cig) {
      if (is.na(cig) || cig == "*") return(NA_integer_)
      ops    <- base::regmatches(cig, base::gregexpr("[0-9]+[MIDNSHPX=]", cig))[[1]]
      widths <- as.integer(sub("[A-Z=]$", "", ops))
      types  <- sub("^[0-9]+",  "", ops)
      sum(widths[types %in% c("M","D","N","X","=")])
    }

    # NM tag: edit distance to the reference (mismatches + indel bases)
    extract_nm <- function(row_fields) {
      nm_tag <- row_fields[startsWith(row_fields, "NM:i:")]
      if (length(nm_tag) == 0L) return(NA_integer_)
      as.integer(sub("^NM:i:", "", nm_tag[1L]))
    }

    # Build optional-field matrix (cols 12+) as character
    opt_cols <- if (ncol(sam_dt) >= 12L) {
      as.matrix(sam_dt[, 12:ncol(sam_dt), with = FALSE])
    } else {
      matrix(character(0), nrow = nrow(sam_dt), ncol = 0L)
    }

    sam_tbl <- sam_dt %>%
      as_tibble() %>%
      filter(rname != "*",           # mapped reads only
             !bitwAnd(flag, 0x4),    # not unmapped
             !bitwAnd(flag, 0x100),  # not secondary
             !bitwAnd(flag, 0x800))  # not supplementary

    if (nrow(sam_tbl) > 0) {
      # Row indices that survive the filter (needed to index opt_cols)
      keep_idx <- which(
        sam_dt$rname != "*" &
        !bitwAnd(sam_dt$flag, 0x4) &
        !bitwAnd(sam_dt$flag, 0x100) &
        !bitwAnd(sam_dt$flag, 0x800))

      ref_aln_len <- vapply(sam_tbl$cigar,  cigar_ref_len, integer(1L))
      nm_vals <- vapply(seq_along(keep_idx), function(i) {
        extract_nm(as.character(opt_cols[keep_idx[i], ]))
      }, integer(1L))

      sam_orient <- sam_tbl %>%
        mutate(
          strand = if_else(bitwAnd(flag, 0x10) != 0, "-", "+"),
          orientation = strand,
          ref_start = as.integer(pos),
          ref_end = as.integer(pos) + ref_aln_len - 1L,
          # percent identity: (aligned_bases - edit_distance) / aligned_bases * 100
          percent_identity = ifelse(
            !is.na(ref_aln_len) & ref_aln_len > 0 & !is.na(nm_vals),
            round(100 * (ref_aln_len - nm_vals) / ref_aln_len, 2),
            NA_real_)) %>%
        group_by(qname) %>%
        slice_max(mapq, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        dplyr::select(read_id = qname, strand, orientation,
                      ref_start, ref_end, percent_identity)

      cat(sprintf("  - Orientation from SAM flags: %d reads\n", nrow(sam_orient)))
      pid_ok <- sum(!is.na(sam_orient$percent_identity))
      cat(sprintf("  - Percent identity computed : %d reads (NM tag present)\n", pid_ok))
    } else {
      sam_orient <- tibble()
    }
  }
}

# Merge SAM-derived orientation + percent_identity into seq_lookup
if (nrow(sam_orient) > 0) {
  seq_lookup <- seq_lookup %>%
    left_join(sam_orient,
              by = c("read_id_base" = "read_id"),
              suffix = c("", "_sam")) %>%
    mutate(strand           = coalesce(strand_sam,           strand),
           orientation      = coalesce(orientation_sam,      orientation),
           ref_start        = coalesce(ref_start_sam,        ref_start),
           ref_end          = coalesce(ref_end_sam,          ref_end),
           percent_identity = coalesce(percent_identity_sam, percent_identity)) %>%
    dplyr::select(-ends_with("_sam"))
}

# ---------------------------------------------------------------------------
# 4b. K-mer fallback for reads not resolved by SAM
# ---------------------------------------------------------------------------
# Build a compact set of 21-mers from the first 1 kb (5′) and last 1 kb (3′)
# of the reference.  We sample every 50 bp to get ~40 kmers per end.
# A read is assigned the orientation whose end gives more k-mer hits.

cat("  - Building k-mer index from HIV reference for fallback orientation...\n")
ref_seq <- as.character(hiv_ref[[1]])
ref_len <- nchar(ref_seq)
K <- 21L
SAMPLE_STEP <- 50L
KMER_REGION <- min(1000L, floor(ref_len / 3))

build_kmer_set <- function(region_seq) {
  starts <- seq(1L, nchar(region_seq) - K, by = SAMPLE_STEP)
  kmers <- substring(region_seq, starts, starts + K - 1L)
  unique(kmers[!grepl("[^ACGT]", kmers)])
}

ref_5prime_region <- substr(ref_seq, 1L, KMER_REGION)
ref_3prime_region <- substr(ref_seq, ref_len - KMER_REGION + 1L, ref_len)

kmers_5 <- build_kmer_set(ref_5prime_region)
kmers_3 <- build_kmer_set(ref_3prime_region)

rev_comp <- function(s) {
  as.character(reverseComplement(DNAString(s)))
}
kmers_5_rc <- vapply(kmers_5, rev_comp, character(1))
kmers_3_rc <- vapply(kmers_3, rev_comp, character(1))

count_kmer_hits <- function(seq_str, kmer_set) {
  if (nchar(seq_str) < K) return(0L)
  sum(sapply(kmer_set, function(km) {
    grepl(km, seq_str, fixed = TRUE)
  }))
}

needs_fallback <- seq_lookup %>%
  filter(is.na(orientation),
         !is.na(viral_sequence),
         viral_sequence != "",
         nchar(viral_sequence) >= 500)

cat(sprintf("  - K-mer fallback for %d unresolved sequences\n", nrow(needs_fallback)))

if (nrow(needs_fallback) > 0) {
  fallback_results <- needs_fallback %>%
    rowwise() %>%
    mutate(
      hits_5fwd = count_kmer_hits(viral_sequence, kmers_5),
      hits_3fwd = count_kmer_hits(viral_sequence, kmers_3),
      hits_5rev = count_kmer_hits(viral_sequence, kmers_5_rc),
      hits_3rev = count_kmer_hits(viral_sequence, kmers_3_rc),
      score_fwd = hits_5fwd + hits_3fwd,
      score_rev = hits_5rev + hits_3rev,
      score_diff = abs(score_fwd - score_rev),
      # Require at least 2 more k-mer hits on one strand to make a call;
      # reads with equal or near-zero scores get NA rather than silently
      # piling into "+" and inflating the forward count.
      orientation = case_when(
        pmax(score_fwd, score_rev) == 0 ~ NA_character_,
        score_diff <= 1 ~ NA_character_,
        score_fwd > score_rev ~ "+",
        score_rev > score_fwd ~ "-",
        TRUE ~ NA_character_),
      strand = orientation,
      alignment_score = as.numeric(pmax(score_fwd, score_rev))) %>%
    ungroup() %>%
    dplyr::select(read_id, orientation, strand, alignment_score)

  seq_lookup <- seq_lookup %>%
    rows_update(fallback_results, by = "read_id", unmatched = "ignore")

  cat(sprintf("  - K-mer orientation resolved: %d reads\n",
              sum(!is.na(fallback_results$orientation))))
}

ori_table <- seq_lookup %>% count(orientation) %>% arrange(desc(n))
cat("  - Orientation summary:\n")
for (i in seq_len(nrow(ori_table))) {
  cat(sprintf("      %-12s: %d\n", ori_table$orientation[i], ori_table$n[i]))
}

# ============================================================================
# STEP 5: Merge sequence data with CSV
# ============================================================================
cat("\nStep 5: Merging sequence and orientation data...\n")

df <- data.table(df %>%
  mutate(read_match = case_when(
    READ %in% seq_lookup$read_id ~ READ,
    paste0(READ, "/0") %in% seq_lookup$read_id ~ paste0(READ, "/0"),
    str_remove(READ, "/0$") %in% seq_lookup$read_id_base ~ {
      idx <- match(str_remove(READ, "/0$"), seq_lookup$read_id_base)
      seq_lookup$read_id[idx]
    },
    TRUE ~ NA_character_)) %>%
  left_join(
    seq_lookup %>%
      dplyr::select(read_id, viral_sequence, viral_seq_length,
                    viral_orientation = orientation,
                    viral_strand = strand,
                    alignment_score,
                    ref_start, ref_end, percent_identity),
    by = c("read_match" = "read_id")) %>%
  dplyr::select(-read_match))

# ============================================================================
# STEP 6: Parse genomic locations
# ============================================================================
cat("\nStep 6: Processing integration sites...\n")

parse_location <- function(loc_string) {
  if (is.na(loc_string) || loc_string == "")
    return(tibble(chrom = NA_character_, start = NA_integer_, end = NA_integer_))
  tryCatch({
    parts <- str_split(loc_string, ":", n = 2)[[1]]
    chrom <- parts[1]
    coords <- str_split(parts[2], "-")[[1]]
    start <- as.integer(coords[1])
    end <- as.integer(coords[2])
    tibble(chrom = chrom, start = min(start, end), end = max(start, end))
  }, error = function(e)
    tibble(chrom = NA_character_, start = NA_integer_, end = NA_integer_))
}

get_integration_site <- function(left_flank, right_flank) {
  for (fl in c(left_flank, right_flank)) {
    if (!is.na(fl) && fl != "") {
      loc <- parse_location(fl)
      if (!is.na(loc$chrom))
        return(sprintf("%s:%d-%d", loc$chrom, loc$start, loc$end))
    }
  }
  NA_character_
}

get_viral_region <- function(insert_string) {
  if (is.na(insert_string) || insert_string == "") return(NA_character_)
  tryCatch({
    parts <- str_split(insert_string, ":")[[1]]
    if (length(parts) >= 2) {
      ref <- paste(parts[-length(parts)], collapse = ":")
      coords <- parts[length(parts)]
      cp <- str_split(coords, "-")[[1]]
      p1 <- as.integer(cp[1]); p2 <- as.integer(cp[2])
      sprintf("%s:%d-%d", ref, min(p1, p2), max(p1, p2))
    } else insert_string
  }, error = function(e) insert_string)
}

df <- df %>%
  mutate(integration_site = map2_chr(LEFT_FLANK, RIGHT_FLANK, get_integration_site),
         viral_region = map_chr(INSERT, get_viral_region),
         chromosome = map_chr(integration_site, ~ {
      if (is.na(.x)) { NA_character_ } else { str_split(.x, ":", n = 2)[[1]][1] }
    }),
    sample = sub("\\.csv$", "", basename(input_csv)))

# ============================================================================
# STEP 7: GTF annotation of integration sites
# ============================================================================
cat("\nStep 7: Loading GTF and annotating genomic features...\n")

gtf <- import(gtf_file)
cat(sprintf("  - Loaded %d features from GTF\n", length(gtf)))
gtf_df <- as.data.frame(gtf)

genes_df <- gtf_df %>%
  filter(!is.na(gene_id)) %>%
  group_by(seqnames, gene_id, gene_name, strand) %>%
  summarise(start = min(start), end = max(end), .groups = "drop")

genes <- GRanges(
  seqnames = genes_df$seqnames,
  ranges = IRanges(start = genes_df$start, end = genes_df$end),
  strand = genes_df$strand,
  gene_id = genes_df$gene_id,
  gene_name = genes_df$gene_name)

integration_data <- df %>%
  filter(!is.na(integration_site)) %>%
  mutate(row_idx = row_number()) %>%
  rowwise() %>%
  mutate(loc_data = list(parse_location(integration_site))) %>%
  ungroup() %>%
  unnest(loc_data) %>%
  filter(!is.na(chrom))

if (nrow(integration_data) > 0) {
  integration_gr <- GRanges(
    seqnames = integration_data$chrom,
    ranges = IRanges(start = integration_data$start,
                     end = integration_data$end),
    original_row = integration_data$row_idx)

  overlaps <- findOverlaps(integration_gr, genes)

  if (length(overlaps) > 0) {
    ann <- tibble(
      original_row = integration_gr$original_row[queryHits(overlaps)],
      gene_name = genes$gene_name[subjectHits(overlaps)],
      gene_id = genes$gene_id[subjectHits(overlaps)])

    ann_agg <- ann %>%
      group_by(original_row) %>%
      summarize(gene_name = paste(unique(na.omit(gene_name)), collapse = ";"),
                gene_id = paste(unique(na.omit(gene_id)), collapse = ";"),
                .groups = "drop") %>%
      mutate(gene_name = if_else(gene_name == "", NA_character_, gene_name),
             gene_id = if_else(gene_id == "", NA_character_, gene_id))

    df <- df %>%
      mutate(row_idx = row_number()) %>%
      left_join(ann_agg, by = c("row_idx" = "original_row")) %>%
      dplyr::select(-row_idx)

    cat(sprintf("  - Annotated %d reads with gene information\n",
                sum(!is.na(df$gene_name))))
  } else {
    cat("  - No overlaps found between integration sites and genes\n")
    df <- df %>% mutate(gene_name = NA_character_, gene_id = NA_character_)
  }
} else {
  cat("  - No valid integration sites to annotate\n")
  df <- df %>% mutate(gene_name = NA_character_, gene_id = NA_character_)
}

# ============================================================================
# STEP 8: Sequence similarity matrix (k-mer based)
# ============================================================================
cat("\nStep 8: Sequence similarity analysis (k-mer Jaccard)...\n")

unique_seqs <- df %>%
  filter(!is.na(viral_sequence), viral_seq_length >= 100) %>%
  distinct(viral_sequence, .keep_all = TRUE)

cat(sprintf("  - Unique viral sequences (>= 100 bp): %d\n", nrow(unique_seqs)))

# Save viral sequences regardless
viral_seqs <- DNAStringSet(unique_seqs$viral_sequence)
names(viral_seqs) <- paste0("seq", seq_along(viral_seqs), "_",
                             unique_seqs$viral_orientation, "_",
                             str_sub(unique_seqs$READ, 1, 30))
viral_fasta_file <- paste0(output_prefix, "_viral_sequences.fasta")
writeXStringSet(viral_seqs, viral_fasta_file)
cat(sprintf("  - Saved viral sequences to: %s\n", viral_fasta_file))

compute_kmer_similarity <- function(seqs, k = 8L) {
  # Returns an n x n Jaccard-similarity matrix using k-mers
  n <- length(seqs)
  if (n < 2) return(diag(n))

  kmer_sets <- lapply(seqs, function(s) {
    sv <- as.character(s)
    L <- nchar(sv)
    if (L < k) return(character(0))
    unique(substring(sv, 1:(L - k + 1), k:(L)))
  })

  mat <- diag(n)
  rownames(mat) <- names(seqs)
  colnames(mat) <- names(seqs)

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      si <- kmer_sets[[i]]; sj <- kmer_sets[[j]]
      if (length(si) == 0 || length(sj) == 0) next
      inter <- length(intersect(si, sj))
      union <- length(union(si, sj))
      jac <- if (union > 0) inter / union else 0
      mat[i, j] <- mat[j, i] <- jac
    }
  }
  mat
}

n_seqs <- length(viral_seqs)

if (n_seqs >= 2 && n_seqs <= 200) {
  cat("  - Computing k-mer Jaccard similarity matrix (k=8)...\n")
  sim_mat <- compute_kmer_similarity(viral_seqs, k = 8L)

  sim_file <- paste0(output_prefix, "_similarity_matrix.csv")
  write_csv(as_tibble(sim_mat, rownames = "sequence"), sim_file)
  cat(sprintf("  - Similarity matrix saved: %s\n", sim_file))

  # Diversity: mean pairwise dissimilarity
  off_diag <- sim_mat[lower.tri(sim_mat)]
  cat(sprintf("  - Mean pairwise similarity   : %.3f\n", mean(off_diag)))
  cat(sprintf("  - Mean pairwise dissimilarity: %.3f\n", 1 - mean(off_diag)))

  # Heatmap (base R – no extra deps)
  if (n_seqs <= 50) {
    heatmap_file <- paste0(output_prefix, "_similarity_heatmap.pdf")
    pdf(heatmap_file, width = 12, height = 10)
    heatmap(sim_mat,
            main = "Viral Sequence k-mer Similarity",
            scale = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            margins = c(10, 10))
    dev.off()
    cat(sprintf("  - Heatmap saved: %s\n", heatmap_file))
  }

  # NJ phylogenetic tree
  if (n_seqs >= 3 && n_seqs <= 100) {
    cat("  - Building NJ phylogenetic tree...\n")
    dist_mat <- as.dist(1 - sim_mat)
    tree <- nj(dist_mat)
    tree_file <- paste0(output_prefix, "_tree.nwk")
    write.tree(tree, tree_file)
    tree_pdf <- paste0(output_prefix, "_tree.pdf")
    pdf(tree_pdf, width = 12, height = max(8, n_seqs * 0.4))
    plot(tree, cex = 0.7,
         main = "Phylogenetic Tree – Viral Sequences (NJ, k-mer Jaccard)",
         sub = "Distance = 1 − Jaccard similarity (k=8)")
    add.scale.bar()
    dev.off()
    cat(sprintf("  - Tree saved: %s\n", tree_file))
  }
}

# ============================================================================
# STEP 9: Console summary
# ============================================================================
cat("\n", strrep("-", 70), "\n")
cat("SEQUENCE LENGTH SUMMARY\n")
cat(strrep("-", 70), "\n")

length_stats <- df %>%
  filter(!is.na(viral_seq_length), viral_seq_length > 0) %>%
  summarise(n = n(),
            mean_length = mean(viral_seq_length),
            median_length = median(viral_seq_length),
            min_length = min(viral_seq_length),
            max_length = max(viral_seq_length),
            sd_length = sd(viral_seq_length))

if (length_stats$n > 0) {
  cat(sprintf("N sequences : %d\n",      length_stats$n))
  cat(sprintf("Mean length : %.1f bp\n", length_stats$mean_length))
  cat(sprintf("Median      : %.0f bp\n", length_stats$median_length))
  cat(sprintf("Range       : %.0f - %.0f bp\n",
              length_stats$min_length, length_stats$max_length))
}

cat("\n", strrep("-", 70), "\n")
cat("INTEGRATION SITE STATISTICS\n")
cat(strrep("-", 70), "\n")
cat(sprintf("Total reads            : %d\n", nrow(df)))
cat(sprintf("Unique integration sites: %d\n",
            n_distinct(df$integration_site, na.rm = TRUE)))
if ("gene_name" %in% colnames(df))
  cat(sprintf("Reads in annotated genes: %d\n", sum(!is.na(df$gene_name))))

# ============================================================================
# STEP 10: Save annotated CSV
# ============================================================================
output_csv <- paste0(output_prefix, "_annotated.csv")
write_csv(df, output_csv)
cat(sprintf("\n  Saved annotated CSV: %s\n", output_csv))

# Plain-text report
report_file <- paste0(output_prefix, "_report.txt")
sink(report_file)
cat("VIRAL ORIENTATION AND GENOMIC ANNOTATION REPORT\n")
cat(strrep("=", 70), "\n\n")
cat("Date       :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input CSV  :", input_csv, "\n")
cat("FASTA      :", input_fasta, "\n")
cat("HIV Ref    :", hiv_ref_fasta, "\n")
cat("GTF        :", gtf_file, "\n")
cat("SAM (iter1):", if (!is.null(iter1_sam)) iter1_sam else "not provided", "\n\n")

cat("ORIENTATION METHOD\n", strrep("-", 70), "\n")
cat("Primary  : SAM FLAG bit-16 from iteration-1 alignment\n")
cat("Fallback : k-mer vote (k=21) against 5′/3′ HIV reference regions\n\n")

ori_summary <- df %>% count(viral_orientation) %>% arrange(desc(n))
cat("ORIENTATION COUNTS\n", strrep("-", 70), "\n")
for (i in seq_len(nrow(ori_summary)))
  cat(sprintf("  %-12s: %d\n", ori_summary$viral_orientation[i], ori_summary$n[i]))

if (length_stats$n > 0) {
  cat("\nSEQUENCE LENGTHS\n", strrep("-", 70), "\n")
  cat(sprintf("N=%d  mean=%.1f  median=%.0f  range=[%.0f,%.0f]\n",
              length_stats$n, length_stats$mean_length, length_stats$median_length,
              length_stats$min_length, length_stats$max_length))
}
if (n_seqs >= 2 && n_seqs <= 200 && exists("off_diag")) {
  cat("\nSEQUENCE DIVERSITY (k-mer Jaccard)\n", strrep("-", 70), "\n")
  cat(sprintf("Mean pairwise similarity   : %.3f\n", mean(off_diag)))
  cat(sprintf("Mean pairwise dissimilarity: %.3f\n", 1 - mean(off_diag)))
}
pid_data <- df %>% filter(!is.na(percent_identity))
if (nrow(pid_data) > 0) {
  cat("\nPERCENT IDENTITY TO VIRAL REFERENCE\n", strrep("-", 70), "\n")
  cat(sprintf("N reads with PID  : %d\n",  nrow(pid_data)))
  cat(sprintf("Mean              : %.2f%%\n", mean(pid_data$percent_identity)))
  cat(sprintf("Median            : %.2f%%\n", median(pid_data$percent_identity)))
  cat(sprintf("Range             : %.2f%% - %.2f%%\n",
              min(pid_data$percent_identity), max(pid_data$percent_identity)))
}
sink()
cat(sprintf("  Saved text report : %s\n", report_file))

# ============================================================================
# STEP 11: Write FASTA for external alignment (MAFFT-ready)
#   Header format: >sample|READ|intactness|strand|chromosome|integration_site
#   Sequences have host Ns stripped (pure viral sequence only)
# ============================================================================
cat("\nStep 11: Writing FASTA for MAFFT alignment...\n")

fasta_seqs <- df %>%
  filter(!is.na(viral_sequence),
         nchar(trimws(as.character(viral_sequence))) > 0) %>%
  mutate(
    seq_clean = gsub("N", "", as.character(viral_sequence))) %>%
  filter(nchar(seq_clean) >= 100) %>%
  mutate(
    intact_lbl = ifelse("intactness" %in% colnames(.),
                        gsub(" ", "_", as.character(intactness)),
                        "NA"),
    strand_lbl = replace_na(as.character(viral_orientation), "NA"),
    chr_lbl = replace_na(as.character(chromosome), "NA"),
    site_lbl = replace_na(as.character(integration_site), "NA"),
    fasta_hdr = paste(sample, READ, intact_lbl,
                        strand_lbl, chr_lbl, site_lbl, sep = "|"),
    fasta_hdr = gsub("[[:space:]]", "_", fasta_hdr))

if (nrow(fasta_seqs) > 0) {
  fasta_lines <- unlist(lapply(seq_len(nrow(fasta_seqs)), function(i) {
    c(paste0(">", fasta_seqs$fasta_hdr[i]),
      fasta_seqs$seq_clean[i])
  }))
  fasta_out <- paste0(output_prefix, "_for_mafft.fasta")
  writeLines(fasta_lines, fasta_out)
  cat(sprintf("  - %d sequences written to: %s\n", nrow(fasta_seqs), fasta_out))
  cat(sprintf("  - Suggested command: mafft --auto --thread -1 %s > %s_aligned.fasta\n",
              fasta_out, output_prefix))
} else {
  cat("  - No sequences >= 100 bp available.\n")
}

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

