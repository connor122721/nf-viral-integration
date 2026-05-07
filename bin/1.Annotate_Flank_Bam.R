#!/usr/bin/env Rscript
#
# Viral Integration Orientation and Sequence Analysis 
#
# Usage:
#   Rscript 1.Annotate_Flank_Bam.R \
#       <input.csv> <unmasked.fasta> <hiv_ref.fasta> <gtf> \
#       <output_prefix> [iteration1.sam] [repeatmasker.bed] [regulatory.bed]
#
# Output
# ------
#   <prefix>_annotated.csv              – per-read annotated table (includes clone_id, clone_class)
#   <prefix>_viral_sequences.fasta      – extracted viral sequences (FASTA)
#   <prefix>_clone_viral_sequences.fasta – viral sequences for clone reads only
#   <prefix>_clone_summary.csv          – one row per clone/unique site
#   <prefix>_similarity_matrix.csv      – pairwise % identity matrix (alignment-based)
#   <prefix>_report.txt                 – plain-text summary
#
# Integration site: single base-pair position (junction point of host flank
#   closest to the viral insert).
# Clone calling: reads whose single-bp integration positions fall within
#   +/- 5 bp of each other are grouped as the same clone.
# Naming convention:
#   Clones  : Clone_{sample}_{chr}:{pos}_{N}
#   PCR dups: Dup_{sample}_{chr}:{pos}_{N}
# Chromatin/regulatory annotation: if BED files are supplied, each
#   integration site is annotated with in_repeat (0/1), repeat_name,
#   repeat_class, and in_regulatory_region columns.

library(tidyverse)
library(data.table)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(ape)

# ============================================================================
# 0. Arguments
# ============================================================================
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript simple_annotate_bam_v2.R",
      "<input.csv> <unmasked.fasta> <hiv_ref.fasta> <gtf>",
      "<output_prefix> [iteration1.sam] [repeatmasker.bed] [regulatory.bed]\n")
  quit(status = 1)
}

input_csv <- args[1]
input_fasta <- args[2]
hiv_ref_fasta <- args[3]
gtf_file <- args[4]
output_prefix <- args[5]
iter1_sam <- if (length(args) >= 6 && args[6] != "NA" && args[6] != "") args[6] else NULL
repeatmasker_bed <- if (length(args) >= 7 && args[7] != "NA" && args[7] != "") args[7] else NULL
regulatory_bed <- if (length(args) >= 8 && args[8] != "NA" && args[8] != "") args[8] else NULL

# Test
# input_csv="X:/nf-viral-integration_t2t/test_results/04_final_results/92UG_029.subset.csv"; input_fasta="X:/nf-viral-integration_t2t/test_results/unmasked_sequences/92UG_029.subset.unmasked.fa"; hiv_ref_fasta="X:/nf-viral-integration_t2t/data/hiv_genome_panal/A1.UG..UG031.AB098330.fasta"; gtf_file="X:/nf-viral-integration_t2t/test_results/genome_files/host.gtf"; output_prefix="test1"; iter1_sam="X:/nf-viral-integration_t2t/test_results/02_iterative_masking/92UG_029.subset/92UG_029.subset.dedup.1.sam"

cat("\n", strrep("=", 70), "\n")
cat("VIRAL ORIENTATION AND GENOMIC ANNOTATION ANALYSIS\n")
cat(strrep("=", 70), "\n\n")

# ============================================================================
# STEP 1: Read CSV data
# ============================================================================
cat("Step 1: Reading integration data...\n")
df <- fread(input_csv) %>%
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
      ops <- base::regmatches(cig, base::gregexpr("[0-9]+[MIDNSHPX=]", cig))[[1]]
      widths <- as.integer(sub("[A-Z=]$", "", ops))
      types <- sub("^[0-9]+",  "", ops)
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
    mutate(strand = coalesce(strand_sam, strand),
           orientation = coalesce(orientation_sam, orientation),
           ref_start = coalesce(ref_start_sam, ref_start),
           ref_end = coalesce(ref_end_sam, ref_end),
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
K <- 21
SAMPLE_STEP <- 50
KMER_REGION <- min(1000, floor(ref_len / 3))

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
  cat(sprintf(" %-12s: %d\n", ori_table$orientation[i], ori_table$n[i]))
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
  dplyr::select(-c(read_match, RTF_NUM, HUMAN_GROUP, 
                HIV_DIR_ERR, FLANK_DIR_ERR, HUMAN_MAP_ERR, 
                OVERLAP_ERR, UNMAPPED)))

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
  # Return a single base-pair position: the host coordinate at the junction

  # between the flank and the viral insert.
  # For LEFT_FLANK  the junction is the END   coordinate (rightmost host bp).
  # For RIGHT_FLANK the junction is the START coordinate (leftmost host bp).
  # If both flanks are present, use LEFT_FLANK (upstream junction).
  if (!is.na(left_flank) && left_flank != "") {
    loc <- parse_location(left_flank)
    if (!is.na(loc$chrom))
      return(sprintf("%s:%d", loc$chrom, loc$end))
  }
  if (!is.na(right_flank) && right_flank != "") {
    loc <- parse_location(right_flank)
    if (!is.na(loc$chrom))
      return(sprintf("%s:%d", loc$chrom, loc$start))
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
    integration_pos = map_int(integration_site, ~ {
      if (is.na(.x)) { NA_integer_ } else { as.integer(str_split(.x, ":", n = 2)[[1]][2]) }
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

# Use READ as a stable key — it is unique per row and survives filtering
has_site <- !is.na(df$integration_pos) & !is.na(df$chromosome)
df$gene_name <- NA_character_
df$gene_id   <- NA_character_

if (any(has_site)) {
  query_gr <- GRanges(
    seqnames = df$chromosome[has_site],
    ranges = IRanges(start = df$integration_pos[has_site],
                     end   = df$integration_pos[has_site]))

  overlaps <- findOverlaps(query_gr, genes)

  if (length(overlaps) > 0) {
    # Map query indices back to df row positions
    df_rows_with_site <- which(has_site)

    ann <- tibble(
      df_row  = df_rows_with_site[queryHits(overlaps)],
      g_name  = genes$gene_name[subjectHits(overlaps)],
      g_id    = genes$gene_id[subjectHits(overlaps)])

    ann_agg <- ann %>%
      group_by(df_row) %>%
      summarize(gene_name = paste(unique(na.omit(g_name)), collapse = ";"),
                gene_id   = paste(unique(na.omit(g_id)),   collapse = ";"),
                .groups = "drop") %>%
      mutate(gene_name = if_else(gene_name == "", NA_character_, gene_name),
             gene_id   = if_else(gene_id   == "", NA_character_, gene_id))

    df$gene_name[ann_agg$df_row] <- ann_agg$gene_name
    df$gene_id[ann_agg$df_row]   <- ann_agg$gene_id

    cat(sprintf("  - Annotated %d reads with gene information\n",
                sum(!is.na(df$gene_name))))
  } else {
    cat("  - No overlaps found between integration sites and genes\n")
  }
} else {
  cat("  - No valid integration sites to annotate\n")
}

# ============================================================================
# STEP 7b: RepeatMasker and regulatory region annotation
# ============================================================================
cat("\nStep 7b: Annotating repeats and regulatory regions...\n")

# --- RepeatMasker annotation ------------------------------------------------
# The rmsk BED has col4 = "name#class/family" (e.g. "L1MC3#LINE/L1").
# We annotate: in_repeat (0/1), repeat_name, repeat_class

df$in_repeat     <- 0L
df$repeat_name   <- NA_character_
df$repeat_class  <- NA_character_

if (!is.null(repeatmasker_bed) && file.exists(repeatmasker_bed)) {
  cat(sprintf("  - Loading RepeatMasker BED: %s\n", basename(repeatmasker_bed)))
  rmsk <- fread(repeatmasker_bed, header = FALSE, select = 1:4,
                col.names = c("chrom", "start", "end", "name_class"))
  # Parse name#class/family
  rmsk[, c("rpt_name", "rpt_class") := tstrsplit(name_class, "#", fixed = TRUE, keep = 1:2)]
  rmsk_gr <- GRanges(seqnames = rmsk$chrom,
                     ranges = IRanges(start = rmsk$start + 1L, end = rmsk$end),
                     rpt_name  = rmsk$rpt_name,
                     rpt_class = rmsk$rpt_class)
  cat(sprintf("    %d repeat regions loaded\n", length(rmsk_gr)))

  has_site <- !is.na(df$integration_pos) & !is.na(df$chromosome)
  if (any(has_site)) {
    query_gr <- GRanges(seqnames = df$chromosome[has_site],
                        ranges = IRanges(start = df$integration_pos[has_site],
                                         end   = df$integration_pos[has_site]))
    hits <- findOverlaps(query_gr, rmsk_gr)

    if (length(hits) > 0) {
      df_rows_with_site <- which(has_site)
      hit_df <- tibble(
        df_row    = df_rows_with_site[queryHits(hits)],
        rpt_name  = rmsk_gr$rpt_name[subjectHits(hits)],
        rpt_class = rmsk_gr$rpt_class[subjectHits(hits)]) %>%
        group_by(df_row) %>%
        summarize(repeat_name  = paste(unique(na.omit(rpt_name)),  collapse = ";"),
                  repeat_class = paste(unique(na.omit(rpt_class)), collapse = ";"),
                  .groups = "drop")

      df$in_repeat[hit_df$df_row]     <- 1L
      df$repeat_name[hit_df$df_row]   <- hit_df$repeat_name
      df$repeat_class[hit_df$df_row]  <- hit_df$repeat_class
    }
    cat(sprintf("    %d / %d sites overlap repeats\n",
                sum(df$in_repeat), sum(has_site)))
  }
} else {
  cat("  - RepeatMasker BED not provided, skipping\n")
}

# --- Regulatory region annotation -------------------------------------------
# Binary overlap with regulatory BED (PLS/pELS/dELS from ENCODE cCREs)
df$in_regulatory_region <- 0L

if (!is.null(regulatory_bed) && file.exists(regulatory_bed)) {
  cat(sprintf("  - Loading regulatory BED: %s\n", basename(regulatory_bed)))
  reg <- fread(regulatory_bed, header = FALSE, select = 1:3,
               col.names = c("chrom", "start", "end"))
  reg_gr <- GRanges(seqnames = reg$chrom,
                    ranges = IRanges(start = reg$start + 1L, end = reg$end))
  cat(sprintf("    %d regulatory regions loaded\n", length(reg_gr)))

  has_site <- !is.na(df$integration_pos) & !is.na(df$chromosome)
  if (any(has_site)) {
    query_gr <- GRanges(seqnames = df$chromosome[has_site],
                        ranges = IRanges(start = df$integration_pos[has_site],
                                         end   = df$integration_pos[has_site]))
    hits <- overlapsAny(query_gr, reg_gr)
    df$in_regulatory_region[which(has_site)] <- as.integer(hits)
    cat(sprintf("    %d / %d sites overlap regulatory regions\n",
                sum(hits), sum(has_site)))
  }
} else {
  cat("  - Regulatory BED not provided, skipping\n")
}

# ============================================================================
# STEP 8: Clone / PCR-replicate classification
#
# Integration site: single base-pair position (junction point).
# Clone grouping: reads whose integration positions are within +/- 5 bp
#   on the same chromosome are grouped together.
# Naming convention:
#   Clones  : Clone_{sample}_{chr}:{canonical_pos}_{N}
#   PCR dups: Dup_{sample}_{chr}:{canonical_pos}_{N}
# ============================================================================
cat("\nStep 8: Clone / PCR-replicate classification...\n")

sample_name <- sub("\\.csv$", "", basename(input_csv))
CLONE_WINDOW <- 5L   # +/- 5 bp

# --- 8a. Group integration sites within +/- 5 bp --------------------------
#     Sort by chromosome and position, then greedily merge positions that
#     fall within the window.  The canonical position for each group is the
#     median position (rounded).

if (all(c("integration_pos", "chromosome") %in% colnames(df))) {

  sites <- df %>%
    filter(!is.na(integration_pos), !is.na(chromosome)) %>%
    distinct(chromosome, integration_pos) %>%
    arrange(chromosome, integration_pos)

  # Greedy clustering: walk sorted positions; start a new group when the

  # gap to the previous position exceeds CLONE_WINDOW.
  if (nrow(sites) > 0) {
    sites$group <- 1L
    g <- 1L
    for (i in seq_len(nrow(sites))) {
      if (i == 1L) next
      if (sites$chromosome[i] != sites$chromosome[i - 1L] ||
          (sites$integration_pos[i] - sites$integration_pos[i - 1L]) > CLONE_WINDOW) {
        g <- g + 1L
      }
      sites$group[i] <- g
    }

    # Canonical position = median of group positions (integer)
    site_groups <- sites %>%
      group_by(group) %>%
      mutate(canonical_pos = as.integer(round(median(integration_pos)))) %>%
      ungroup() %>%
      dplyr::select(chromosome, integration_pos, canonical_pos, group)

    df <- df %>%
      left_join(site_groups, by = c("chromosome", "integration_pos"))
  } else {
    df <- df %>% mutate(canonical_pos = NA_integer_, group = NA_integer_)
  }
} else {
  df <- df %>% mutate(canonical_pos = NA_integer_, group = NA_integer_)
}

# --- 8b. Assign clone_id / dup_id with structured naming -------------------
#     site_key is now chr:canonical_pos.  Within each group:
#       - If 1 read  → unique clone  → Clone_{sample}_{site}_{1}
#       - If >1 reads → first read is the clone, rest are PCR dups
#         Clone: Clone_{sample}_{site}_{1}
#         Dups : Dup_{sample}_{site}_{2}, Dup_{sample}_{site}_{3}, ...

df <- df %>%
  mutate(
    site_key = case_when(
      !is.na(canonical_pos) & !is.na(chromosome) ~
        paste0(as.character(chromosome), ":", canonical_pos),
      TRUE ~ NA_character_))

# Number reads within each site_key group (order by viral_seq_length desc
# so the longest read gets _1 = clone representative)
df <- df %>%
  group_by(site_key) %>%
  mutate(
    reads_at_site = sum(!is.na(site_key)),
    rank_in_site = rank(-replace_na(viral_seq_length, 0L),
                        ties.method = "first")) %>%
  ungroup()

df <- df %>%
  mutate(
    clone_id = case_when(
      is.na(site_key) ~ NA_character_,
      rank_in_site == 1L ~
        paste0("Clone_", sample_name, "_", site_key, "_", rank_in_site),
      TRUE ~
        paste0("Dup_", sample_name, "_", site_key, "_", rank_in_site)),
    clone_class = case_when(
      is.na(site_key) ~ "no host flanks",
      reads_at_site == 1L ~ "unique integration",
      rank_in_site == 1L ~ "clone representative",
      TRUE ~ "PCR duplicate"),
    is_pcr_replicate = (!is.na(reads_at_site) & reads_at_site > 1L &
                          rank_in_site > 1L))

n_unique_clones <- n_distinct(df$site_key, na.rm = TRUE)
n_pcr_reads <- sum(df$is_pcr_replicate, na.rm = TRUE)
cat(sprintf("  - Unique integration sites (clones): %d\n", n_unique_clones))
cat(sprintf("  - PCR duplicate reads (same site):   %d\n", n_pcr_reads))

# --- 8c. Alignment-based sequence similarity within site groups ----------
#     For each site with >1 read, align all reads to the first read in the
#     group (representative).  This confirms PCR-dup status and flags
#     divergent reads at the same site.

cat("  - Computing within-site alignment similarities...\n")

# Helper: pairwise % identity via Biostrings local alignment
safe_pairwise_pid <- function(s1, s2, max_len = 5000L) {
  # Truncate to max_len for speed (SMRTCap reads can be long)
  s1_t <- subseq(DNAString(s1), 1L, min(nchar(s1), max_len))
  s2_t <- subseq(DNAString(s2), 1L, min(nchar(s2), max_len))
  tryCatch({
    aln <- pairwiseAlignment(s1_t, s2_t,
                              type = "local",
                              substitutionMatrix = nucleotideSubstitutionMatrix(
                                match = 1, mismatch = -1, baseOnly = TRUE),
                              gapOpening = 5, gapExtension = 2)
    # pid = (aligned length - mismatches - indels) / aligned length
    aln_len <- nchar(aln)
    if (aln_len == 0) return(0)
    round(pwalign::pid(aln), 2)
  }, error = function(e) NA_real_)
}

# Within-site comparison
within_site_results <- tibble()
multi_read_sites <- df %>%
  filter(!is.na(site_key), reads_at_site > 1) %>%
  distinct(site_key) %>%
  pull(site_key)

if (length(multi_read_sites) > 0) {
  ws_list <- list()
  for (sk in multi_read_sites) {
    grp <- df %>%
      filter(site_key == sk, !is.na(viral_sequence), nchar(viral_sequence) >= 100)
    if (nrow(grp) < 2) next

    ref_seq <- grp$viral_sequence[1]
    ref_read <- grp$READ[1]

    for (i in 2:nrow(grp)) {
      p <- safe_pairwise_pid(ref_seq, grp$viral_sequence[i])
      ws_list[[length(ws_list) + 1]] <- tibble(
        site_key = sk,
        read_A = ref_read,
        read_B = grp$READ[i],
        pct_identity = p,
        same_site = TRUE,
        classification = case_when(
          is.na(p) ~ "undetermined",
          p >= 95 ~ "PCR duplicate",
          p >= 90 ~ "same-site divergent",
          TRUE ~ "same-site highly divergent"))
    }
  }
  within_site_results <- bind_rows(ws_list)
  cat(sprintf("  - Within-site comparisons: %d pairs\n", nrow(within_site_results)))

  # Update clone_class for divergent reads at same site
  divergent_reads <- within_site_results %>%
    filter(!is.na(pct_identity), pct_identity < 97) %>%
    pull(read_B)

  if (length(divergent_reads) > 0) {
    df <- df %>%
      mutate(clone_class = if_else(
        READ %in% divergent_reads & clone_class == "PCR duplicate",
        "same-site divergent",
        clone_class))
    cat(sprintf("  - Same-site divergent reads reclassified: %d\n",
                length(divergent_reads)))
  }
}

# --- 8d. Cross-site clonal expansion detection ---------------------------
#     Among unique-site representative reads, compute pairwise alignment
#     similarity.  Pairs with ≥95% identity at *different* sites may
#     indicate clonal expansion from the same founder virus.

cat("  - Checking for cross-site clonal expansion...\n")

# Take one representative read per clone (longest viral sequence)
clone_reps <- df %>%
  filter(!is.na(clone_id), !is.na(viral_sequence), nchar(viral_sequence) >= 500) %>%
  group_by(clone_id) %>%
  slice_max(viral_seq_length, n = 1, with_ties = FALSE) %>%
  ungroup()

n_reps <- nrow(clone_reps)
cross_site_results <- tibble()

# Only compute if tractable (≤100 clones → ≤4950 pairs)
CROSS_SITE_MAX_CLONES <- 100L
CROSS_SITE_IDENTITY_THRESHOLD <- 95.0

if (n_reps >= 2 && n_reps <= CROSS_SITE_MAX_CLONES) {
  cat(sprintf("  - Comparing %d clone representatives pairwise...\n", n_reps))
  cs_list <- list()
  for (i in seq_len(n_reps - 1)) {
    for (j in (i + 1):n_reps) {
      p <- safe_pairwise_pid(clone_reps$viral_sequence[i],
                              clone_reps$viral_sequence[j])
      if (!is.na(p) && p >= CROSS_SITE_IDENTITY_THRESHOLD) {
        cs_list[[length(cs_list) + 1]] <- tibble(
          clone_A = clone_reps$clone_id[i],
          clone_B = clone_reps$clone_id[j],
          site_A = clone_reps$site_key[i],
          site_B = clone_reps$site_key[j],
          pct_identity = p,
          note = "possible clonal expansion (same virus, different site)")
      }
    }
  }
  cross_site_results <- bind_rows(cs_list)
  if (nrow(cross_site_results) > 0) {
    cat(sprintf("  - Cross-site clonal candidates: %d pairs (>= %.0f%% identity)\n",
                nrow(cross_site_results), CROSS_SITE_IDENTITY_THRESHOLD))
  } else {
    cat("  - No cross-site clonal expansion detected.\n")
  }
} else if (n_reps > CROSS_SITE_MAX_CLONES) {
  cat(sprintf("  - Skipping cross-site comparison (%d clones > %d limit)\n",
              n_reps, CROSS_SITE_MAX_CLONES))
}

# ============================================================================
# STEP 8e: Within-group PID + reference-coordinate similarity matrix
#
# Two-part approach:
#   A) Within-group PID:  For each clone/PCR group (small N), compute all
#      pairwise alignment PIDs and store the group mean.  This is the gold
#      standard for confirming PCR duplicates vs divergent reads.
#   B) Full similarity matrix (≤200 seqs):  Use reference-coordinate overlap
#      for speed.  Two sequences that map to overlapping regions on the HIV
#      reference are compared by k-mer identity over their shared window.
#      Sequences with no reference-coordinate overlap get NA (not 0 — they
#      simply can't be compared meaningfully).
# ============================================================================
cat("\nStep 8e: Within-group PID and similarity matrix...\n")

# --- A) Within-group (clone/PCR) pairwise alignment PID ------------------
#     For each site_key with ≥2 reads, compute all pairwise PIDs and store
#     the mean.  This uses safe_pairwise_pid (local alignment, accurate but
#     slow — fine here because groups are small).

cat("  - Computing within-group pairwise PIDs...\n")

within_group_pids <- tibble()
all_site_keys <- df %>%
  filter(!is.na(site_key), !is.na(viral_sequence), nchar(viral_sequence) >= 100) %>%
  count(site_key) %>%
  filter(n >= 2) %>%
  pull(site_key)

if (length(all_site_keys) > 0) {
  wg_list <- list()
  for (sk in all_site_keys) {
    grp <- df %>%
      filter(site_key == sk, !is.na(viral_sequence), nchar(viral_sequence) >= 100)
    n_grp <- nrow(grp)
    if (n_grp < 2) next

    pair_pids <- numeric(0)
    for (i in seq_len(n_grp - 1)) {
      for (j in (i + 1):min(n_grp, i + 10)) {  # cap at 10 pairs per read for speed
        p <- safe_pairwise_pid(grp$viral_sequence[i], grp$viral_sequence[j])
        if (!is.na(p)) pair_pids <- c(pair_pids, p)
      }
    }

    if (length(pair_pids) > 0) {
      wg_list[[length(wg_list) + 1]] <- tibble(
        site_key = sk,
        clone_id = grp$clone_id[1],
        n_reads_in_group = n_grp,
        n_pairs_compared = length(pair_pids),
        mean_within_pid = round(mean(pair_pids), 2),
        min_within_pid = round(min(pair_pids), 2),
        max_within_pid = round(max(pair_pids), 2))
    }
  }
  within_group_pids <- bind_rows(wg_list)

  if (nrow(within_group_pids) > 0) {
    cat(sprintf("  - Within-group PID computed for %d clone groups\n",
                nrow(within_group_pids)))
    cat(sprintf("  - Overall mean within-group PID: %.1f%%\n",
                mean(within_group_pids$mean_within_pid)))

    # Merge mean_within_pid back into df via site_key
    df <- df %>%
      left_join(within_group_pids %>%
                  dplyr::select(site_key, mean_within_pid),
                by = "site_key")
  } else {
    df <- df %>% mutate(mean_within_pid = NA_real_)
  }
} else {
  cat("  - No multi-read sites to compare.\n")
  df <- df %>% mutate(mean_within_pid = NA_real_)
}

# --- B) Full similarity matrix -------------------------------------------
#     Reference-coordinate-aware k-mer identity.
#     Strategy: for each pair of sequences, check if their viral reference
#     coordinates overlap.  If yes, extract the overlapping sub-sequences
#     and compute k-mer containment (proportion of k-mers shared).
#     If no overlap (different viral regions), mark as NA.
#     This is O(n²) in comparisons but each comparison is fast (no alignment).

unique_seqs <- df %>%
  filter(!is.na(viral_sequence), viral_seq_length >= 100) %>%
  distinct(viral_sequence, .keep_all = TRUE)

cat(sprintf("  - Unique viral sequences (>= 100 bp): %d\n", nrow(unique_seqs)))

# Save viral sequences
viral_seqs <- DNAStringSet(unique_seqs$viral_sequence)
names(viral_seqs) <- paste0(
  ifelse(!is.na(unique_seqs$clone_id), unique_seqs$clone_id, "no_host_flanks"),
  "_", unique_seqs$viral_orientation,
  "_", str_sub(unique_seqs$READ, 1, 30))
viral_fasta_file <- paste0(output_prefix, "_viral_sequences.fasta")
writeXStringSet(viral_seqs, viral_fasta_file)
cat(sprintf("  - Saved viral sequences to: %s\n", viral_fasta_file))

# Clone-specific FASTA: only reads classified as clone representatives
clone_seqs_df <- df %>%
  filter(grepl("^Clone_", clone_id),
         !is.na(viral_sequence),
         nchar(trimws(as.character(viral_sequence))) > 0) %>%
  mutate(seq_clean = gsub("N", "", as.character(viral_sequence))) %>%
  filter(nchar(seq_clean) >= 100)

if (nrow(clone_seqs_df) > 0) {
  clone_dss <- DNAStringSet(clone_seqs_df$seq_clean)
  names(clone_dss) <- paste0(
    clone_seqs_df$clone_id, "|",
    replace_na(as.character(clone_seqs_df$viral_orientation), "NA"), "|",
    replace_na(as.character(clone_seqs_df$integration_site), "NA"))
  clone_fasta_file <- paste0(output_prefix, "_clone_viral_sequences.fasta")
  writeXStringSet(clone_dss, clone_fasta_file)
  cat(sprintf("  - Clone-specific viral FASTA: %d sequences written to %s\n",
              length(clone_dss), clone_fasta_file))
} else {
  cat("  - No clone representatives with viral sequence >= 100 bp\n")
}

n_seqs <- length(viral_seqs)
SIM_MATRIX_MAX <- 1000

# k-mer containment function (fast, no alignment)
kmer_identity <- function(s1, s2, k = 12L) {
  if (nchar(s1) < k || nchar(s2) < k) return(NA_real_)
  kmers1 <- unique(substring(s1, 1:(nchar(s1) - k + 1), k:(nchar(s1))))
  kmers2 <- unique(substring(s2, 1:(nchar(s2) - k + 1), k:(nchar(s2))))
  shared <- length(intersect(kmers1, kmers2))
  # Containment: fraction of the smaller set found in the larger
  denom <- min(length(kmers1), length(kmers2))
  if (denom == 0) return(0)
  round(100 * shared / denom, 2)
}

if (n_seqs >= 2 && n_seqs <= SIM_MATRIX_MAX) {
  cat(sprintf("  - Computing %d x %d similarity matrix (k-mer containment, k=12)...\n",
              n_seqs, n_seqs))

  # Extract ref coordinates for overlap checking
  ref_starts <- unique_seqs$ref_start
  ref_ends   <- unique_seqs$ref_end
  seq_chars  <- as.character(viral_seqs)

  sim_mat <- diag(100, n_seqs)
  rownames(sim_mat) <- names(viral_seqs)
  colnames(sim_mat) <- names(viral_seqs)

  for (i in seq_len(n_seqs - 1)) {
    for (j in (i + 1):n_seqs) {
      # Check reference-coordinate overlap
      has_coords <- !is.na(ref_starts[i]) && !is.na(ref_ends[i]) &&
                    !is.na(ref_starts[j]) && !is.na(ref_ends[j])

      if (has_coords) {
        overlap_start <- max(ref_starts[i], ref_starts[j])
        overlap_end   <- min(ref_ends[i], ref_ends[j])
        overlap_len   <- overlap_end - overlap_start

        if (overlap_len < 100) {
          # No meaningful overlap — different viral regions
          sim_mat[i, j] <- sim_mat[j, i] <- NA
          next
        }
      }

      # Compute k-mer containment identity
      p <- kmer_identity(seq_chars[i], seq_chars[j], k = 12L)
      sim_mat[i, j] <- p
      sim_mat[j, i] <- p
    }
    if (n_seqs > 50 && i %% 20 == 0)
      cat(sprintf("    ... row %d / %d\n", i, n_seqs))
  }

  sim_file <- paste0(output_prefix, "_similarity_matrix.csv")
  write_csv(as_tibble(sim_mat, rownames = "sequence"), sim_file)
  cat(sprintf("  - Similarity matrix saved: %s\n", sim_file))

  # Diversity stats (excluding NA pairs — different-region comparisons)
  off_diag <- sim_mat[lower.tri(sim_mat)]
  off_diag_valid <- off_diag[!is.na(off_diag)]
  if (length(off_diag_valid) > 0) {
    cat(sprintf("  - Comparable pairs : %d / %d (%.0f%%)\n",
                length(off_diag_valid), length(off_diag),
                100 * length(off_diag_valid) / length(off_diag)))
    cat(sprintf("  - Mean pairwise identity : %.1f%%\n", mean(off_diag_valid)))
    cat(sprintf("  - Median pairwise identity : %.1f%%\n", median(off_diag_valid)))
  }

  # Heatmap with clone annotations
  if (n_seqs <= 50) {
    heatmap_file <- paste0(output_prefix, "_similarity_heatmap.pdf")

    # Replace NA with 0 for heatmap visualisation (mark as grey)
    sim_mat_plot <- sim_mat
    sim_mat_plot[is.na(sim_mat_plot)] <- 0

    pdf(heatmap_file, width = 12, height = 10)
    heatmap(sim_mat_plot,
            symm = TRUE,
            scale = "none",
            col = colorRampPalette(c("#e9ecef", "#2166ac", "#f7f7f7", "#b2182b"))(100),
            main = "Viral Sequence Similarity (k-mer Containment, k=12)\nGrey = no reference overlap (different viral regions)",
            margins = c(10, 10))
    dev.off()
    cat(sprintf("  - Heatmap saved: %s\n", heatmap_file))
  }

  # NJ tree (use only comparable pairs; set NA distances to max observed + 10)
  if (n_seqs >= 3 && n_seqs <= 100) {
    cat("  - Building NJ phylogenetic tree...\n")
    dist_vals <- 100 - sim_mat
    max_dist <- max(dist_vals, na.rm = TRUE) + 10
    dist_vals[is.na(dist_vals)] <- max_dist
    diag(dist_vals) <- 0
    dist_mat <- as.dist(dist_vals)
    tree <- tryCatch({
      nj(dist_mat)
    }, error = function(e) {
      cat(sprintf("  - NJ tree failed: %s\n", conditionMessage(e)))
      NULL
    })
    if (!is.null(tree)) {
      tree_file <- paste0(output_prefix, "_tree.nwk")
      write.tree(tree, tree_file)
      tree_pdf <- paste0(output_prefix, "_tree.pdf")
      pdf(tree_pdf, width = 12, height = max(8, n_seqs * 0.4))
      plot(tree, cex = 0.7,
           main = "Phylogenetic Tree – Viral Sequences (NJ, k-mer Containment Distance)",
           sub = "Distance = 100 - % k-mer containment (k=12); NA pairs set to max distance")
      add.scale.bar()
      dev.off()
      cat(sprintf("  - Tree saved: %s\n", tree_file))
    }
  }
} else if (n_seqs > SIM_MATRIX_MAX) {
  cat(sprintf("  - Skipping full matrix (%d seqs > %d limit). ",
              n_seqs, SIM_MATRIX_MAX))
  cat("Clone/PCR classification still applied via site-based grouping.\n")
}

# ============================================================================
# STEP 9: Clone summary table
# ============================================================================
cat("\nStep 9: Building clone summary...\n")

clone_summary <- df %>%
  filter(!is.na(clone_id)) %>%
  group_by(site_key) %>%
  summarise(
    clone_id = dplyr::first(clone_id[grepl("^Clone_", clone_id)]),
    n_reads = n(),
    n_pcr_duplicates = sum(grepl("^Dup_", clone_id)),
    representative_read = READ[which.max(viral_seq_length)],
    chromosome = dplyr::first(na.omit(chromosome)),
    integration_site = dplyr::first(na.omit(integration_site)),
    canonical_pos = dplyr::first(na.omit(canonical_pos)),
    gene_name = dplyr::first(na.omit(gene_name)),
    viral_orientation = dplyr::first(na.omit(viral_orientation)),
    mean_insert_len = round(mean(INSERT_LEN, na.rm = TRUE)),
    mean_pct_identity_to_ref = round(mean(percent_identity, na.rm = TRUE), 1),
    mean_within_group_pid = round(dplyr::first(na.omit(mean_within_pid)), 1),
    in_repeat = max(in_repeat, na.rm = TRUE),
    repeat_name = dplyr::first(na.omit(repeat_name)),
    repeat_class = dplyr::first(na.omit(repeat_class)),
    in_regulatory_region = max(in_regulatory_region, na.rm = TRUE),
    .groups = "drop") %>%
  mutate(in_repeat = replace(in_repeat, is.infinite(in_repeat), 0L),
         in_regulatory_region = replace(in_regulatory_region, is.infinite(in_regulatory_region), 0L)) %>%
  arrange(desc(n_reads))

clone_csv <- paste0(output_prefix, "_clone_summary.csv")
write_csv(clone_summary, clone_csv)
cat(sprintf("  - Clone summary: %d clones written to %s\n",
            nrow(clone_summary), clone_csv))

# Write cross-site results if any
if (nrow(cross_site_results) > 0) {
  cross_csv <- paste0(output_prefix, "_cross_site_clonal.csv")
  write_csv(cross_site_results, cross_csv)
  cat(sprintf("  - Cross-site clonal candidates: %s\n", cross_csv))
}

# Write within-site comparison results
if (nrow(within_site_results) > 0) {
  ws_csv <- paste0(output_prefix, "_within_site_comparisons.csv")
  write_csv(within_site_results, ws_csv)
  cat(sprintf("  - Within-site comparisons: %s\n", ws_csv))
}

# ============================================================================
# STEP 10: Console summary
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
cat(sprintf("Total reads              : %d\n", nrow(df)))
cat(sprintf("Unique integration sites : %d\n", n_unique_clones))
cat(sprintf("PCR replicate reads      : %d\n", n_pcr_reads))
if ("gene_name" %in% colnames(df))
  cat(sprintf("Reads in annotated genes : %d\n", sum(!is.na(df$gene_name))))
if (nrow(cross_site_results) > 0)
  cat(sprintf("Cross-site clonal pairs  : %d (>= %.0f%% seq identity)\n",
              nrow(cross_site_results), CROSS_SITE_IDENTITY_THRESHOLD))

# ============================================================================
# STEP 11: Save annotated CSV
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
cat("SAM (iter1):", if (!is.null(iter1_sam)) iter1_sam else "not provided", "\n")
cat("RepeatMask :", if (!is.null(repeatmasker_bed)) repeatmasker_bed else "not provided", "\n")
cat("Regulatory :", if (!is.null(regulatory_bed)) regulatory_bed else "not provided", "\n\n")

cat("ORIENTATION METHOD\n", strrep("-", 70), "\n")
cat("Primary  : SAM FLAG bit-16 from iteration-1 alignment\n")
cat("Fallback : k-mer vote (k=21) against 5'/3' HIV reference regions\n\n")

ori_summary <- df %>% count(viral_orientation) %>% arrange(desc(n))
cat("ORIENTATION COUNTS\n", strrep("-", 70), "\n")
for (i in seq_len(nrow(ori_summary)))
  cat(sprintf("  %-12s: %d\n", ori_summary$viral_orientation[i], ori_summary$n[i]))

cat("\nCLONE / PCR-REPLICATE CLASSIFICATION\n", strrep("-", 70), "\n")
cat("Integration site: single base-pair (junction point)\n")
cat("Clone grouping  : +/- 5 bp window\n")
cat("Naming          : Clone_{sample}_{chr}:{pos}_{N} / Dup_{sample}_{chr}:{pos}_{N}\n")
cat(sprintf("Unique integration sites (clones): %d\n", n_unique_clones))
cat(sprintf("PCR replicate reads (same site)  : %d\n", n_pcr_reads))
clone_class_tbl <- df %>% count(clone_class) %>% arrange(desc(n))
for (i in seq_len(nrow(clone_class_tbl)))
  cat(sprintf("  %-35s: %d\n", clone_class_tbl$clone_class[i], clone_class_tbl$n[i]))
if (nrow(cross_site_results) > 0) {
  cat(sprintf("\nCross-site clonal expansion candidates: %d pairs\n",
              nrow(cross_site_results)))
  for (i in seq_len(min(nrow(cross_site_results), 10))) {
    cat(sprintf("  %s <-> %s  (%.1f%% identity)\n",
                cross_site_results$clone_A[i],
                cross_site_results$clone_B[i],
                cross_site_results$pct_identity[i]))
  }
}

if (length_stats$n > 0) {
  cat("\nSEQUENCE LENGTHS\n", strrep("-", 70), "\n")
  cat(sprintf("N=%d  mean=%.1f  median=%.0f  range=[%.0f,%.0f]\n",
              length_stats$n, length_stats$mean_length, length_stats$median_length,
              length_stats$min_length, length_stats$max_length))
}
if (n_seqs >= 2 && n_seqs <= SIM_MATRIX_MAX && exists("off_diag_valid") && length(off_diag_valid) > 0) {
  cat("\nSEQUENCE DIVERSITY (k-mer containment, k=12)\n", strrep("-", 70), "\n")
  cat(sprintf("Comparable pairs        : %d\n", length(off_diag_valid)))
  cat(sprintf("Mean pairwise identity  : %.1f%%\n", mean(off_diag_valid)))
  cat(sprintf("Median pairwise identity: %.1f%%\n", median(off_diag_valid)))
}
if (nrow(within_group_pids) > 0) {
  cat("\nWITHIN-GROUP PID (clone/PCR groups)\n", strrep("-", 70), "\n")
  cat(sprintf("Groups with >= 2 reads  : %d\n", nrow(within_group_pids)))
  cat(sprintf("Mean within-group PID   : %.1f%%\n", mean(within_group_pids$mean_within_pid)))
  cat(sprintf("Min within-group PID    : %.1f%%\n", min(within_group_pids$min_within_pid)))
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
# STEP 12: Write FASTA for external alignment (MAFFT-ready)
#   Header: >clone_id|sample|READ|strand|chromosome|integration_site
# ============================================================================
cat("\nStep 12: Writing FASTA for MAFFT alignment...\n")

fasta_seqs <- df %>%
  filter(!is.na(viral_sequence),
         nchar(trimws(as.character(viral_sequence))) > 0) %>%
  mutate(seq_clean = gsub("N", "", as.character(viral_sequence))) %>%
  filter(nchar(seq_clean) >= 100) %>%
  mutate(
    clone_lbl = replace_na(clone_id, "no_host_flank"),
    strand_lbl = replace_na(as.character(viral_orientation), "NA"),
    chr_lbl = replace_na(as.character(chromosome), "NA"),
    site_lbl = replace_na(as.character(integration_site), "NA"),
    fasta_hdr = paste(clone_lbl, sample, READ,
                      strand_lbl, chr_lbl, site_lbl, sep = "|"),
    fasta_hdr = gsub("[[:space:]]", "_", fasta_hdr))

if (nrow(fasta_seqs) > 0) {
  fasta_lines <- unlist(lapply(seq_len(nrow(fasta_seqs)), function(i) {
    c(paste0(">", fasta_seqs$fasta_hdr[i]), fasta_seqs$seq_clean[i])
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
