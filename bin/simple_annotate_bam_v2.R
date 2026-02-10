#!/usr/bin/env Rscript
# ssh bioinfo003
# module load miniconda3; source activate smrtENV; R;
#
# Viral Integration Orientation and Sequence Analysis
# Determines viral orientation by alignment to reference genome
# Includes genomic feature annotation from GTF/GFF
#
# Usage: Rscript viral_orientation_gff_analysis.R <input.csv> <unmasked.fasta> <hiv_ref.fasta> <gtf> <output_prefix>

# Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
  library(Biostrings)
  library(GenomicRanges)
  library(rtracklayer)
  library(ape)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  cat("Usage: Rscript viral_orientation_gff_analysis.R <input.csv> <unmasked.fasta> <hiv_ref.fasta> <gtf> <output_prefix>\n\n")
  cat("Arguments:\n")
  cat("  input.csv       - Integration site CSV file\n")
  cat("  unmasked.fasta  - FASTA with viral sequences (N = host, ATCG = viral)\n")
  cat("  hiv_ref.fasta   - HIV reference genome for orientation determination\n")
  cat("  gtf             - GTF/GFF3 annotation file for genomic features\n")
  cat("  output_prefix   - Prefix for output files\n\n")
  cat("Example:\n")
  cat("  Rscript viral_orientation_gff_analysis.R \\\n")
  cat("    data.csv sequences.fa hiv_ref.fa chm13v2.gtf output\n")
  quit(status = 1)
}

# Input arguments
input_csv <- args[1]
input_fasta <- args[2]
hiv_ref_fasta <- args[3]
gtf_file <- args[4]
output_prefix <- args[5]

# For testing
# input_csv="~/nf-viral-integration_t2t/output/04_final_results/m84248_250919_192326_s1.hifi_reads.univ_v3_bc1002.csv"; input_fasta="~/nf-viral-integration_t2t/output/04_final_results/m84248_250919_192326_s1.hifi_reads.univ_v3_bc1002.unmasked.fa"; hiv_ref_fasta="/home/c0murr09/viral_genome/hiv_genome_panal/Ref.B.FR.83.HXB2_LAI_IIIB_BRU.K03455.fasta"; gtf_file="~/human_genome/t2t/chm13v2.0_RefSeq_Liftoff_v5.2.gtf"; output_prefix="test2"

cat("\n")
cat(strrep("=", 70), "\n")
cat("VIRAL ORIENTATION AND GENOMIC ANNOTATION ANALYSIS\n")
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

# Extract viral sequences (remove N's which represent host sequences)
extract_viral <- function(seq_string) {
  str_replace_all(seq_string, "N", "")
}

# Create sequence lookup table with tidyverse
seq_lookup <- tibble(
  read_id = names(fasta),
  full_sequence = as.character(fasta)) %>%
  mutate(viral_sequence = map_chr(full_sequence, extract_viral),
         viral_seq_length = nchar(viral_sequence),
         read_id_base = str_remove(read_id, "/0$"))

# ============================================================================
# STEP 3: Read HIV reference genome
# ============================================================================
cat("\nStep 3: Reading HIV reference genome...\n")
hiv_ref <- readDNAStringSet(hiv_ref_fasta)
cat(sprintf("  - Reference: %s\n", names(hiv_ref)[1]))
cat(sprintf("  - Length: %d bp\n", width(hiv_ref)[1]))

# ============================================================================
# STEP 4: Determine viral orientation by alignment to reference
# ============================================================================
cat("\nStep 4: Determining viral orientation by alignment to HIV reference...\n")

# Function to determine orientation by aligning to reference
determine_orientation_by_alignment <- function(viral_seq, ref_seq, min_length = 500) {
  if (is.na(viral_seq) || viral_seq == "" || nchar(viral_seq) < min_length) {
    return(tibble(
      orientation = NA_character_,
      strand = NA_character_,
      alignment_score_fwd = NA_real_,
      alignment_score_rev = NA_real_,
      ref_start = NA_integer_,
      ref_end = NA_integer_,
      percent_identity = NA_real_
    ))
  }
  
  # Convert to DNAString
  query <- DNAString(viral_seq)
  
  # Align forward strand
  aln_fwd <- pwalign::pairwiseAlignment(query, ref_seq, type = "local")
  score_fwd <- score(aln_fwd)
  
  # Align reverse complement
  query_rc <- reverseComplement(query)
  aln_rc <- pwalign::pairwiseAlignment(query_rc, ref_seq, type = "local")
  score_rc <- score(aln_rc)
  
  # Determine orientation based on which alignment is better
  if (score_fwd > score_rc) {
    orientation <- "5'"
    strand <- "+"
    best_aln <- aln_fwd
  } else {
    orientation <- "3'"
    strand <- "-"
    best_aln <- aln_rc
  }

  print(paste0("Best mapping score - ", orientation))
  
  # Get alignment coordinates on reference
  ref_start <- start(subject(best_aln))
  ref_end <- end(subject(best_aln))
  
  tibble(
    orientation = orientation,
    strand = strand,
    alignment_score_fwd = score_fwd,
    alignment_score_rev = score_rc,
    ref_start = ref_start,
    ref_end = ref_end,
    percent_identity = pwalign::pid(best_aln))
}

# Apply orientation determination (this may take a while)
cat("  - Aligning sequences to reference (this may take a few minutes)...\n")
orientation_results <- seq_lookup %>%
  mutate(orientation_data = map(viral_sequence, ~determine_orientation_by_alignment(.x, hiv_ref[[1]])))
seq_lookup <- data.table(orientation_results, rbindlist(orientation_results$orientation_data))
seq_lookup$orientation_data = NULL

#print(head(seq_lookup))

# ============================================================================
# STEP 5: Merge with sequence data
# ============================================================================
cat("\nStep 5: Merging sequence and orientation data...\n")

# Match READ IDs (handle potential /0 suffix mismatch) using tidyverse
df <- df %>%
  mutate(read_match = case_when(READ %in% seq_lookup$read_id ~ READ,
                                paste0(READ, "/0") %in% seq_lookup$read_id ~ paste0(READ, "/0"),
                                str_remove(READ, "/0$") %in% seq_lookup$read_id_base ~ {
        idx <- match(str_remove(READ, "/0$"), seq_lookup$read_id_base)
        seq_lookup$read_id[idx]
      }, TRUE ~ NA_character_)) %>%
  left_join(seq_lookup %>% 
      dplyr::select(read_id, viral_sequence, viral_seq_length, viral_orientation=orientation,
                    viral_strand=strand, alignment_score_fwd, alignment_score_rev,
                    ref_start, ref_end, percent_identity),
      by = c("read_match" = "read_id")) %>%
  select(-read_match)

#print(head(df))

# ============================================================================
# STEP 6: Parse genomic locations
# ============================================================================
cat("\nStep 6: Processing integration sites...\n")

# Parse genomic location using tidyverse
parse_location <- function(loc_string) {
  if (is.na(loc_string) || loc_string == "") {
    return(tibble(chrom = NA_character_, start = NA_integer_, end = NA_integer_))
  }
  
  tryCatch({
    parts <- str_split(loc_string, ":", n = 2)[[1]]
    chrom <- parts[1]
    coords <- str_split(parts[2], "-")[[1]]
    start <- as.integer(coords[1])
    end <- as.integer(coords[2])
    
    tibble(
      chrom = chrom,
      start = min(start, end),
      end = max(start, end)
    )
  }, error = function(e) {
    tibble(chrom = NA_character_, start = NA_integer_, end = NA_integer_)
  })
}

# Get integration site
get_integration_site <- function(left_flank, right_flank) {
  if (!is.na(left_flank) && left_flank != "") {
    loc <- parse_location(left_flank)
    if (!is.na(loc$chrom)) {
      return(sprintf("%s:%d-%d", loc$chrom, loc$start, loc$end))
    }
  }
  
  if (!is.na(right_flank) && right_flank != "") {
    loc <- parse_location(right_flank)
    if (!is.na(loc$chrom)) {
      return(sprintf("%s:%d-%d", loc$chrom, loc$start, loc$end))
    }
  }
  
  return(NA_character_)
}

# Get viral region
get_viral_region <- function(insert_string) {
  if (is.na(insert_string) || insert_string == "") {
    return(NA_character_)
  }
  
  tryCatch({
    parts <- str_split(insert_string, ":")[[1]]
    if (length(parts) >= 2) {
      ref <- paste(parts[-length(parts)], collapse = ":")
      coords <- parts[length(parts)]
      coord_parts <- str_split(coords, "-")[[1]]
      start <- as.integer(coord_parts[1])
      end <- as.integer(coord_parts[2])
      
      pos1 <- min(start, end)
      pos2 <- max(start, end)
      return(sprintf("%s:%d-%d", ref, pos1, pos2))
    }
    return(insert_string)
  }, error = function(e) {
    return(insert_string)
  })
}

# Add integration site information 
df <- df %>%
  mutate(integration_site = map2_chr(LEFT_FLANK, RIGHT_FLANK, get_integration_site),
         viral_region = map_chr(INSERT, get_viral_region),
         chromosome = map_chr(integration_site, ~{
      if (is.na(.x)) return(NA_character_)
      str_split(.x, ":", n = 2)[[1]][1]}),
    sample = basename(input_csv))

# ============================================================================
# STEP 7: Load GTF and annotate genomic features
# ============================================================================
cat("\nStep 7: Loading GTF annotation and annotating genomic features...\n")

# Read GTF using rtracklayer
gtf <- import(gtf_file)
cat(sprintf("  - Loaded %d features from GTF\n", length(gtf)))

# Convert to data frame for easier manipulation
gtf_df <- as.data.frame(gtf)
  
# Create gene-level ranges by taking min/max coordinates per gene
genes_df <- gtf_df %>%
  filter(!is.na(gene_id)) %>%
  group_by(seqnames, gene_id, gene_name, strand) %>%
  summarise(start = min(start),
            end = max(end),.groups = "drop")
  
# Convert back to GRanges
genes <- GRanges(seqnames = genes_df$seqnames,
                ranges = IRanges(start = genes_df$start, end = genes_df$end),
                strand = genes_df$strand,
                gene_id = genes_df$gene_id,
                gene_name = genes_df$gene_name)

# Create GRanges from integration sites - FIXED VERSION
integration_data <- df %>%
  filter(!is.na(integration_site)) %>%
  mutate(row_idx = row_number()) %>%
  rowwise() %>%
  mutate(loc_data = list(parse_location(integration_site))) %>%
  ungroup() %>%
  unnest(loc_data) %>%
  filter(!is.na(chrom))

if (nrow(integration_data) > 0) {
  integration_gr <- GRanges(seqnames = integration_data$chrom,
                            ranges = IRanges(start = integration_data$start,
                                             end = integration_data$end),
                            original_row = integration_data$row_idx)
  
  if (length(integration_gr) > 0) {
    
    # Find overlaps with genes
    overlaps <- findOverlaps(integration_gr, genes)
    
    if (length(overlaps) > 0) {
      # Create annotation mapping using tidyverse
      ann <- tibble(query_idx = queryHits(overlaps),
                    original_row = integration_gr$original_row[queryHits(overlaps)],
                    gene_name = genes$gene_name[subjectHits(overlaps)],
                    gene_id = genes$gene_id[subjectHits(overlaps)])
      
      # Aggregate multiple genes per integration site
      ann_agg <- ann %>%
        group_by(original_row) %>%
        summarise(gene_name = paste(unique(na.omit(gene_name)), collapse = ";"),
                  gene_id = paste(unique(na.omit(gene_id)), collapse = ";"),
                  .groups = "drop") %>%
        mutate(gene_name = if_else(gene_name == "", NA_character_, gene_name),
               gene_id = if_else(gene_id == "", NA_character_, gene_id))
      
      # Add annotations to dataframe using left_join
      df <- df %>%
        mutate(row_idx = row_number()) %>%
        left_join(ann_agg, by = c("row_idx" = "original_row")) %>%
        select(-row_idx)
      
      cat(sprintf("  - Annotated %d reads with gene information\n", 
                  sum(!is.na(df$gene_name))))
    } else {
      cat("  - No overlaps found between integration sites and genes\n")
      df <- df %>%
        mutate(gene_name = NA_character_,
               gene_id = NA_character_)
    }
  }
} else {
  cat("  - No valid integration sites to annotate\n")
  df <- df %>%
    mutate(gene_name = NA_character_,
           gene_id = NA_character_)
}

# Sequence length summary
cat("\n", strrep("-", 70), "\n")
cat("VIRAL SEQUENCE LENGTH SUMMARY\n")
cat(strrep("-", 70), "\n")

length_stats <- df %>%
  filter(!is.na(viral_seq_length), viral_seq_length > 0) %>%
  summarise(
    n = n(),
    mean_length = mean(viral_seq_length),
    median_length = median(viral_seq_length),
    min_length = min(viral_seq_length),
    max_length = max(viral_seq_length),
    sd_length = sd(viral_seq_length)
  )

if (length_stats$n > 0) {
  cat(sprintf("Number of sequences: %d\n", length_stats$n))
  cat(sprintf("Mean length: %.1f bp\n", length_stats$mean_length))
  cat(sprintf("Median length: %.0f bp\n", length_stats$median_length))
  cat(sprintf("Range: %.0f - %.0f bp\n", length_stats$min_length, length_stats$max_length))
  cat(sprintf("Std dev: %.1f bp\n", length_stats$sd_length))
}

# Alignment quality summary
cat("\n", strrep("-", 70), "\n")
cat("ALIGNMENT QUALITY SUMMARY\n")
cat(strrep("-", 70), "\n")

alignment_stats <- df %>%
  filter(!is.na(percent_identity)) %>%
  summarise(n = n(),
            mean_pid = mean(percent_identity),
            median_pid = median(percent_identity),
            min_pid = min(percent_identity),
            max_pid = max(percent_identity))

if (alignment_stats$n > 0) {
  cat(sprintf("Mean percent identity to reference: %.1f%%\n", alignment_stats$mean_pid))
  cat(sprintf("Median percent identity: %.1f%%\n", alignment_stats$median_pid))
  cat(sprintf("Range: %.1f%% - %.1f%%\n", alignment_stats$min_pid, alignment_stats$max_pid))
}

# Phylogenetic and similarity analysis
cat("\n", strrep("-", 70), "\n")
cat("SEQUENCE SIMILARITY ANALYSIS\n")
cat(strrep("-", 70), "\n")

# Get unique sequences
unique_seqs <- df %>%
  filter(!is.na(viral_sequence), viral_seq_length >= 100) %>%
  distinct(viral_sequence, .keep_all = TRUE)

cat(sprintf("Number of unique viral sequences (>= 100 bp): %d\n", nrow(unique_seqs)))

if (nrow(unique_seqs) >= 2) {
  cat("\nCalculating pairwise sequence similarities...\n")
  
  # Create DNAStringSet for viral sequences
  viral_seqs <- DNAStringSet(unique_seqs$viral_sequence)
  names(viral_seqs) <- paste0("seq", 1:length(viral_seqs), "_",
                              unique_seqs$viral_orientation, "_",
                              str_sub(unique_seqs$READ, 1, 30))
  
  # Write viral sequences to FASTA
  viral_fasta_file <- paste0(output_prefix, "_viral_sequences.fasta")
  writeXStringSet(viral_seqs, viral_fasta_file)
  cat(sprintf("  - Saved viral sequences to: %s\n", viral_fasta_file))
  
  # Calculate pairwise alignment similarities
  n_seqs <- length(viral_seqs)
  similarity_matrix <- matrix(0, nrow = n_seqs, ncol = n_seqs)
  rownames(similarity_matrix) <- names(viral_seqs)
  colnames(similarity_matrix) <- names(viral_seqs)
  
  # Calculate similarities
  for (i in 1:n_seqs) {
    similarity_matrix[i, i] <- 1.0
    if (i < n_seqs) {
      for (j in (i+1):n_seqs) {
        # Simple pairwise alignment
        aln <- pairwiseAlignment(viral_seqs[i], viral_seqs[j])
        pid <- pwalign::pid(aln) / 100  # Convert to proportion
        similarity_matrix[i, j] <- pid
        similarity_matrix[j, i] <- pid
      }
    }
  }
  
  # Save similarity matrix
  sim_file <- paste0(output_prefix, "_similarity_matrix.csv")
  write_csv(as_tibble(similarity_matrix, rownames = "sequence"), sim_file)
  
  # Summary statistics
  upper_tri <- similarity_matrix[upper.tri(similarity_matrix)]
  
  # Create heatmap
  if (n_seqs <= 50) {
    heatmap_file <- paste0(output_prefix, "_similarity_heatmap.pdf")
    pdf(heatmap_file, width = 12, height = 10)
    heatmap(similarity_matrix, 
            main = "Viral Sequence Similarity Matrix",
            xlab = "Sequences", ylab = "Sequences",
            scale = "none",
            col = colorRampPalette(c("blue", "white", "red"))(100),
            margins = c(10, 10))
    dev.off()
  }
  
  # Build phylogenetic tree
  if (n_seqs >= 3 && n_seqs <= 100) {
    cat("\nBuilding phylogenetic tree...\n")
    
    # Convert to distance matrix
    dist_mat <- as.dist(1 - similarity_matrix)
    
    # Build neighbor-joining tree
    tree <- nj(dist_mat)
    
    # Save tree
    tree_file <- paste0(output_prefix, "_tree.nwk")
    write.tree(tree, tree_file)
    
    # Plot tree
    tree_pdf <- paste0(output_prefix, "_tree.pdf")
    pdf(tree_pdf, width = 12, height = max(8, n_seqs * 0.4))
    plot(tree, 
         cex = 0.7,
         main = "Phylogenetic Tree of Viral Sequences (Neighbor-Joining)",
         sub = "Distance based on sequence similarity")
    add.scale.bar()
    dev.off()
  }
}

# Summary statistics
cat("\n", strrep("-", 70), "\n")
cat("INTEGRATION SITE STATISTICS\n")
cat(strrep("-", 70), "\n")

cat(sprintf("Total reads: %d\n", nrow(df)))
cat(sprintf("Unique integration sites: %d\n", 
            n_distinct(df$integration_site, na.rm = TRUE)))

if ("gene_name" %in% colnames(df)) {
  cat(sprintf("Reads in annotated genes: %d\n", sum(!is.na(df$gene_name))))
  unique_genes <- df %>%
    filter(!is.na(gene_name)) %>%
    pull(gene_name) %>%
    str_split(";") %>%
    unlist() %>%
    unique()
  cat(sprintf("Unique genes hit: %d\n", length(unique_genes)))
}

# Chromosome distribution
cat("\nTop 10 chromosomes by integration count:\n")
chr_summary <- df %>%
  filter(!is.na(chromosome)) %>%
  count(chromosome, sort = TRUE) %>%
  head(10)

for (i in 1:nrow(chr_summary)) {
  cat(sprintf("  %s: %d\n", chr_summary$chromosome[i], chr_summary$n[i]))
}

# Save annotated output
cat("\n", strrep("-", 70), "\n")
cat("SAVING RESULTS\n")
cat(strrep("-", 70), "\n")

output_csv <- paste0(output_prefix, "_annotated.csv")
write_csv(df, output_csv)

# Write summary report
report_file <- paste0(output_prefix, "_report.txt")
sink(report_file)

cat("VIRAL ORIENTATION AND GENOMIC ANNOTATION REPORT\n")
cat(strrep("=", 70), "\n\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input CSV:", input_csv, "\n")
cat("Input FASTA:", input_fasta, "\n")
cat("HIV Reference:", hiv_ref_fasta, "\n")
cat("GTF Annotation:", gtf_file, "\n\n")

cat("ORIENTATION SUMMARY\n")
cat(strrep("-", 70), "\n")

cat("\n5' = Forward strand (same orientation as reference)\n")
cat("3' = Reverse strand (reverse complement of reference)\n")

if (length_stats$n > 0) {
  cat("\n\nSEQUENCE LENGTH SUMMARY\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Number of sequences: %d\n", length_stats$n))
  cat(sprintf("Mean length: %.1f bp\n", length_stats$mean_length))
  cat(sprintf("Median length: %.0f bp\n", length_stats$median_length))
  cat(sprintf("Range: %.0f - %.0f bp\n", length_stats$min_length, length_stats$max_length))
  cat(sprintf("Std dev: %.1f bp\n", length_stats$sd_length))
}

if (alignment_stats$n > 0) {
  cat("\n\nALIGNMENT QUALITY\n")
  cat(strrep("-", 70), "\n")
  cat(sprintf("Mean percent identity to reference: %.1f%%\n", alignment_stats$mean_pid))
  cat(sprintf("Median percent identity: %.1f%%\n", alignment_stats$median_pid))
  cat(sprintf("Range: %.1f%% - %.1f%%\n", alignment_stats$min_pid, alignment_stats$max_pid))
}

cat("\n\nINTEGRATION SITE STATISTICS\n")
cat(strrep("-", 70), "\n")
cat(sprintf("Total reads: %d\n", nrow(df)))
cat(sprintf("Unique integration sites: %d\n", 
            n_distinct(df$integration_site, na.rm = TRUE)))
if ("gene_name" %in% colnames(df)) {
  cat(sprintf("Reads in annotated genes: %d\n", sum(!is.na(df$gene_name))))
  cat(sprintf("Unique genes hit: %d\n", length(unique_genes)))
}

cat("\n\nMETHODOLOGY\n")
cat(strrep("-", 70), "\n")
cat("Orientation Determination:\n")
cat("  1. Viral sequences extracted from unmasked FASTA (N = host, ATCG = viral)\n")
cat("  2. Each viral sequence aligned to HIV reference genome (both strands)\n")
cat("  3. Orientation assigned based on which strand gives better alignment\n")
cat("  4. 5' orientation = forward strand (same as reference)\n")
cat("  5. 3' orientation = reverse complement strand\n\n")
cat("Genomic Annotation:\n")
cat("  1. Integration sites extracted from LEFT_FLANK and RIGHT_FLANK fields\n")
cat("  2. Sites overlapped with genes from GTF annotation\n")
cat("  3. Multiple overlapping genes are semicolon-separated\n")

sink()
cat(sprintf("  - Saved report to: %s\n", report_file))

cat("\n", strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE!\n")
cat(strrep("=", 70), "\n\n")

cat("Output files:\n")
cat("  1.", output_csv, "\n")
if (exists("viral_fasta_file")) cat("  2.", viral_fasta_file, "\n")
if (exists("sim_file")) cat("  3.", sim_file, "\n")
if (exists("heatmap_file")) cat("  4.", heatmap_file, "\n")
if (exists("tree_file")) cat("  5.", tree_file, "\n")
if (exists("tree_pdf")) cat("  6.", tree_pdf, "\n")
cat("  7.", report_file, "\n\n")
