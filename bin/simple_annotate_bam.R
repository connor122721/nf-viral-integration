#!/usr/bin/env Rscript
#
# Annotate integration sites for PCR duplicates, clonal IDs, and genomic features
# PCR duplicates = same integration site + same viral insert region + same shear sites +/- 5bps
# Clonal ID = unique integration site (regardless of viral sequence) + shear sites +/- 5bps

# Libraries
library(tidyverse)
library(data.table)
library(GenomicRanges)
library(rtracklayer)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("Usage: annotate_integration_sites.R <input.csv> <annotation.gff3> [output.csv]\n")
  cat("\nExample:\n")
  cat("  Rscript annotate_integration_sites.R data.csv chm13v2.0_RefSeq_Liftoff_v5.1.gff3.gz output.csv\n")
  quit(status = 1)
}

# Input files argument
input_file <- args[1]
gtf_file <- args[2]
output_file <- if (length(args) >= 3) args[3] else sub("\\.csv$", "_annotated.csv", input_file)

# TESTING
# input_file="Z:/nf-viral-integration_hg38/output/04_final_results/m84248_250919_192326_s1.hifi_reads.univ_v3_bc1002.trim.csv"; gtf_file="Z:/human_genome/t2t/chm13v2.0_RefSeq_Liftoff_v5.2.gtf"
# input_file="~/nf-viral-integration_t2t/output/04_final_results/m84248_250919_192326_s1.hifi_reads.univ_v3_bc1002.csv"; gtf_file="~/human_genome/t2t/chm13v2.0_RefSeq_Liftoff_v5.2.gff3"

# Read in summary
df <- data.table(fread(input_file) %>% 
                   filter(INSERT_LEN >= 200))

# Parse genomic location from LEFT_FLANK or RIGHT_FLANK
parse_location <- function(loc_string) {
  # loc_string="chr9:131412592-131413034"
  
  if (is.na(loc_string) || loc_string == "") {
    return(list(chrom = NA, start = NA, end = NA))
  }
  
  tryCatch({
    parts <- str_split(loc_string, ":", n = 2)[[1]]
    chrom <- parts[1]
    coords <- str_split(parts[2], "-")[[1]]
    start <- as.integer(coords[1])
    end <- as.integer(coords[2])
    
    # Normalize to start < end
    list(
      chrom = chrom,
      start = min(start, end),
      end = max(start, end)
    )
  }, error = function(e) {
    list(chrom = NA, start = NA, end = NA)
  })
}

# Extract integration site
get_integration_site <- function(left_flank, right_flank) {
  # Try left flank first
  if (!is.na(left_flank) && left_flank != "") {
    loc <- parse_location(left_flank)
    if (!is.na(loc$chrom)) {
      return(sprintf("%s:%d-%d", loc$chrom, loc$start, loc$end))
    }
  }
  
  # Try right flank
  if (!is.na(right_flank) && right_flank != "") {
    loc <- parse_location(right_flank)
    if (!is.na(loc$chrom)) {
      return(sprintf("%s:%d-%d", loc$chrom, loc$start, loc$end))
    }
  }
  
  return(NA_character_)
}

# Extract and normalize viral region from INSERT field
get_viral_region <- function(insert_string) {
  if (is.na(insert_string) || insert_string == "") {
    return(NA_character_)
  }
  
  tryCatch({
    # Split on last colon to get coordinates
    parts <- str_split(insert_string, ":")
    if (length(parts[[1]]) >= 2) {
      ref <- paste(parts[[1]][-length(parts[[1]])], collapse = ":")
      coords <- parts[[1]][length(parts[[1]])]
      coord_parts <- str_split(coords, "-")[[1]]
      start <- as.integer(coord_parts[1])
      end <- as.integer(coord_parts[2])
      
      # Normalize coordinates
      pos1 <- min(start, end)
      pos2 <- max(start, end)
      return(sprintf("%s:%d-%d", ref, pos1, pos2))
    }
    return(insert_string)
  }, error = function(e) {
    return(insert_string)
  })
}

cat("Extracting integration sites and viral regions...\n")

# Add integration site and viral region columns
df <- df %>%
  rowwise() %>%
  mutate(integration_site = get_integration_site(LEFT_FLANK, RIGHT_FLANK),
         viral_region = get_viral_region(INSERT),
         chromosome = if_else(!is.na(integration_site), str_split(integration_site, ":", n = 2)[[1]][1], NA_character_),
         sample = unique(basename(input_file))) %>%
  ungroup()
dt <- as.data.table(df)

# Parse LEFT_FLANK coordinates for shear sites
dt[, c('Left_Prior', 'fivePrimeInsertion') := tstrsplit(LEFT_FLANK, '-', fixed = TRUE)]
dt[, c('Left_Omit', 'fivePrimeShear') := tstrsplit(Left_Prior, ':', fixed = TRUE)]

# Parse RIGHT_FLANK coordinates for shear sites
dt[, c('Right_Prior', 'threePrimeShear') := tstrsplit(RIGHT_FLANK, '-', fixed = TRUE)]
dt[, c('Right_Omit', 'threePrimeInsertion') := tstrsplit(Right_Prior, ':', fixed = TRUE)]

# Convert to numeric
dt[, `:=`(
  fivePrimeShear = as.numeric(fivePrimeShear),
  fivePrimeInsertion = as.numeric(fivePrimeInsertion),
  threePrimeInsertion = as.numeric(threePrimeInsertion),
  threePrimeShear = as.numeric(threePrimeShear))]

# Remove temporary columns
dt[, c('Left_Prior', 'Left_Omit', 'Right_Prior', 'Right_Omit') := NULL]

# Name flanking of viral integration
dt1 <- data.table(dt %>% 
   mutate(flanking = case_when(
        !LEFT_FLANK=="" & !RIGHT_FLANK=="" ~ "Dually-Flanked",
        !LEFT_FLANK=="" & RIGHT_FLANK=="" ~ "5'-Flanked",
        LEFT_FLANK=="" & !RIGHT_FLANK=="" ~ "3'-Flanked",
        TRUE ~ "Non-Flanked")))

# Set tolerance for integration site grouping
tolerance <- 5

# Sort by chromosome and positions
setorder(dt, flanking, chromosome, viral_region, fivePrimeInsertion, threePrimeInsertion, na.last = TRUE)

# Parameters
tolerance <- 5
sample_name <- tools::file_path_sans_ext(basename(input_file))

cat("Assigning clonal IDs...\n")

# Convert to tibble for easier manipulation
df <- as_tibble(dt)

# Assign clonal IDs
df <- df %>%
  # Sort data
  arrange(chromosome, viral_region, fivePrimeInsertion, threePrimeInsertion) %>%
  
  # Group integration sites within tolerance
  group_by(chromosome, viral_region) %>%
  mutate(
    integration_group = if_else(
      is.na(fivePrimeInsertion) | is.na(threePrimeInsertion),
      NA_integer_,
      cumsum(
        row_number() == 1 |
        abs(fivePrimeInsertion - lag(fivePrimeInsertion, default = -9999)) > tolerance |
        abs(threePrimeInsertion - lag(threePrimeInsertion, default = -9999)) > tolerance
      )
    )
  ) %>%
  ungroup() %>%
  
  # Number different shear patterns within each integration group
  group_by(chromosome, viral_region, integration_group) %>%
  mutate(
    clone_number = dense_rank(paste(fivePrimeShear, threePrimeShear))
  ) %>%
  ungroup() %>%
  
  # Create clonal ID
  mutate(
    clonal_id = if_else(
      !is.na(integration_group),
      paste(sample_name, viral_region, integration_group, clone_number, sep = "_"),
      NA_character_
    )
  )

# Convert back to data.table
dt <- as.data.table(df)

# Report
cat("Assigned", n_distinct(dt$clonal_id, na.rm = TRUE), "unique clonal IDs\n")







dt <- as_tibble(dt) %>%
  # Sort data
  arrange(flanking, chromosome, viral_region, fivePrimeInsertion, threePrimeInsertion) %>%
  
  # Assign integration groups (only for valid rows)
  group_by(chromosome, viral_region) %>%
  mutate(
    integration_group = if_else(
      is.na(fivePrimeInsertion) | is.na(threePrimeInsertion),
      NA_integer_,
      cumsum(
        row_number() == 1 |
        abs(fivePrimeInsertion - lag(fivePrimeInsertion, default = -999)) > tolerance |
        abs(threePrimeInsertion - lag(threePrimeInsertion, default = -999)) > tolerance
      )
    )
  ) %>%
  ungroup() %>%
  
  # Assign clone numbers based on shear sites
  group_by(chromosome, viral_region, integration_group) %>%
  mutate(
    clone_number = dense_rank(paste(fivePrimeShear, threePrimeShear))
  ) %>%
  ungroup() %>%
  
  # Create clonal ID
  mutate(
    clonal_id = if_else(
      !is.na(integration_group),
      paste(sample, viral_region, integration_group, clone_number, sep = "_"),
      NA_character_
    )
  ) %>%
  
  # Clean up temporary columns
  select(-integration_group, -clone_number) %>%
  
  # Convert back to data.table
  as.data.table()




# Create integration groups with tolerance
dt[, integration_group := {
  n <- .N
  groups <- integer(n)
  
  if (n == 0) {
    return(groups)
  }
  
  # Initialize first group
  groups[1] <- 1
  current_group <- 1
  
  if (n > 1) {
    for (i in 2:n) {
      # Handle NA values - each NA gets its own group
      if (is.na(fivePrimeInsertion[i]) || is.na(threePrimeInsertion[i]) ||
          is.na(fivePrimeInsertion[i-1]) || is.na(threePrimeInsertion[i-1])) {
        current_group <- current_group + 1
        groups[i] <- current_group
      } else if (abs(fivePrimeInsertion[i] - fivePrimeInsertion[i-1]) <= tolerance &&
                 abs(threePrimeInsertion[i] - threePrimeInsertion[i-1]) <= tolerance) {
        # Within tolerance - same group
        groups[i] <- current_group
      } else {
        # Outside tolerance - new group
        current_group <- current_group + 1
        groups[i] <- current_group
      }
    }
  }
  
  groups  # Return the groups vector
}, by = .(chromosome, viral_region)]

# Assign clone numbers within each integration group
dt[, clone_number := {
  if (.N == 0 || all(is.na(fivePrimeShear)) || all(is.na(threePrimeShear))) {
    rep(NA_integer_, .N)
  } else {
    frank(paste(fivePrimeShear, threePrimeShear), ties.method = "dense", na.last = "keep")
  }
}, by = .(chromosome, viral_region, integration_group)]

# Create final clonal ID
dt[, clonal_id := fifelse(
  is.na(integration_group) | is.na(clone_number),
  NA_character_,
  paste(sample, viral_region, integration_group, clone_number, sep = "_")
)]

# Clean up temporary columns
dt[, c('integration_group', 'clone_number') := NULL]



# Assign clone numbers within each integration group (only for non-NA)
dt[, clone_number := {
  if (all(is.na(fivePrimeShear)) || all(is.na(threePrimeShear))) {
    rep(NA_integer_, .N)
  } else {
    frank(paste(fivePrimeShear, threePrimeShear), ties.method = "dense", na.last = "keep")
  }
}, by = .(chromosome, viral_region, integration_group)]

# Create final clonal ID
dt[, clonal_id := fcase(
  is.na(integration_site), NA_character_,
  is.na(fivePrimeInsertion) | is.na(threePrimeInsertion), NA_character_,
  default = paste(sample_name, viral_region, integration_group, clone_number, sep = "_")
)]

# Clean up
dt[, c('integration_group', 'clone_number') := NULL]





# Assign clonal IDs (unique integration sites + shear sites +/- 5bps)
cat("Assigning clonal IDs...\n")
dt1 <- data.table(dt %>%
  group_by(integration_site) %>%
  mutate(
    # Check if shear sites differ within same integration site
    n_unique_shears = n_distinct(paste(fivePrimeShear, threePrimeShear), na.rm = TRUE),
    # Create sequential clone number within each integration site
    clone_number = paste0("N", dense_rank(paste(fivePrimeShear, threePrimeShear)))) %>%
  ungroup() %>%
  mutate(clonal_id = if_else(!is.na(integration_site),
                             paste(sample, viral_region, clone_number, sep = "_"), NA_character_)) %>%
  select(-clone_number))


dt1 <- data.table(dt %>%
  group_by(integration_site) %>%
  mutate(clonal_id = cur_group_id()) %>%
  ungroup() %>%
  mutate(clonal_id = if_else(is.na(integration_site), NA_integer_, clonal_id)))

# Identify PCR duplicates (same integration site + same viral region + same shear sites)
cat("Identifying PCR duplicates...\n")
dt1 <- data.table(dt1 %>%
  group_by(integration_site, viral_region) %>%
  mutate(pcr_duplicate_group = cur_group_id(),
         n_pcr_duplicates = n(),
         is_pcr_duplicate = n() > 1) %>%
  ungroup() %>%
  mutate(pcr_duplicate_group = if_else(is.na(integration_site) | is.na(viral_region), 
      NA_integer_, 
      pcr_duplicate_group),
    is_pcr_duplicate = if_else(
      is.na(integration_site) | is.na(viral_region),
      NA,
      is_pcr_duplicate)))

# Load GTF/GFF and annotate genomic features
cat("Loading GTF/GFF annotation:", gtf_file, "\n")
gtf <- import(gtf_file)

# Filter for relevant features (genes, exons, etc.)
# Adjust feature types based on your GTF format
genes <- gtf[gtf$type == "gene"]

cat("Annotating genomic features...\n")

# Create GRanges from integration sites
integration_gr <- dt1 %>%
  filter(!is.na(integration_site)) %>%
  rowwise() %>%
  mutate(loc = list(parse_location(integration_site)),
         chrom = loc$chrom,
         start = loc$start,
         end = loc$end) %>%
  filter(!is.na(chrom)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE,
                           seqnames.field = "chrom",
                           start.field = "start",
                           end.field = "end")

# Find overlaps with genes
if (length(integration_gr) > 0) {
  overlaps <- findOverlaps(integration_gr, genes)
  
  # Create annotation mapping
  ann <- tibble(query_idx = queryHits(overlaps),
                gene_name = genes$gene_name[subjectHits(overlaps)],
                gene_id = genes$ID[subjectHits(overlaps)],
                gene_type = genes$gene_biotype[subjectHits(overlaps)]) %>%
    group_by(query_idx) %>%
    summarize(gene_name = paste(unique(gene_name), collapse = ";"),
              gene_id = paste(unique(gene_id), collapse = ";"),
              gene_type = paste(unique(gene_type), collapse = ";"))
  
  # Map back to original dataframe indices
  integration_gr$original_idx <- which(!is.na(dt1$integration_site) & 
                                         !is.na(dt1 %>% 
                                                  rowwise() %>% 
                                                  mutate(loc = list(parse_location(integration_site))) %>% 
                                                  pull(loc) %>% 
                                                  sapply(function(x) x$chrom)))
  
  # Add annotations to dataframe
  dt1 <- data.table(dt1 %>%
    mutate(gene_name = NA_character_,
           gene_id = NA_character_,
           gene_type = NA_character_))
  
  for (i in seq_len(nrow(ann))) {
    orig_idx <- integration_gr$original_idx[ann$query_idx[i]]
    dt1$gene_name[orig_idx] <- ann$gene_name[i]
    dt1$gene_id[orig_idx] <- ann$gene_id[i]
    dt1$gene_type[orig_idx] <- ann$gene_type[i]
  }
}

# Summary statistics
cat("\n=== Summary Statistics ===\n")
cat("Total reads:", nrow(dt1), "\n")
cat("Unique integration sites (clones):", n_distinct(dt1$clonal_id, na.rm = TRUE), "\n")
cat("Reads with PCR duplicates:", sum(dt1$is_pcr_duplicate, na.rm = TRUE), "\n")
cat("PCR duplicate groups:", n_distinct(dt1$pcr_duplicate_group, na.rm = TRUE), "\n")
if ("gene_name" %in% colnames(dt1)) {
  cat("Reads in annotated genes:", sum(!is.na(dt1$gene_name)), "\n")
  cat("Unique genes hit:", n_distinct(dt1$gene_name[!is.na(dt1$gene_name)]), "\n")
}

# Write output
cat("\nWriting output to:", output_file, "\n")
write_csv(dt1, output_file)

cat("Done!\n")