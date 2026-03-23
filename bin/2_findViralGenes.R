#!/usr/bin/env Rscript

############################################
## AUTHOR (original Perl):                ##
##               Eric C. Rouchka          ##
##               University of Louisville ##
## LAST UPDATED: 6/5/2024                 ##
##                                        ##
## R CONVERSION + TIDYVERSE REFACTOR:     ##
##               Connor S. Murray         ##
##               University of Louisville ##
##                                        ##
## SCRIPT GOAL:                           ##
##   Annotate the completeness of HIV/SIV ##
##   proviral genomes from long-read HiFi ##
##   sequencing. For each CCS read, BLAST ##
##   identifies which viral gene segments ##
##   are present (LTR, GAG, POL, VIF,     ##
##   VPR, VPU, ENV, NEF) and what         ##
##   fraction of each gene is recovered.  ##
##   A 9-character match string encodes   ##
##   coverage across the 5'LTR→3'LTR      ##
##   axis (0=absent, 1/2=partial,         ##
##   3=full), and the read is classified  ##
##   as: INTACT, PUTATIVELY INTACT,       ##
##   INDETERMINATE, TRUNCATED, or         ##
##   INTERNAL DELETION. IPDA probe hits   ##
##   (PSI, RRE, LTR-GAG, POL, ENV) are   ##
##   scored in parallel. An EPISOME flag  ##
##   is raised when gene ordering is      ##
##   inconsistent with a linear provirus, ##
##   suggesting a 2-LTR circle.           ##
##                                        ##
## V8: Added code to save gene sequences  ##
## V5: Added 0,1,2,3 designation          ##
##     Added length histograms            ##
## V4: Added IPDA probes (PMC9525950)     ##
############################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

########################################################
## GLOBAL VARIABLES

DB_BASE <- "data/HIV_GENES/BLAST_DB"
BLAST_CMD <- "blastn"
TMPDIR <- "tmp"
IPDA_BLASTDB <- "data/HIV_GENES/BLAST_DB/HIV/IPDA/IPDA"
IPDA_V2_BLASTDB <- "data/HIV_GENES/BLAST_DB/HIV/IPDA/IPDA_PMC9525950"
PCT_CUTOFF <- 0.7    # threshold for val=3 (full coverage)
PCT_CUTOFF_2 <- 0.5    # threshold for val=2 (partial coverage)
MIN_MATCH_LEN <- 30    # minimum BLAST hit length to retain

BLAST_COLS <- c("qseqid","qstart","qend","sseqid","sstrand",
                "sstart","send","evalue","bitscore","score",
                "length","pident","gaps")

########################################################
## ARGUMENT PARSING

args <- commandArgs(trailingOnly = TRUE)

parse_arg <- function(flag, args, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0 || idx == length(args)) return(default)
  args[idx + 1]
}

inFN      <- parse_arg("--in",        args)
OUT_FN    <- parse_arg("--out",       args)
reference <- parse_arg("--reference", args)
virusType <- parse_arg("--virus",     args)
PREFIX    <- parse_arg("--prefix",    args, default = "")

if (!is.null(PREFIX) && nchar(PREFIX) > 0) {
  PREFIX <- str_remove(PREFIX, "/$")
  DB_BASE <- file.path(PREFIX, DB_BASE)
  TMPDIR <- file.path(PREFIX, TMPDIR)
  IPDA_BLASTDB <- file.path(PREFIX, IPDA_BLASTDB)
  IPDA_V2_BLASTDB <- file.path(PREFIX, IPDA_V2_BLASTDB)
}

########################################################
## UTILITY FUNCTIONS

print_usage <- function() {
  cat("USAGE: --in <IN_FILE> --out <OUT_FILE> --reference <REF> --virus <HIV|SIV>\n")
  quit(status = 0)
}

get_tmp_path <- function(extension) {
  repeat {
    fn <- file.path(TMPDIR, str_c(as.numeric(Sys.time()), "_", runif(1), extension))
    if (!file.exists(fn)) return(fn)
  }
}

safe_rm <- function(fn) { if (file.exists(fn)) invisible(file.remove(fn)) }

reverse_complement <- function(s) {
  s %>% toupper() %>% chartr("ATCG", "TAGC", x = .) %>%
    str_split("") %>% pluck(1) %>% rev() %>% str_c(collapse = "")
}

overall_strand <- function(n_minus, n_plus) if (n_minus > n_plus) "minus" else "plus"

########################################################
## VALIDATE PARAMETERS

validate_parameters <- function(inFN, OUT_FN, reference, virusType) {
  if (any(is.null(c(inFN, OUT_FN, reference, virusType)))) print_usage()
  if (!file.exists(inFN))                  stop("Input file does not exist: ", inFN)
  if (!virusType %in% c("HIV", "SIV"))     stop("--virus must be HIV or SIV")
  pat     <- file.path(DB_BASE, virusType, reference, "allGenes.fa")
  matches <- Sys.glob(pat)
  if (length(matches) == 0) stop("No BLAST database found matching: ", pat)
  if (length(matches) >  1) stop("Multiple databases match: ",         pat)
  str_remove(matches[1], "\\.fa$")
}

########################################################
## READ FASTA  →  tibble(header, seq)

read_fasta <- function(path) {
  lines   <- read_lines(path)
  hdr_idx <- which(str_starts(lines, ">"))
  seq_end <- c(hdr_idx[-1] - 1L, length(lines))
  tibble(
    header = lines[hdr_idx],
    seq    = map2_chr(hdr_idx + 1L, seq_end,
                      ~str_c(str_squish(lines[.x:.y]), collapse = ""))
  )
}

########################################################
## GET REFERENCE SEQUENCES  →  named list (gene → FASTA text)

get_reference_sequences <- function(blast_fa) {
  read_fasta(blast_fa) %>%
    mutate(gene = str_remove(str_extract(header, "^>[^_]+"), "^>")) %>%
    group_by(gene) %>%
    summarise(fasta_text = str_c(str_c(header, "\n", seq, "\n"), collapse = ""),
              .groups = "drop") %>%
    deframe()   # named character vector / list
}

########################################################
## GET GENE LENGTHS  →  named integer vector

get_gene_lengths <- function(reference, virusType) {
  genes <- c("LTR","GAG","POL","VIF","VPR","VPU","ENV","NEF")
  map_int(set_names(genes), function(gene) {
    pat <- file.path(DB_BASE, virusType, reference, str_c(gene, ".fa"))
    f   <- Sys.glob(pat)
    if (length(f) == 0) stop("Cannot find gene file: ", pat)
    fa  <- read_fasta(f[1])
    nchar(fa$seq[1])
  })
}

########################################################
## GET SEQUENCE COUNT

get_sequence_count <- function(path) {
  read_lines(path) %>% str_starts(">") %>% sum()
}

########################################################
## RUN BLAST  →  data.table

run_blast <- function(query_fa, db, word_size = 12,
                      gapopen = 5, gapextend = 2, reward = 2, penalty = -3) {
  out_fn <- get_tmp_path(".blast")
  cmd <- str_glue(
    "{BLAST_CMD} -db {shQuote(db)} -query {shQuote(query_fa)} \\
     -word_size {word_size} -gapopen {gapopen} -gapextend {gapextend} \\
     -reward {reward} -penalty {penalty} \\
     -outfmt '6 qseqid qstart qend sseqid sstrand sstart send \\
              evalue bitscore score length pident gaps' \\
     > {shQuote(out_fn)}"
  )
  system(cmd)
  dt <- if (file.exists(out_fn) && file.info(out_fn)$size > 0)
    fread(out_fn, header = FALSE, sep = "\t", col.names = BLAST_COLS)
  else
    data.table()[, (BLAST_COLS) := character()]
  safe_rm(out_fn)
  dt
}

########################################################
## FIND IPDA (PSI / RRE)

find_ipda <- function(query_fa) {
  out_fn <- get_tmp_path(".IPDA.blast")
  cmd <- str_glue(
    "{BLAST_CMD} -db {shQuote(IPDA_BLASTDB)} -query {shQuote(query_fa)} \\
     -word_size 11 \\
     -outfmt '6 qseqid qstart qend sseqid sstrand sstart send \\
              evalue bitscore score length pident gaps' \\
     > {shQuote(out_fn)}"
  )
  system(cmd)
  dt <- if (file.exists(out_fn) && file.info(out_fn)$size > 0)
    fread(out_fn, header = FALSE, sep = "\t", col.names = BLAST_COLS)
  else data.table()
  safe_rm(out_fn)
  list(
    psi = as.integer(any(dt$sseqid == "PSI")),
    rre = as.integer(any(dt$sseqid == "RRE"))
  )
}

########################################################
## FIND IPDA V2 (LTR_GAG / POL / ENV)
## Scoring: 0=not found, 1=FWD only, 2=REV only, 3=both strands

ipda_v2_score <- function(regions, fwd_tag, rev_tag) {
  has_fwd <- fwd_tag %in% regions
  has_rev <- rev_tag %in% regions
  case_when(
    has_fwd & has_rev  ~ 3L,
    has_fwd & !has_rev ~ 1L,
    !has_fwd & has_rev ~ 2L,
    TRUE               ~ 0L
  )
}

find_ipda_v2 <- function(query_fa) {
  out_fn <- get_tmp_path(".IPDA_V2.blast")
  cmd <- str_glue(
    "{BLAST_CMD} -db {shQuote(IPDA_V2_BLASTDB)} -query {shQuote(query_fa)} \\
     -word_size 9 \\
     -outfmt '6 qseqid qstart qend sseqid sstrand sstart send \\
              evalue bitscore score length pident gaps' \\
     > {shQuote(out_fn)}"
  )
  system(cmd)
  dt <- if (file.exists(out_fn) && file.info(out_fn)$size > 0)
    fread(out_fn, header = FALSE, sep = "\t", col.names = BLAST_COLS)
  else data.table()
  safe_rm(out_fn)
  regions <- dt$sseqid
  list(
    ltr_gag = ipda_v2_score(regions, "LTR_FWD", "LTR_REV"),
    pol     = ipda_v2_score(regions, "POL_FWD", "POL_REV"),
    env     = ipda_v2_score(regions, "ENV_FWD", "ENV_REV")
  )
}

########################################################
## MERGE BLAST HITS FOR ONE GENE
## Collapses overlapping / near-adjacent (<100 bp gap) intervals.

merge_hits <- function(beg_vec, end_vec) {
  # Ensure all intervals are [min, max]
  dt <- data.table(beg = pmin(beg_vec, end_vec),
                   end = pmax(beg_vec, end_vec)) %>%
    setorder(beg)
  merged <- dt[1]
  for (i in seq_len(nrow(dt))[-1]) {
    gap <- dt$beg[i] - merged$end[nrow(merged)]
    if (gap < 100) {
      merged[nrow(merged), end := max(end, dt$end[i])]
    } else {
      merged <- rbind(merged, dt[i])
    }
  }
  merged
}

########################################################
## PROCESS BLAST HITS FOR ONE READ
## Returns list: match_ranges (gene → data.table of intervals),
##               match_lens   (gene → total merged length),
##               n_plus, n_minus, first_loc, last_loc,
##               ltr_hit_dt   (LTR subject coords)

process_blast_hits <- function(blast_dt, seq) {
  blast_dt <- blast_dt[length >= MIN_MATCH_LEN]
  if (nrow(blast_dt) == 0) {
    return(list(match_ranges = list(), match_lens = list(),
                n_plus = 0L, n_minus = 0L,
                first_loc = 1e6L, last_loc = -1L,
                ltr_hit_dt = data.table(), gene_seqs = list()))
  }

  blast_dt[, gene := str_extract(sseqid, "^[^_]+")]
  blast_dt[, match_seq := {
    mb <- qstart; me <- qend
    raw <- if (mb <= me) substr(seq, mb, me) else substr(seq, me, mb)
    if (sstrand == "minus") reverse_complement(raw) else raw
  }, by = seq_len(nrow(blast_dt))]

  n_plus  <- sum(blast_dt$sstrand == "plus")
  n_minus <- sum(blast_dt$sstrand == "minus")
  first_loc <- min(c(blast_dt$qstart, blast_dt$qend))
  last_loc  <- max(c(blast_dt$qstart, blast_dt$qend))

  # Merge hits per gene
  match_ranges <- blast_dt[, {
    m <- merge_hits(qstart, qend)
    .(intervals = list(m), total_len = sum(m$end - m$beg + 1L))
  }, by = gene]

  gene_seqs <- blast_dt[, .(seqs = list(match_seq)), by = gene] %>%
    { set_names(.$seqs, .$gene) }

  ltr_hits <- blast_dt[gene == "LTR", .(sstart, send)]

  list(
    match_ranges = set_names(match_ranges$intervals, match_ranges$gene),
    match_lens   = set_names(as.list(match_ranges$total_len), match_ranges$gene),
    n_plus       = n_plus,
    n_minus      = n_minus,
    first_loc    = first_loc,
    last_loc     = last_loc,
    ltr_hit_dt   = ltr_hits,
    gene_seqs    = gene_seqs
  )
}

########################################################
## GENE ORDER CONSISTENCY CHECK
## Returns episome_flag and pruned match_ranges/match_lens.

check_gene_order <- function(match_ranges, match_lens, seq, strand, virusType) {
  gene_list <- c("LTR","GAG","POL","VIF","VPR","VPU","ENV","NEF","LTR")
  if (virusType == "SIV") gene_list[6] <- "VPX"

  get_pos <- function(gene, idx) {
    if (is.null(match_ranges[[gene]])) return(-1)
    coords <- match_ranges[[gene]]   # already a data.table with $beg / $end columns
    if (strand == "plus") coords$beg[1]
    else nchar(seq) - coords$end[1]
  }

  positions <- map_dbl(seq_along(gene_list), ~get_pos(gene_list[.x], .x))

  last_found   <- -1L
  n_out        <- 0L
  gene_out_idx <- -1L

  for (g in 2:(length(gene_list) - 1)) {
    pos <- positions[g]
    if (pos == -1) next
    if (pos < last_found) {
      if (gene_list[g] == "NEF") {
        match_ranges[["NEF"]] <- NULL
        match_lens[["NEF"]]   <- NULL
      } else {
        n_out        <- n_out + 1L
        gene_out_idx <- g
        last_found   <- pos
      }
    } else {
      if (last_found == -1 && gene_list[g] == "NEF" && positions[1] != -1) {
        match_ranges[["NEF"]] <- NULL
        match_lens[["NEF"]]   <- NULL
      } else {
        last_found <- pos
      }
    }
  }

  episome_flag <- if (n_out == 1)
    str_c("EPISOME BREAK IN ", gene_list[gene_out_idx - 1], "-", gene_list[gene_out_idx])
  else ""

  list(match_ranges = match_ranges, match_lens = match_lens,
       episome_flag = episome_flag)
}

########################################################
## ASSIGN LTR TO 5' / 3'

assign_ltr <- function(match_ranges, match_lens, min_gene, max_gene,
                       min_len, max_len, episome_flag) {
  if (is.null(match_ranges[["LTR"]])) return(list(match_ranges = match_ranges,
                                                   match_lens   = match_lens))
  intervals <- match_ranges[["LTR"]]   # data.table with $beg / $end columns
  has_multi  <- nrow(intervals) > 1

  if (has_multi) {
    if (min_gene == "LTR") { match_ranges[["5LTR"]] <- intervals[1]; match_lens[["5LTR"]] <- min_len }
    if (max_gene == "LTR") { match_ranges[["3LTR"]] <- intervals[nrow(intervals)]; match_lens[["3LTR"]] <- max_len }
    if (nchar(episome_flag) > 0) {
      ord <- order(intervals$beg)
      match_ranges[["5LTR"]] <- intervals[ord[1]]
      match_ranges[["3LTR"]] <- intervals[ord[length(ord)]]
      match_lens[["5LTR"]]   <- match_lens[["LTR"]] %||% min_len
      match_lens[["3LTR"]]   <- match_lens[["LTR"]] %||% max_len
    }
  } else {
    if (nchar(episome_flag) > 0) {
      match_ranges[["5LTR"]] <- match_ranges[["LTR"]]
      match_lens[["5LTR"]]   <- match_lens[["LTR"]]
    } else if (min_gene == "LTR") {
      match_ranges[["5LTR"]] <- match_ranges[["LTR"]]
      match_lens[["5LTR"]]   <- min_len
    } else if (max_gene == "LTR") {
      match_ranges[["3LTR"]] <- match_ranges[["LTR"]]
      match_lens[["3LTR"]]   <- max_len
    }
  }
  list(match_ranges = match_ranges, match_lens = match_lens)
}

########################################################
## GET MIN / MAX GENE POSITIONS

get_min_max_gene <- function(match_ranges, strand) {
  genes   <- names(match_ranges)
  # Put LTR last so non-LTR genes can win ties
  genes   <- c(setdiff(genes, "LTR"), intersect(genes, "LTR"))

  best <- tibble(gene = genes) %>%
    mutate(data = map(gene, function(g) {
      iv  <- match_ranges[[g]]   # data.table with $beg / $end columns
      n   <- if (g != "LTR") 1L else nrow(iv)
      iv  <- iv[seq_len(n)]
      tibble(beg = iv$beg, end = iv$end,
             len = abs(iv$end - iv$beg) + 1L, gene = g)
    })) %>%
    unnest(data)

  if (strand == "minus") {
    best <- best %>% mutate(pos = -beg) %>% arrange(pos)
  } else {
    best <- best %>% arrange(beg)
  }

  list(
    min_gene = best$gene[1],   max_gene = best$gene[nrow(best)],
    min_len  = best$len[1],    max_len  = best$len[nrow(best)]
  )
}

########################################################
## TRUNCATION FLAGS

get_trunc_flags <- function(first_loc, last_loc, seq_len, strand) {
  t5 <- 0L; t3 <- 0L
  if (first_loc < 100)            { if (strand == "plus") t5 <- 1L else t3 <- 1L }
  if (last_loc  > seq_len - 100)  { if (strand == "plus") t3 <- 1L else t5 <- 1L }
  list(t5 = t5, t3 = t3)
}

########################################################
## UPDATE MATCH STRING (replace leading/trailing 0s with -)

update_match_string <- function(s, beg_flag, end_flag) {
  if (beg_flag) {
    n   <- str_length(str_extract(s, "^0+") %||% "")
    s   <- str_replace(s, "^0+", str_c(rep("-", n), collapse = ""))
  }
  if (end_flag) {
    n   <- str_length(str_extract(s, "0+$") %||% "")
    s   <- str_replace(s, "0+$", str_c(rep("-", n), collapse = ""))
  }
  s
}

########################################################
## BUILD MATCH STRING

build_match_string <- function(match_lens, gene_len_hash, t5, t3, virusType) {
  test_genes <- c("5LTR","GAG","POL","VIF","VPR","VPU","ENV","NEF","3LTR")
  if (virusType == "SIV") test_genes[6] <- "VPX"

  vals <- map_chr(test_genes, function(gene) {
    if (is.null(match_lens[[gene]])) return("0")
    ref_gene <- if (str_detect(gene, "LTR")) "LTR" else gene
    g_len    <- gene_len_hash[[ref_gene]] %||% 1L
    m_len    <- match_lens[[gene]]
    pct      <- m_len / g_len
    if      (pct > PCT_CUTOFF)   "3"
    else if (pct > PCT_CUTOFF_2) "2"
    else                          "1"
  })
  update_match_string(str_c(vals, collapse = ""), t5, t3)
}

########################################################
## APPLY MATCH STRING POST-CORRECTIONS (port of Perl logic)

correct_match_string <- function(match_str, t5, t3) {
  # [123]--------  (LTR only, rest missing) → collapse to single gene + zeros
  if (str_detect(match_str, "[123]-{8}") && !(t5 && t3)) {
    first <- str_sub(match_str, 1, 1)
    match_str <- str_c(first, "00000000")
  }
  # --------[123]
  if (str_detect(match_str, "-{8}[123]") && !(t5 && t3)) {
    last  <- str_sub(match_str, 9, 9)
    match_str <- str_c("00000000", last)
  }
  # [123]0000000[123]  (gene at both ends, nothing in between → likely chimera)
  if (str_detect(match_str, "[123]0{7}[123]")) {
    first <- str_sub(match_str, 1, 1)
    match_str <- str_c(first, "00000000")
    both  <- t5 && t3
    match_str <- update_match_string(match_str, both, both)
  }
  match_str
}

########################################################
## CLASSIFY PROVIRAL INTEGRITY

find_match_type <- function(p) {
  case_when(
    p == "333333333"                                 ~ "INTACT",
    str_detect(p, "[0123]3{7}[0123]")               |
      p %in% c("33333333-","-33333333","-3333333-",
               "33333332-","-23333333","-2333333-",
               "33333331-","-13333333","-1333333-",
               "-1333332-","-1333331-","-2333332-",
               "-2333331-","--3333333","3333333--")  ~ "PUTATIVELY INTACT",
    str_detect(p, "^0{1,}[123]3*[123]-+$") |
      str_detect(p, "^-+[123]3*[123]0{1,}$") |
      str_detect(p, "^0{1,}3*$") | str_detect(p, "^3*0{1,}$") |
      str_detect(p, "^0{4,}") | str_detect(p, "0{4,}$") |
      str_detect(p, "^-+[12]0{2,}$")                ~ "TRUNCATED",
    str_detect(p, "^-{2,}[123]3*$") |
      str_detect(p, "^3*[123]-{2,}$") |
      str_detect(p, "^-{2,}[123]3*[123]-{2,}$") |
      str_detect(p, "^-+[12]{2}-+$") |
      str_detect(p, "^-+[123]3*[123]-+$") |
      str_detect(p, "^[12]3*[12]-+$") |
      str_detect(p, "^-+[12]3*[12]$")               ~ "INDETERMINATE",
    TRUE                                             ~ "INTERNAL DELETION"
  )
}

########################################################
## LTR POSITION COVERAGE  (vectorised update)

update_ltr_coverage <- function(ltr_pos_vec, ltr_hit_dt, ltr_len) {
  if (nrow(ltr_hit_dt) == 0) return(ltr_pos_vec)
  ltr_hit_dt[, `:=`(s = pmin(sstart, send) - 1L,
                    e = pmax(sstart, send) - 1L)]
  ltr_hit_dt <- ltr_hit_dt[(e - s + 1L) > 12]
  for (i in seq_len(nrow(ltr_hit_dt))) {
    idx <- (ltr_hit_dt$s[i] + 1L):(ltr_hit_dt$e[i] + 1L)
    idx <- idx[idx >= 1L & idx <= ltr_len]
    ltr_pos_vec[idx] <- ltr_pos_vec[idx] + 1L
  }
  ltr_pos_vec
}

########################################################
## MAKE LTR COVERAGE PLOT

make_ltr_graph <- function(ltr_len, ltr_pos_vec) {
  dt <- tibble(pos = seq_len(ltr_len) - 1L, count = ltr_pos_vec[seq_len(ltr_len)])
  write_tsv(dt, "LTRPosCnts.txt", col_names = FALSE)
  p <- ggplot(dt, aes(pos, count)) +
    geom_line() + xlab("LTR Position") + ylab("Count") +
    theme_bw()
  ggsave("LTRUsage.png", p, width = 10, height = 5, units = "in")
  system("chmod g+w LTRPosCnts.txt")
}

########################################################
## MAKE LTR LENGTH HISTOGRAMS

make_ltr_histogram <- function(ltr5_lens, ltr3_lens) {
  make_df <- function(lens_list, utr_label) {
    k <- sort(as.integer(names(lens_list)))
    tibble(UTR   = utr_label,
           BASES = k * 25L,
           COUNT = map_int(as.character(k), ~lens_list[[.x]] %||% 0L))
  }
  df <- bind_rows(make_df(ltr5_lens, "5LTR"), make_df(ltr3_lens, "3LTR"))
  write_tsv(df, "LTRHistogramData.txt")

  make_bar <- function(data, colour, fill, title) {
    ggplot(data, aes(BASES, COUNT)) +
      geom_col(color = colour, fill = fill) +
      xlab("Number of Bases") + ylab("Count") + ggtitle(title) +
      theme_bw() + theme(plot.title = element_text(hjust = 0.5))
  }

  p5   <- make_bar(filter(df, UTR == "5LTR"), "red",  "orange",     "5' LTR Lengths")
  p3   <- make_bar(filter(df, UTR == "3LTR"), "blue", "dodgerblue", "3' LTR Lengths")
  both <- p5 / p3   # patchwork stacking operator

  ggsave("LTR5_Histogram.png",    p5,   width = 10, height = 5,  units = "in")
  ggsave("LTR3_Histogram.png",    p3,   width = 10, height = 5,  units = "in")
  ggsave("LTRBOTH_Histogram.png", both, width = 10, height = 10, units = "in")
  system("chmod g+w LTRHistogramData.txt")
}

########################################################
## PRINT SUMMARY TABLE

print_summary <- function(summary_dt, virusType) {
  # summary_dt is a data.table with cols: match_str, count
  cat("\n\n")
  summary_dt %>% arrange(match_str) %>%
    { walk2(.$match_str, .$count, ~cat(.x, "\t", .y, "\n")) }

  cat("\n\n5               3\n'               '\n")
  cat("L G P V V V E N L\nT A O I P P N E T\n")
  cat(if (virusType == "HIV") "R G L F R U V F R\n" else "R G L F R X V F R\n")
  cat("=========================\n")
  summary_dt %>% arrange(match_str) %>%
    { walk2(.$match_str, .$count, ~{
        spaced <- str_c(str_split(.x, "")[[1]], collapse = " ")
        cat(spaced, "\t", .y, "\n")
      })
    }
}

########################################################
## MAIN

blast_db      <- validate_parameters(inFN, OUT_FN, reference, virusType)
blast_fa      <- str_c(blast_db, ".fa")
ref_gene_seqs <- get_reference_sequences(blast_fa)
gene_len_hash <- get_gene_lengths(reference, virusType)
ltr_len       <- gene_len_hash[["LTR"]]

message(get_sequence_count(inFN), " sequences to process")

# Read all sequences at once
fasta_tbl <- read_fasta(inFN) %>%
  mutate(seq_id = str_remove(header, "^>"))

# Accumulator structures
ltr_pos_vec <- integer(ltr_len)
num_bins    <- as.integer(ltr_len / 25) + 1L
ltr5_lens   <- setNames(as.list(integer(num_bins)), as.character(seq_len(num_bins) - 1L))
ltr3_lens   <- setNames(as.list(integer(num_bins)), as.character(seq_len(num_bins) - 1L))
gene_match_hash  <- list()   # gene → FASTA text for alignment output
result_rows      <- vector("list", nrow(fasta_tbl))
summary_counts   <- data.table(match_str = character(), count = integer())

for (i in seq_len(nrow(fasta_tbl))) {
  hdr    <- fasta_tbl$header[i]
  seq    <- fasta_tbl$seq[i]
  seq_id <- fasta_tbl$seq_id[i]
  seq_len <- nchar(seq)

  ## --- Write temp query FASTA ---
  tmp_fa <- get_tmp_path(".fa")
  write_lines(c(hdr, seq), tmp_fa)

  ## --- BLAST ---
  blast_dt <- run_blast(tmp_fa, blast_db)

  ## --- IPDA ---
  psi_val <- 0L; rre_val <- 0L
  ltr_gag_val <- 0L; pol_val <- 0L; env_val <- 0L
  if (virusType == "HIV") {
    ipda    <- find_ipda(tmp_fa)
    psi_val <- ipda$psi; rre_val <- ipda$rre
    ipda2   <- find_ipda_v2(tmp_fa)
    ltr_gag_val <- ipda2$ltr_gag; pol_val <- ipda2$pol; env_val <- ipda2$env
  }
  safe_rm(tmp_fa)

  ## --- Process hits ---
  hits <- process_blast_hits(blast_dt, seq)
  strand <- overall_strand(hits$n_minus, hits$n_plus)

  # Accumulate gene match sequences for downstream alignment
  for (g in names(hits$gene_seqs)) {
    entry <- str_c(map_chr(hits$gene_seqs[[g]], ~str_c(hdr, "\n", .x, "\n")),
                   collapse = "")
    gene_match_hash[[g]] <- str_c(gene_match_hash[[g]] %||% "", entry)
  }

  ## --- LTR coverage ---
  ltr_pos_vec <- update_ltr_coverage(ltr_pos_vec, hits$ltr_hit_dt, ltr_len)

  ## --- Gene ordering / episome check ---
  order_res    <- check_gene_order(hits$match_ranges, hits$match_lens,
                                   seq, strand, virusType)
  match_ranges <- order_res$match_ranges
  match_lens   <- order_res$match_lens
  episome_flag <- order_res$episome_flag

  ## --- Min/Max gene positions ---
  mm <- if (length(match_ranges) > 0)
    get_min_max_gene(match_ranges, strand)
  else list(min_gene="", max_gene="", min_len=0, max_len=0)

  ## --- Assign LTRs ---
  ltr_res      <- assign_ltr(match_ranges, match_lens,
                             mm$min_gene, mm$max_gene,
                             mm$min_len,  mm$max_len, episome_flag)
  match_ranges <- ltr_res$match_ranges
  match_lens   <- ltr_res$match_lens

  ## --- Truncation flags ---
  trunc <- get_trunc_flags(hits$first_loc, hits$last_loc, seq_len, strand)

  ## --- Build + correct match string ---
  match_str <- build_match_string(match_lens, gene_len_hash,
                                  trunc$t5, trunc$t3, virusType)
  match_str <- correct_match_string(match_str, trunc$t5, trunc$t3)

  ## --- LTR histogram bins ---
  ltr5_len_val <- match_lens[["5LTR"]] %||% 0
  ltr3_len_val <- match_lens[["3LTR"]] %||% 0

  if (!str_starts(match_str, "-")) {
    bin <- as.character(if (ltr5_len_val == 0) 0L else as.integer(ltr5_len_val/25) + 1L)
    ltr5_lens[[bin]] <- (ltr5_lens[[bin]] %||% 0L) + 1L
  }
  if (!str_ends(match_str, "-")) {
    bin <- as.character(if (ltr3_len_val == 0) 0L else as.integer(ltr3_len_val/25) + 1L)
    ltr3_lens[[bin]] <- (ltr3_lens[[bin]] %||% 0L) + 1L
  }

  ## --- Classify ---
  curr_type <- find_match_type(match_str)

  # Track match string counts
  existing <- summary_counts[match_str == match_str]
  if (nrow(existing) == 0) {
    summary_counts <- rbind(summary_counts,
                            data.table(match_str = match_str, count = 1L))
  } else {
    summary_counts[match_str == match_str, count := count + 1L]
  }

  ## --- Build output row ---
  row <- tibble(
    CCS_READ_ID        = seq_id,
    STRAND             = strand,
    GENE_MATCH_STRING  = match_str,
    MATCH_TYPE         = curr_type,
    EPISOME_FLAG       = episome_flag
  )
  if (virusType == "HIV") {
    row <- row %>% mutate(
      IPDA_PSI      = psi_val,
      IPDA_RRE      = rre_val,
      IPDA_INTACT   = as.integer(psi_val & rre_val),
      IPDA_LTR_GAG  = ltr_gag_val,
      IPDA_POL      = pol_val,
      IPDA_ENV      = env_val,
      IPDA_V2_INTACT = as.integer(ltr_gag_val >= 3 & pol_val >= 3 & env_val >= 3)
    )
  }
  result_rows[[i]] <- row

  if (i %% 100 == 0) message(i, " sequences processed")
}

## --- Combine results and write output ---
results_tbl <- bind_rows(result_rows)

# Reorder columns to match original output spec
col_order <- c("CCS_READ_ID","STRAND","GENE_MATCH_STRING","MATCH_TYPE")
if (virusType == "HIV") {
  col_order <- c(col_order, "IPDA_PSI","IPDA_RRE","IPDA_INTACT",
                 "IPDA_LTR_GAG","IPDA_POL","IPDA_ENV","IPDA_V2_INTACT")
}
col_order <- c(col_order, "EPISOME_FLAG")
results_tbl <- select(results_tbl, any_of(col_order))

write_tsv(results_tbl, OUT_FN, na = "")

## --- Console summary ---
type_counts <- results_tbl %>% count(MATCH_TYPE) %>% deframe()
cat("\n\n")
walk(c("INTACT","PUTATIVELY INTACT","INDETERMINATE",
       "HEAVILY TRUNCATED","TRUNCATED","INTERNAL DELETION"),
     ~cat(.x, ":", type_counts[[.x]] %||% 0L, "\n"))

## --- Gene match FASTA + (optional) alignment ---
walk(names(gene_match_hash), function(gene) {
  match_fn <- str_c(OUT_FN, "_", gene, "_matches.fa")
  ref_txt  <- ref_gene_seqs[[gene]] %||% ""
  write_lines(str_c(ref_txt, gene_match_hash[[gene]]), match_fn)
  # system(str_c("muscle -in ", shQuote(match_fn), " -out ",
  #              shQuote(str_c(OUT_FN, "_", gene, "_ALIGNMENT.fa"))))
})

## --- Plots ---
print_summary(summary_counts, virusType)
cat("MIN LTR LENGTH:", min(ltr_pos_vec[ltr_pos_vec > 0] %>% { seq_along(.) }), "\n")
make_ltr_graph(ltr_len, ltr_pos_vec)
make_ltr_histogram(ltr5_lens, ltr3_lens)
