#!/usr/bin/env Rscript

###############################################################################
# Permutation-based empirical FDR (HIT-COUNT LEVEL)
# Randomize: 2 GR gene lists (columns) ONLY
# Fixed: phenotype gene lists (rows)
#
# Output:
#  - observed_hits (from real GR lists)
#  - random_hits distribution (n_iter permutations; each uses 2 new random GR lists)
#  - empirical_p = P(random_hits >= observed_hits)
#  - FDR_hat = E[random_hits] / observed_hits (cap at 1)
###############################################################################

message("============================================================")
message("🚀 START: Permutation FDR (hit-count) | RANDOMIZE GR LISTS")
message("   • Randomized: 2 GR gene lists (length-matched)")
message("   • Fixed     : phenotype gene lists")
message("============================================================")

###############################################################################
# [1/9] LOAD LIBRARIES
###############################################################################
message("[1/9] 📦 Loading libraries")

suppressPackageStartupMessages({
  library(tidyverse)
  library(glue)
})

message("   ✔ Packages loaded")

###############################################################################
# [2/9] SOURCE REQUIRED FUNCTIONS
###############################################################################
message("[2/9] 📦 Sourcing required functions")

source("preprocessing_v2/src/randomGeneLists_utils/generate_random_gene_lists_by_lengths.R")

source("preprocessing_v2/src/overlap_utils/perform_chi2_tests.R")
source("preprocessing_v2/src/overlap_utils/perform_chi2_tests.R")
source("preprocessing_v2/src/overlap_utils/processing_overlap_results.R")
source("preprocessing_v2/src/overlap_utils/visualization-utils/heatmap_overlap_log2OR_ggplot.R")
source("preprocessing_v2/src/overlap_utils/visualization-utils/heatmap_overlap_logCHI2_ggplot.R")
source("preprocessing_v2/src/overlap_utils/run_full_overlap_analysis.R")

message("   ✔ All functions sourced")

###############################################################################
# [3/9] LOAD INPUT RDATA
###############################################################################
message("[3/9] 📦 Loading RData input package")

# Ten plik powinien zawierać:
#  - flat_allGrSignatures_31.10.2025
#  - hgnc_symbols_vector_v110
#  - sig_names (2 elementy!)
#  - gene_list (named list phenotype -> gengrepl("(^|_)(G3x|F[0-9]x)(_|$)", names(gene_list)) &
!grepl("F[0-9]x_.*F[0-9]x", names(gene_list))
es)
rdata_path <- "data/run_RData_inputs/permutation_fdr/rdata2permutationFDR_systemicGR_allSources.RData"
load(rdata_path)

message("   ✔ RData loaded: ", rdata_path)

stopifnot(exists("flat_allGrSignatures_31.10.2025"))
stopifnot(exists("hgnc_symbols_vector_v110"))
stopifnot(exists("sig_names"))
stopifnot(exists("gene_list"))

stopifnot(length(sig_names) == 2)
stopifnot(all(sig_names %in% names(flat_allGrSignatures_31.10.2025)))

message("   • background genes      : ", length(hgnc_symbols_vector_v110))
message("   • GR signatures (n)     : ", length(sig_names), "  (must be 2)")
message("   • phenotype lists (n)   : ", length(gene_list))

###############################################################################
# [4/9] PARAMETERS + HIT CRITERIA
###############################################################################
message("[4/9] ⚙ Defining parameters & hit criteria")

# permutations
n_iter <- 2L
seed   <- 123L

# hit definition (identyczna jak u Ciebie)
p_thr <- 0.05
min_overlap_count <- 3L         # gene_overlap_count > 2  <=> >=3
require_or_positive <- TRUE     # odds_ratio > 0

# outputs
out_dir <- "data/results/permutation_fdr/hitcount_randomGR_2lists/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_dir, glue("permFDR_hitcount_randomGR2_n{n_iter}.rds"))
out_tsv <- file.path(out_dir, glue("permFDR_hitcount_randomGR2_n{n_iter}.tsv"))

message("   • n_iter              : ", n_iter)
message("   • seed                : ", seed)
message("   • hit criteria        : overlap>2, OR>0, p<0.05")
message("   • out_dir             : ", out_dir)

set.seed(seed)

###############################################################################
# [5/9] PREP: GR lengths (these are what we randomize)
###############################################################################
message("[5/9] 📏 Preparing GR list lengths (to generate 2 random GR lists/iter)")

gr_lengths <- vapply(
  flat_allGrSignatures_31.10.2025[sig_names],
  length,
  integer(1)
)

message("   ✔ GR lengths ready")
message("   • ", sig_names[1], " length: ", gr_lengths[[1]])
message("   • ", sig_names[2], " length: ", gr_lengths[[2]])

message("   • Random GR lists per iteration : ", length(sig_names))
message("   • Total random GR lists         : ", length(sig_names) * n_iter)

###############################################################################
# [6/9] OBSERVED RUN → observed_hits (real GR vs fixed phenotypes)
###############################################################################
message("[6/9] 🧾 Running observed overlap analysis (real GR lists)")

obs <- run_full_overlap_analysis(
  gene_lists = c(
    flat_allGrSignatures_31.10.2025[sig_names],  # REAL GR
    gene_list                                 # FIXED phenotypes
  ),
  total_genes = hgnc_symbols_vector_v110,
  cols_to_filter = sig_names,
  rows_to_filter = names(gene_list),
  plot_title_or = "",
  triangle_mode = "full",
  fdr_threshold = 1,
  data_type = "original_data",
  verbose = FALSE,
  palette_or = c("#c6d3e3", "white", "darkred"),
  text_contrast_range_or = c(-30, 4.9)
)

obs_df <- obs$processed$original_data$df %>%
  mutate(row = Var1, col = Var2) %>%
  filter(row %in% names(gene_list), col %in% sig_names)

# observed hits with SAME criteria
obs_hits_df <- obs_df %>%
  filter(gene_overlap_count >= min_overlap_count) %>%
  { if (require_or_positive) filter(., odds_ratio > 0) else . } %>%
  filter(p_value < p_thr)

observed_hits <- nrow(obs_hits_df)

message("   ✔ Observed hits: ", observed_hits)
message("   • (criteria: overlap>2, OR>0, p<0.05)")

###############################################################################
# [7/9] PERMUTATION LOOP: randomize GR lists only
###############################################################################
message("[7/9] 🔁 Starting permutation loop (randomize GR only)")

random_hits <- integer(n_iter)

t_start <- Sys.time()
iter_times <- numeric(n_iter)

log_every <- 25L

for (i in seq_len(n_iter)) {
  
  t0 <- Sys.time()
  
  # 7.1) Generate 2 RANDOM GR lists (length-matched to real sigs)
  random_gr <- generate_random_gene_lists_by_lengths(
    gene_list = hgnc_symbols_vector_v110,
    lengths   = gr_lengths
  )
  names(random_gr) <- sig_names
  
  # 7.2) Run overlap analysis: RANDOM GR vs FIXED phenotypes
  perm <- run_full_overlap_analysis(
    gene_lists = c(
      random_gr,   # RANDOM GR (2 lists)
      gene_list    # FIXED phenotype lists
    ),
    total_genes = hgnc_symbols_vector_v110,
    cols_to_filter = sig_names,
    rows_to_filter = names(gene_list),
    plot_title_or = "",
    triangle_mode = "full",
    fdr_threshold = 1,
    data_type = "original_data",
    verbose = FALSE,
    palette_or = c("#c6d3e3", "white", "darkred"),
    text_contrast_range_or = c(-30, 4.9)
  )
  
  perm_df <- perm$processed$original_data$df %>%
    mutate(row = Var1, col = Var2) %>%
    filter(row %in% names(gene_list), col %in% sig_names)
  
  # 7.3) Count hits in this permutation (same criteria)
  perm_hits <- perm_df %>%
    filter(gene_overlap_count >= min_overlap_count) %>%
    { if (require_or_positive) filter(., odds_ratio > 0) else . } %>%
    filter(p_value < p_thr) %>%
    nrow()
  
  random_hits[i] <- perm_hits
  iter_times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  
  if (i %% log_every == 0L || i == 1L || i == n_iter) {
    
    elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    avg_sec <- mean(iter_times[seq_len(i)])
    eta_sec <- avg_sec * (n_iter - i)
    
    message(sprintf(
      "⏳ %d/%d | iter=%.2fs | avg=%.2fs | elapsed=%s | ETA=%s | hits=%d",
      i, n_iter,
      iter_times[i],
      avg_sec,
      format(as.POSIXct(elapsed_sec, origin="1970-01-01", tz="UTC"), "%H:%M:%S"),
      format(as.POSIXct(eta_sec, origin="1970-01-01", tz="UTC"), "%H:%M:%S"),
      random_hits[i]
    ))
  }
}

###############################################################################
# [8/9] SUMMARY: empirical p + FDR_hat
###############################################################################
message("[8/9] 📊 Computing empirical p-value and FDR estimate")

empirical_p <- mean(random_hits >= observed_hits)
expected_fp <- mean(random_hits)
FDR_hat <- if (observed_hits == 0) NA_real_ else min(1, expected_fp / observed_hits)

q95 <- as.numeric(quantile(random_hits, 0.95))
q99 <- as.numeric(quantile(random_hits, 0.99))

message("============================================================")
message("✅ PERMUTATION (RANDOM GR) COMPLETED")
message("------------------------------------------------------------")
message("• observed_hits                  : ", observed_hits)
message("• E[random_hits] (expected FP)   : ", round(expected_fp, 3))
message("• empirical_p (>= observed_hits) : ", signif(empirical_p, 4))
message("• FDR_hat = E[random]/observed   : ", ifelse(is.na(FDR_hat), "NA", signif(FDR_hat, 4)))
message("• random_hits 95th percentile    : ", q95)
message("• random_hits 99th percentile    : ", q99)
message("============================================================")

###############################################################################
# [9/9] SAVE OUTPUTS
###############################################################################
message("[9/9] 💾 Saving outputs")

result <- list(
  params = list(
    n_iter = n_iter,
    seed = seed,
    p_thr = p_thr,
    min_overlap_count = min_overlap_count,
    require_or_positive = require_or_positive,
    randomized = "GR lists only",
    randomized_lists_per_iter = length(sig_names),
    total_randomized_lists = length(sig_names) * n_iter,
    gr_lengths = gr_lengths
  ),
  observed_hits = observed_hits,
  random_hits = random_hits,
  empirical_p = empirical_p,
  expected_fp = expected_fp,
  FDR_hat = FDR_hat,
  random_hits_q95 = q95,
  random_hits_q99 = q99
)

saveRDS(result, out_rds)

summary_df <- tibble(
  observed_hits = observed_hits,
  expected_fp = expected_fp,
  empirical_p = empirical_p,
  FDR_hat = FDR_hat,
  random_hits_q95 = q95,
  random_hits_q99 = q99,
  n_iter = n_iter,
  randomized_lists_per_iter = length(sig_names),
  total_randomized_lists = length(sig_names) * n_iter,
  gr_length_1 = unname(gr_lengths[[1]]),
  gr_length_2 = unname(gr_lengths[[2]])
)

readr::write_tsv(summary_df, out_tsv)

message("💾 Saved:")
message("   • RDS: ", out_rds)
message("   • TSV: ", out_tsv)
message("🏁 DONE")