###############################################################################
# Permutation-based FDR estimation
# GWAS Catalog – Mental Health
###############################################################################

message("============================================================")
message("🚀 START: Permutation FDR analysis")
message("   • Dataset : GWAS Catalog – Mental Health")
message("   • Method  : permutation-based empirical FDR")
message("============================================================")


###############################################################################
# [1/7] LOAD INPUT DATA
###############################################################################
message("[1/7] 📦 Loading input data")
message("   • Loading RData object")

load("data/run_RData_inputs/permutation_fdr/rdata2permuationFDR_GWASCatalogMentalHealth.RData")

message("   ✔ RData loaded successfully")


###############################################################################
# [2/7] SOURCE REQUIRED FUNCTIONS
###############################################################################
message("[2/7] 📦 Sourcing required functions")
message("   • Random gene list generation")
message("   • Overlap & chi² utilities")
message("   • Visualization helpers")
message("   • Full overlap pipeline")

source("preprocessing_v2/src/randomGeneLists_utils/generate_random_gene_lists_by_lengths.R")

source("preprocessing_v2/src/overlap_utils/perform_chi2_tests.R")
source("preprocessing_v2/src/overlap_utils/perform_chi2_tests.R")
source("preprocessing_v2/src/overlap_utils/processing_overlap_results.R")
source("preprocessing_v2/src/overlap_utils/visualization-utils/heatmap_overlap_log2OR_ggplot.R")
source("preprocessing_v2/src/overlap_utils/visualization-utils/heatmap_overlap_logCHI2_ggplot.R")
source("preprocessing_v2/src/overlap_utils/run_full_overlap_analysis.R")

message("   ✔ All functions sourced")


###############################################################################
# [3/7] LOAD LIBRARIES
###############################################################################
message("[3/7] 📦 Loading standard analysis packages")

message("   • tidyverse (dplyr, tidyr, ggplot2, readr, purrr, tibble)")
message("   • stringr / forcats")
message("   • scales / patchwork")
message("   • data.table / glue")

suppressPackageStartupMessages({
  
  library(tidyverse)    # dplyr, tidyr, ggplot2, readr, purrr, tibble
  library(stringr)      # string operations
  library(forcats)      # factor handling
  library(scales)       # scales for ggplot
  library(patchwork)    # combining ggplots
  library(data.table)   # fast data manipulation
  library(glue)         # clean messages / file names
  
})

message("   ✔ All packages loaded successfully")

###############################################################################
# [4/7] DEFINE PARAMETERS
###############################################################################
message("[4/7] ⚙ Defining permutation parameters")

lengths_vec <- c(124, 208, 147, 76, 132, 86, 208, 119)

n_iter <- 1000L
target <- 31L

random_hits <- integer(n_iter)

message("   • Number of permutations      : ", n_iter)
message("   • Target (observed hits)      : ", target)
message("   • Random list sizes (n = 8)   : ", paste(lengths_vec, collapse = ", "))


###############################################################################
# [5/7] MAIN PERMUTATION LOOP
###############################################################################
message("[5/7] 🔁 Starting permutation loop")

# ---- timing helpers ----
t_start <- Sys.time()
iter_times <- numeric(n_iter)

# jak często logować (np. co 1 iterację = zawsze; co 10 = rzadziej)
log_every <- 1L

for (i in seq_len(n_iter)) {
  
  t_iter_start <- Sys.time()
  
  # 5.1) Generate random gene lists
  random_8geneLists <- generate_random_gene_lists_by_lengths(
    gene_list = hgnc_symbols_vector_v110,
    lengths   = lengths_vec
  )
  
  # 5.2) Run overlap analysis
  res <- run_full_overlap_analysis(
    gene_lists = c(
      flat_allGrSignatures_31.10.2025[sig_names],
      random_8geneLists,
      gene_list
    ),
    total_genes = hgnc_symbols_vector_v110,
    cols_to_filter = c(sig_names, names(random_8geneLists)),
    rows_to_filter = names(gene_list),
    plot_title_or = "",
    triangle_mode = "full",
    fdr_threshold = 1,
    data_type = "original_data",
    verbose = FALSE,
    palette_or = c("#c6d3e3", "white", "darkred"),
    text_contrast_range_or = c(-30, 4.9)
  )
  
  # 5.3) Count significant random hits
  random_hits[i] <- res$processed$original_data$df %>%
    filter(gene_overlap_count > 2) %>%
    filter(odds_ratio > 0) %>%
    filter(p_value < 0.05) %>%
    filter(Var2 %in% names(random_8geneLists)) %>%
    nrow()
  
  # ---- timing stats ----
  iter_times[i] <- as.numeric(difftime(Sys.time(), t_iter_start, units = "secs"))
  
  # ---- progress log (minimal) ----
  if (i %% log_every == 0L || i == 1L || i == n_iter) {
    
    elapsed_sec <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
    avg_sec <- mean(iter_times[seq_len(i)])
    remaining <- n_iter - i
    eta_sec <- avg_sec * remaining
    
    message(sprintf(
      "⏳ %d/%d | iter=%.1fs | avg=%.1fs | elapsed=%s | ETA=%s | hits=%d",
      i, n_iter,
      iter_times[i],
      avg_sec,
      format(as.POSIXct(elapsed_sec, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S"),
      format(as.POSIXct(eta_sec, origin = "1970-01-01", tz = "UTC"), "%H:%M:%S"),
      random_hits[i]
    ))
  }
}

###############################################################################
# [6/7] SUMMARY STATISTICS
###############################################################################
message("[6/7] 📊 Computing empirical statistics")

n_better_or_equal <- sum(random_hits >= target)
empirical_pvalue  <- mean(random_hits >= target)

message("   • Permutations ≥ target : ", n_better_or_equal)
message("   • Empirical p-value     : ", round(empirical_pvalue, 5))

message("random hits:")
print("#######################################################################")
print(random_hits)
print("#######################################################################")
###############################################################################
# [7/7] FINAL SUMMARY
###############################################################################
message("============================================================")
message("✅ PERMUTATION ANALYSIS COMPLETED")
message("------------------------------------------------------------")
message("• Observed target hits        : ", target)
message("• Permutations ≥ target hits : ", n_better_or_equal)
message("• Empirical p-value           : ", round(empirical_pvalue, 5))
message("============================================================")
message("🏁 DONE")
