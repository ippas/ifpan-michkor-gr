library(dplyr)
library(tibble)
library(purrr)

# -----------------------------
# data
# -----------------------------
df <- tibble::tribble(
  ~source,       ~signatureType, ~n_signif_phenotypes, ~n_ns_phenotypes,
  
  "DisGeNET",    "systemic",     19, 119,
  "GWASCatalog", "systemic",      1, 173,
  "genebass",    "systemic",      9, 157,
  
  "DisGeNET",    "neural",        6, 132,
  "GWASCatalog", "neural",        2, 172,
  "genebass",    "neural",        6, 160,
  
  "DisGeNET",    "blood",         5, 133,
  "GWASCatalog", "blood",         1, 173,
  "genebass",    "blood",         6, 160,
  
  "DisGeNET",    "lung",         10, 128,
  "GWASCatalog", "lung",          6, 168,
  "genebass",    "lung",          5, 161
)

# wymuszenie kolejności
df <- df %>%
  mutate(
    source = factor(source, levels = c("DisGeNET", "GWASCatalog", "genebass")),
    signatureType = factor(signatureType, levels = c("systemic", "neural", "blood", "lung"))
  )

# -----------------------------
# function for chi2 test
# -----------------------------
run_chi2_test <- function(df_subset) {
  
  cont_table <- as.matrix(
    df_subset %>%
      select(n_signif_phenotypes, n_ns_phenotypes)
  )
  
  rownames(cont_table) <- as.character(df_subset$source)
  colnames(cont_table) <- c("significant", "not_significant")
  
  chi <- chisq.test(cont_table)
  
  tibble(
    chi_square = unname(chi$statistic),
    df = unname(chi$parameter),
    p_value = chi$p.value
  )
}

# -----------------------------
# function for posthoc test
# -----------------------------
run_posthoc <- function(df_subset) {
  
  x <- df_subset$n_signif_phenotypes
  n <- df_subset$n_signif_phenotypes + df_subset$n_ns_phenotypes
  
  names(x) <- as.character(df_subset$source)
  names(n) <- as.character(df_subset$source)
  
  pairwise.prop.test(
    x = x,
    n = n,
    p.adjust.method = "BH"
  )
}

# -----------------------------
# convert posthoc result to df
# -----------------------------
posthoc_to_df <- function(posthoc_obj, signature_name) {
  
  pmat <- posthoc_obj$p.value
  
  as.data.frame(as.table(pmat), stringsAsFactors = FALSE) %>%
    filter(!is.na(Freq)) %>%
    rename(
      source_1 = Var1,
      source_2 = Var2,
      p_value_adjusted = Freq
    ) %>%
    mutate(
      signatureType = signature_name,
      comparison = paste(source_1, "vs", source_2)
    ) %>%
    select(signatureType, comparison, source_1, source_2, p_value_adjusted)
}

# -----------------------------
# observed vs expected table
# -----------------------------
run_observed_expected_table <- function(df_subset) {
  
  cont_table <- as.matrix(
    df_subset %>%
      select(n_signif_phenotypes, n_ns_phenotypes)
  )
  
  rownames(cont_table) <- as.character(df_subset$source)
  colnames(cont_table) <- c("significant", "not_significant")
  
  chi <- chisq.test(cont_table)
  
  observed <- chi$observed
  expected <- chi$expected
  
  contrib <- (observed - expected)^2 / expected
  
  tibble(
    source = rownames(observed),
    observed_significant = observed[, "significant"],
    observed_not_significant = observed[, "not_significant"],
    expected_significant = expected[, "significant"],
    expected_not_significant = expected[, "not_significant"],
    observed_significant_prop = observed[, "significant"] / rowSums(observed),
    expected_significant_prop = expected[, "significant"] / rowSums(expected),
    chi2_contrib_significant = contrib[, "significant"],
    chi2_contrib_not_significant = contrib[, "not_significant"],
    chi2_contrib_total = rowSums(contrib)
  )
}

# -----------------------------
# chi2 results
# -----------------------------
chi2_results <- df %>%
  group_by(signatureType) %>%
  group_modify(~ run_chi2_test(.x)) %>%
  ungroup()

chi2_results

# -----------------------------
# posthoc results
# -----------------------------
df_split <- split(df, df$signatureType)

posthoc_results <- lapply(df_split, run_posthoc)

posthoc_results$systemic
posthoc_results$neural
posthoc_results$blood
posthoc_results$lung

# -----------------------------
# posthoc results as one table
# -----------------------------
posthoc_results_df <- imap_dfr(posthoc_results, posthoc_to_df)

posthoc_results_df

# -----------------------------
# observed vs expected results
# -----------------------------
observed_expected_results <- df %>%
  group_by(signatureType) %>%
  group_modify(~ run_observed_expected_table(.x)) %>%
  ungroup()

observed_expected_results

# -----------------------------
# optional: merged summary
# -----------------------------
summary_results <- observed_expected_results %>%
  left_join(chi2_results, by = "signatureType") %>%
  relocate(signatureType, source, chi_square, df, p_value)

summary_results