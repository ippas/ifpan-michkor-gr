library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(ggplot2)

# =============================
# DATA
# =============================
df <- tibble::tribble(
  ~source,       ~signatureType, ~n_signif_genes, ~n_ns_genes,
  
  "DisGeNET",    "blood",        16, 3379,
  "GWASCatalog", "blood",         7, 4107,
  "genebass",    "blood",        41, 13272,
  
  "DisGeNET",    "lung",         41, 3354,
  "GWASCatalog", "lung",         52, 4062,
  "genebass",    "lung",         59, 13254,
  
  "DisGeNET",    "neural",       26, 3369,
  "GWASCatalog", "neural",       12, 4102,
  "genebass",    "neural",       50, 13263,
  
  "DisGeNET",    "systemic",     75, 3320,
  "GWASCatalog", "systemic",     10, 4104,
  "genebass",    "systemic",     62, 13251
)

signature_levels <- c("systemic", "neural", "blood", "lung")
source_levels <- c("DisGeNET", "GWASCatalog", "genebass")

df <- df %>%
  mutate(
    source = factor(source, levels = source_levels),
    signatureType = factor(signatureType, levels = signature_levels)
  )

# =============================
# COLORS
# =============================
signature_colors <- c(
  systemic = "#335C67",
  neural   = "#E09F3E",
  blood    = "#9E2A2B",
  lung     = "#540B0E"
)

# =============================
# CHI2 FUNCTION
# =============================
run_chi2_test <- function(df_subset) {
  cont_table <- as.matrix(
    df_subset %>%
      select(n_signif_genes, n_ns_genes)
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

# =============================
# POSTHOC FUNCTION
# =============================
run_posthoc <- function(df_subset) {
  x <- df_subset$n_signif_genes
  n <- df_subset$n_signif_genes + df_subset$n_ns_genes
  
  names(x) <- as.character(df_subset$source)
  names(n) <- as.character(df_subset$source)
  
  pairwise.prop.test(
    x = x,
    n = n,
    p.adjust.method = "BH"
  )
}

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
      comparison = paste(source_1, "vs", source_2),
      significance = case_when(
        p_value_adjusted < 0.001 ~ "***",
        p_value_adjusted < 0.01  ~ "**",
        p_value_adjusted < 0.05  ~ "*",
        TRUE                     ~ "n.s."
      )
    ) %>%
    select(signatureType, comparison, source_1, source_2, p_value_adjusted, significance)
}

# =============================
# OBSERVED / EXPECTED FUNCTION
# =============================
run_observed_expected_table <- function(df_subset) {
  cont_table <- as.matrix(
    df_subset %>%
      select(n_signif_genes, n_ns_genes)
  )
  
  rownames(cont_table) <- as.character(df_subset$source)
  colnames(cont_table) <- c("significant", "not_significant")
  
  chi <- chisq.test(cont_table)
  
  observed <- chi$observed
  expected <- chi$expected
  contrib  <- (observed - expected)^2 / expected
  
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

# =============================
# CHI2 RESULTS
# =============================
chi2_results <- df %>%
  group_by(signatureType) %>%
  group_modify(~ run_chi2_test(.x)) %>%
  ungroup() %>%
  mutate(
    signatureType = factor(as.character(signatureType), levels = signature_levels),
    global_label = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE            ~ "n.s."
    )
  ) %>%
  arrange(signatureType)

chi2_results

# =============================
# POSTHOC RESULTS
# =============================
df_split <- split(df, df$signatureType)

posthoc_results <- lapply(df_split, run_posthoc)

posthoc_results_df <- imap_dfr(posthoc_results, posthoc_to_df) %>%
  mutate(
    signatureType = factor(as.character(signatureType), levels = signature_levels)
  ) %>%
  arrange(signatureType, source_1, source_2)

posthoc_results_df

# =============================
# OBSERVED / EXPECTED RESULTS
# =============================
observed_expected_results <- df %>%
  group_by(signatureType) %>%
  group_modify(~ run_observed_expected_table(.x)) %>%
  ungroup() %>%
  mutate(
    signatureType = factor(as.character(signatureType), levels = signature_levels),
    source = factor(source, levels = source_levels)
  ) %>%
  arrange(signatureType, source)

observed_expected_results

# =============================
# OPTIONAL SUMMARY TABLE
# =============================
summary_results <- observed_expected_results %>%
  left_join(chi2_results, by = "signatureType") %>%
  relocate(signatureType, source, chi_square, df, p_value, global_label)

summary_results

# =============================
# DATA FOR PLOT
# =============================
plot_df <- observed_expected_results %>%
  select(signatureType, source, observed_significant, expected_significant) %>%
  pivot_longer(
    cols = c(observed_significant, expected_significant),
    names_to = "value_type",
    values_to = "count"
  ) %>%
  mutate(
    signatureType = factor(as.character(signatureType), levels = signature_levels),
    source = factor(source, levels = source_levels),
    value_type = factor(
      value_type,
      levels = c("observed_significant", "expected_significant"),
      labels = c("Observed", "Expected")
    ),
    fill_color = signature_colors[as.character(signatureType)],
    fill_color = ifelse(
      value_type == "Observed",
      scales::alpha(fill_color, 0.8),
      scales::alpha(fill_color, 0.4)
    )
  ) %>%
  arrange(signatureType, source, value_type)

# =============================
# GLOBAL LABEL POSITIONS
# =============================
sig_df <- plot_df %>%
  group_by(signatureType) %>%
  summarise(
    y = max(count) * 1.38,
    .groups = "drop"
  ) %>%
  mutate(signatureType = factor(as.character(signatureType), levels = signature_levels)) %>%
  left_join(chi2_results, by = "signatureType") %>%
  arrange(signatureType)

# =============================
# PAIRWISE TESTS FOR ALL SIGNIFICANT PANELS
# =============================
pair_tests <- list(
  c("DisGeNET", "GWASCatalog"),
  c("DisGeNET", "genebass"),
  c("GWASCatalog", "genebass")
)

significant_signature_types <- chi2_results %>%
  filter(p_value < 0.05) %>%
  pull(signatureType) %>%
  as.character()

pairwise_plot_df <- purrr::map_dfr(significant_signature_types, function(sig_type) {
  df_sig <- df %>%
    filter(signatureType == sig_type) %>%
    mutate(source = as.character(source))
  
  purrr::map_dfr(pair_tests, function(pair) {
    s1 <- pair[1]
    s2 <- pair[2]
    
    sub <- df_sig %>%
      filter(source %in% c(s1, s2))
    
    x <- sub$n_signif_genes
    n <- sub$n_signif_genes + sub$n_ns_genes
    
    test <- prop.test(x = x, n = n, correct = FALSE)
    
    tibble(
      signatureType = sig_type,
      source_1 = s1,
      source_2 = s2,
      p_value = test$p.value
    )
  })
}) %>%
  mutate(
    signatureType = factor(signatureType, levels = signature_levels),
    p_value_adj = ave(
      p_value,
      signatureType,
      FUN = function(x) p.adjust(x, method = "BH")
    ),
    label = case_when(
      p_value_adj < 0.001 ~ "***",
      p_value_adj < 0.01  ~ "**",
      p_value_adj < 0.05  ~ "*",
      TRUE                ~ "n.s."
    ),
    x1 = case_when(
      source_1 == "DisGeNET"    ~ 1,
      source_1 == "GWASCatalog" ~ 2,
      source_1 == "genebass"    ~ 3
    ),
    x2 = case_when(
      source_2 == "DisGeNET"    ~ 1,
      source_2 == "GWASCatalog" ~ 2,
      source_2 == "genebass"    ~ 3
    )
  )

# maksima osobno dla każdego panelu
panel_max_df <- plot_df %>%
  group_by(signatureType) %>%
  summarise(
    panel_max = max(count),
    .groups = "drop"
  )

pairwise_plot_df <- pairwise_plot_df %>%
  left_join(panel_max_df, by = "signatureType") %>%
  group_by(signatureType) %>%
  mutate(
    y = panel_max * c(1.08, 1.18, 1.28),
    y_label = y + panel_max * 0.03
  ) %>%
  ungroup()

pairwise_plot_df

# =============================
# PLOT
# =============================
p <- ggplot(plot_df, aes(x = source, y = count, group = value_type)) +
  geom_col(
    aes(fill = fill_color),
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  geom_text(
    data = sig_df,
    aes(x = 2, y = y, label = global_label),
    inherit.aes = FALSE,
    size = 6
  ) +
  geom_segment(
    data = pairwise_plot_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = pairwise_plot_df,
    aes(x = x1, xend = x1, y = y, yend = y - panel_max * 0.025),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = pairwise_plot_df,
    aes(x = x2, xend = x2, y = y, yend = y - panel_max * 0.025),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = pairwise_plot_df,
    aes(x = (x1 + x2) / 2, y = y_label, label = label),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  facet_wrap(~ signatureType, nrow = 1, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(
    x = "Source",
    y = "Number of significant genes"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p