library(dplyr)
library(tibble)
library(purrr)
library(tidyr)
library(ggplot2)

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

# -----------------------------
# factor order
# -----------------------------
signature_levels <- c("systemic", "neural", "blood", "lung")
source_levels <- c("DisGeNET", "GWASCatalog", "genebass")

df <- df %>%
  mutate(
    source = factor(source, levels = source_levels),
    signatureType = factor(signatureType, levels = signature_levels)
  )

# -----------------------------
# colors
# -----------------------------
signature_colors <- c(
  systemic = "#335C67",
  neural   = "#E09F3E",
  blood    = "#9E2A2B",
  lung     = "#540B0E"
)

# -----------------------------
# chi2 test
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

# -----------------------------
# observed + expected
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
  
  tibble(
    source = rownames(observed),
    observed_significant = observed[, "significant"],
    expected_significant = expected[, "significant"]
  )
}

observed_expected_results <- df %>%
  group_by(signatureType) %>%
  group_modify(~ run_observed_expected_table(.x)) %>%
  ungroup() %>%
  mutate(
    source = factor(source, levels = source_levels),
    signatureType = factor(as.character(signatureType), levels = signature_levels)
  ) %>%
  arrange(signatureType, source)

# -----------------------------
# long format for plotting
# -----------------------------
plot_df <- observed_expected_results %>%
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

# -----------------------------
# y positions for global labels
# -----------------------------
sig_df <- plot_df %>%
  group_by(signatureType) %>%
  summarise(
    y = max(count) * 1.38,
    .groups = "drop"
  ) %>%
  mutate(signatureType = factor(as.character(signatureType), levels = signature_levels)) %>%
  left_join(chi2_results, by = "signatureType") %>%
  arrange(signatureType)

# -----------------------------
# systemic pairwise tests
# -----------------------------
systemic_df <- df %>%
  filter(signatureType == "systemic") %>%
  mutate(source = as.character(source))

pair_tests <- list(
  c("DisGeNET", "GWASCatalog"),
  c("DisGeNET", "genebass"),
  c("GWASCatalog", "genebass")
)

pairwise_systemic_df <- purrr::map_dfr(pair_tests, function(pair) {
  s1 <- pair[1]
  s2 <- pair[2]
  
  sub <- systemic_df %>%
    filter(source %in% c(s1, s2))
  
  x <- sub$n_signif_phenotypes
  n <- sub$n_signif_phenotypes + sub$n_ns_phenotypes
  
  test <- prop.test(x = x, n = n, correct = FALSE)
  
  tibble(
    signatureType = "systemic",
    source_1 = s1,
    source_2 = s2,
    p_value = test$p.value
  )
}) %>%
  mutate(
    signatureType = factor(signatureType, levels = signature_levels),
    p_value_adj = p.adjust(p_value, method = "BH"),
    label = case_when(
      p_value_adj < 0.001 ~ "***",
      p_value_adj < 0.01  ~ "**",
      p_value_adj < 0.05  ~ "*",
      TRUE                ~ "n.s."
    ),
    x1 = c(1, 1, 2),
    x2 = c(2, 3, 3)
  )

systemic_max <- plot_df %>%
  filter(signatureType == "systemic") %>%
  summarise(max_y = max(count)) %>%
  pull(max_y)

pairwise_systemic_df <- pairwise_systemic_df %>%
  mutate(
    y = systemic_max * c(1.08, 1.18, 1.28),
    y_label = y + systemic_max * 0.03
  )

# -----------------------------
# plot
# -----------------------------
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
    data = pairwise_systemic_df,
    aes(x = x1, xend = x2, y = y, yend = y),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = pairwise_systemic_df,
    aes(x = x1, xend = x1, y = y, yend = y - systemic_max * 0.025),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_segment(
    data = pairwise_systemic_df,
    aes(x = x2, xend = x2, y = y, yend = y - systemic_max * 0.025),
    inherit.aes = FALSE,
    linewidth = 0.5
  ) +
  geom_text(
    data = pairwise_systemic_df,
    aes(x = (x1 + x2) / 2, y = y_label, label = label),
    inherit.aes = FALSE,
    size = 4.5
  ) +
  facet_wrap(~ signatureType, nrow = 1, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.24))) +
  labs(
    x = "Source",
    y = "Number of significant phenotypes"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

p