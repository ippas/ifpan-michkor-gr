library(dplyr)
library(tidyr)
library(ggplot2)

panels <- c("systemic", "neural", "blood", "lung")

sig_labels_by_panel <- grSystemicTissues_all %>%
  mutate(
    label = paste0(phenotype, "_", source),
    panel = sub("(Up|Down)$", "", grSignature)
  ) %>%
  filter(panel %in% panels) %>%
  filter(p_value < 0.05, observed_overlap >= 3) %>%
  distinct(panel, label)

plot_df <- grSystemicTissues_all %>%
  mutate(
    label = paste0(phenotype, "_", source),
    panel = sub("(Up|Down)$", "", grSignature),
    direction = ifelse(grepl("Up$", grSignature), "Up", "Down")
  ) %>%
  filter(panel %in% panels) %>%
  inner_join(sig_labels_by_panel, by = c("panel", "label")) %>%
  mutate(
    panel = factor(panel, levels = panels),
    direction = factor(direction, levels = c("Up", "Down"))  # <-- tutaj zmiana
  )

p <- ggplot(plot_df, aes(x = direction, y = combine_score)) +
  geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = direction),
              width = 0.12,
              size = 2,
              alpha = 0.7) +
  scale_color_manual(values = c("Up" = "firebrick3",
                                "Down" = "navy")) +
  facet_wrap(~ panel, ncol = 4) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "combine_score")

p
