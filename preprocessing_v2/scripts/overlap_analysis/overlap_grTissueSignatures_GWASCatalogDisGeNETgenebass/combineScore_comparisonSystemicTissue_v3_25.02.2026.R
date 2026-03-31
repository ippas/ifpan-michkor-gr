library(dplyr)
library(tidyr)
library(ggplot2)

panels <- c("systemic", "neural", "blood", "lung")

# ---- 1. Significant labels per panel ----
sig_labels_by_panel <- grSystemicTissues_all %>%
  mutate(
    label = paste0(phenotype, "_", source),
    panel = sub("(Up|Down)$", "", grSignature)
  ) %>%
  filter(panel %in% panels,
         p_value < 0.05,
         observed_overlap >= 3) %>%
  distinct(panel, label)

# ---- 2. Data for plotting ----
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
    direction = factor(direction, levels = c("Up", "Down"))
  )

# ---- 3. Wilcoxon paired per panel ----
stats_df <- plot_df %>%
  select(panel, label, direction, combine_score) %>%
  pivot_wider(names_from = direction, values_from = combine_score) %>%
  drop_na() %>%
  group_by(panel) %>%
  summarise(
    p_value = wilcox.test(Up, Down, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label_txt = paste0("p = ",
                            format.pval(p_value, digits = 2, eps = 1e-300)))

# ---- 4. Global (frozen) y positions for bracket + text ----
y_min_global <- min(plot_df$combine_score, na.rm = TRUE)
y_max_global <- max(plot_df$combine_score, na.rm = TRUE)
y_range_global <- y_max_global - y_min_global

stats_df <- stats_df %>%
  mutate(
    y_bracket = y_max_global + 0.06 * y_range_global,
    y_text    = y_max_global + 0.11 * y_range_global,
    y_tick    = y_max_global + 0.03 * y_range_global
  )

# ---- 5. Plot ----
p <- ggplot(plot_df, aes(x = direction, y = combine_score)) +
  geom_boxplot(alpha = 0.4, width = 0.5, outlier.shape = NA) +
  geom_jitter(aes(color = direction), width = 0.12, size = 2, alpha = 0.7) +
  scale_color_manual(values = c("Up" = "firebrick3",
                                "Down" = "navy")) +
  facet_wrap(~ panel, ncol = 4) +   # <-- 4 kolumny, wspólna oś Y
  # bracket (pozioma linia)
  geom_segment(data = stats_df,
               aes(x = 1, xend = 2, y = y_bracket, yend = y_bracket),
               inherit.aes = FALSE, linewidth = 0.6) +
  # "zawiasy"
  geom_segment(data = stats_df,
               aes(x = 1, xend = 1, y = y_bracket, yend = y_tick),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_segment(data = stats_df,
               aes(x = 2, xend = 2, y = y_bracket, yend = y_tick),
               inherit.aes = FALSE, linewidth = 0.6) +
  # p-value
  geom_text(data = stats_df,
            aes(x = 1.5, y = y_text, label = label_txt),
            inherit.aes = FALSE,
            size = 3.5) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(x = NULL, y = "combine_score") +
  coord_cartesian(ylim = c(y_min_global, y_max_global + 0.16 * y_range_global))


svg("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/wilcoxonPaired_allGr_UpVsDown.svg",
    width = 6,
    height = 3)
p
dev.off()
