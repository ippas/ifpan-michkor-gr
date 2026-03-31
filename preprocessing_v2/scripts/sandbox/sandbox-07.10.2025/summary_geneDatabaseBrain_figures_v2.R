library(tidyverse)
library(patchwork)

# 🔧 Przygotowanie podstawowe
summary_table <- summary_table %>%
  unique %>% 
  arrange(desc(n_records)) %>%
  mutate(
    percent_both = both_n / n_genes * 100,
    dataset_factor = factor(dataset, levels = dataset)
  )

# 🟦 1️⃣ Wykres rekordów
p_records <- summary_table %>%
  mutate(
    label_info = paste0(dataset, " (", n_publications, ", ", n_geneLists, ", ", n_geneLists_up, ", ", n_geneLists_down, ")"),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(n_records, n_records_up, n_records_down),
    names_to = "type", values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = type, group = type)) +
  geom_line(linewidth = 1) + geom_point(size = 2.8) +
  scale_color_manual(
    values = c(n_records = "black", n_records_up = "firebrick", n_records_down = "darkblue"),
    labels = c("All records", "Up-regulated", "Down-regulated")
  ) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Number of records", color = NULL, title = "A. Records") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "bottom")

# 🟩 2️⃣ Wykres genów (All / Up / Down / Both)
p_genes <- summary_table %>%
  mutate(
    label_info = paste0(dataset, " (", n_publications, ", ", n_geneLists, ", ", n_geneLists_up, ", ", n_geneLists_down, ")"),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(n_genes, only_up_n, only_down_n, both_n),
    names_to = "type", values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = type, group = type)) +
  geom_line(linewidth = 1) + geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      n_genes = "black", only_up_n = "firebrick",
      only_down_n = "darkblue", both_n = "darkgreen"
    ),
    labels = c("All genes", "Up-regulated", "Down-regulated", "Both directions")
  ) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Number of genes", color = NULL, title = "B. Genes") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "bottom")

# 🟨 3️⃣ Wykres genów unikalnych (freq = 1)
p_unique <- summary_table %>%
  mutate(
    label_info = paste0(dataset, " (", n_publications, ", ", n_geneLists, ", ", n_geneLists_up, ", ", n_geneLists_down, ")"),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(n_genes_freq1, n_genes_freq1_up, n_genes_freq1_down),
    names_to = "type", values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = type, group = type)) +
  geom_line(linewidth = 1) + geom_point(size = 2.8) +
  scale_color_manual(
    values = c(n_genes_freq1 = "gray40", n_genes_freq1_up = "darkred", n_genes_freq1_down = "darkblue"),
    labels = c("Unique genes (all)", "Up-regulated", "Down-regulated")
  ) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Genes in one list", color = NULL, title = "C. Unique genes (freq = 1)") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "bottom")

# 🟧 4️⃣ Statystyki genów na listę
p_stats <- summary_table %>%
  mutate(
    label_info = paste0(dataset, " (", n_publications, ", ", n_geneLists, ", ", n_geneLists_up, ", ", n_geneLists_down, ")"),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(min_n_genes_list, q1_n_genes_list, median_n_genes_list, mean_n_genes_list,
             q3_n_genes_list, max_n_genes_list, sd_n_genes_list),
    names_to = "statistic", values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = statistic, group = statistic)) +
  geom_line(linewidth = 1) + geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      min_n_genes_list = "gray40", q1_n_genes_list = "skyblue3",
      median_n_genes_list = "darkorange", mean_n_genes_list = "firebrick",
      q3_n_genes_list = "skyblue4", max_n_genes_list = "black",
      sd_n_genes_list = "purple4"
    ),
    labels = c("Min", "Q1", "Median", "Mean", "Q3", "Max", "SD")
  ) +
  theme_classic(base_size = 13) +
  labs(x = NULL, y = "Genes per list", color = NULL, title = "D. Distribution stats per list") +
  theme(axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        legend.position = "bottom")

# 🧩 Połączenie 4 wykresów
final_plot <- (p_records + p_genes) / (p_unique + p_stats) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

final_plot

