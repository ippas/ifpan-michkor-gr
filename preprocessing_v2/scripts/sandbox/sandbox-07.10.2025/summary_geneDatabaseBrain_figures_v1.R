summary_table %>% 
  arrange(desc(n_records)) %>% 
  as.data.frame() %>% 
  mutate(percent_both = both_n/n_genes *100) %>% 
  mutate(dataset_factor = factor(dataset, levels = dataset)) -> summary_table




summary_table %>%
  arrange(desc(n_records)) %>%
  mutate(
    label_info = paste0(
      dataset, " (",
      n_publications, ", ",
      n_geneLists, ", ",
      n_geneLists_up, ", ",
      n_geneLists_down, ")"
    ),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(n_records, n_records_up, n_records_down),
    names_to = "type",
    values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = type, group = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      n_records = "black",
      n_records_up = "firebrick",
      n_records_down = "darkblue"
    ),
    labels = c(
      n_records = "All records",
      n_records_up = "Up-regulated",
      n_records_down = "Down-regulated"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Dataset (n_publications, n_geneLists, n_geneLists_up, n_geneLists_down)",
    y = "Number of records",
    color = "Category"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "bottom"
  )


# 1️⃣ główny wykres
p_genes_main <- summary_table %>%
  arrange(desc(n_genes)) %>%
  mutate(
    label_info = paste0(
      dataset, " (",
      n_publications, ", ",
      n_geneLists, ", ",
      n_geneLists_up, ", ",
      n_geneLists_down, ")"
    ),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(n_genes, only_up_n, only_down_n, both_n),
    names_to = "type",
    values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = type, group = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      n_genes = "black",
      only_up_n = "firebrick",
      only_down_n = "darkblue",
      both_n = "darkgreen"
    ),
    labels = c(
      n_genes = "All genes",
      only_up_n = "Up-regulated",
      only_down_n = "Down-regulated",
      both_n = "Both directions"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Dataset (n_publications, n_geneLists, n_geneLists_up, n_geneLists_down)",
    y = "Number of genes",
    color = "Category"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# 2️⃣ wykres dla genów występujących tylko w jednej liście
p_genes_unique <- summary_table %>%
  arrange(desc(n_genes)) %>%
  mutate(
    label_info = paste0(
      dataset, " (",
      n_publications, ", ",
      n_geneLists, ", ",
      n_geneLists_up, ", ",
      n_geneLists_down, ")"
    ),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(n_genes_freq1, n_genes_freq1_up, n_genes_freq1_down),
    names_to = "type",
    values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = type, group = type)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      n_genes_freq1 = "gray30",
      n_genes_freq1_up = "darkred",
      n_genes_freq1_down = "darkblue"
    ),
    labels = c(
      n_genes_freq1 = "Genes occurring in one list",
      n_genes_freq1_up = "Up-regulated (one list)",
      n_genes_freq1_down = "Down-regulated (one list)"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Dataset (n_publications, n_geneLists, n_geneLists_up, n_geneLists_down)",
    y = "Number of genes occurring in one list",
    color = "Category"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "bottom"
  )

# 🔄 wspólna skala osi Y
ymax <- max(
  max(summary_table$n_genes, na.rm = TRUE),
  max(summary_table$n_genes_freq1, na.rm = TRUE)
)
p_genes_main <- p_genes_main + coord_cartesian(ylim = c(0, ymax))
p_genes_unique <- p_genes_unique + coord_cartesian(ylim = c(0, ymax))

# 🔗 połączenie z patchwork
p_genes_main + p_genes_unique +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.title = element_blank())


p_genes_main + p_genes_freq1


summary_table %>%
  arrange(desc(median_n_genes_list)) %>%
  mutate(
    label_info = paste0(
      dataset, " (",
      n_publications, ", ",
      n_geneLists, ", ",
      n_geneLists_up, ", ",
      n_geneLists_down, ")"
    ),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(
      min_n_genes_list,
      q1_n_genes_list,
      median_n_genes_list,
      mean_n_genes_list,
      q3_n_genes_list,
      max_n_genes_list,
      sd_n_genes_list
    ),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = statistic, group = statistic)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      min_n_genes_list = "gray40",
      q1_n_genes_list = "skyblue3",
      median_n_genes_list = "darkorange",
      mean_n_genes_list = "firebrick",
      q3_n_genes_list = "skyblue4",
      max_n_genes_list = "black",
      sd_n_genes_list = "purple4"
    ),
    labels = c(
      min_n_genes_list = "Min",
      q1_n_genes_list = "Q1",
      median_n_genes_list = "Median",
      mean_n_genes_list = "Mean",
      q3_n_genes_list = "Q3",
      max_n_genes_list = "Max",
      sd_n_genes_list = "SD"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Dataset (n_publications, n_geneLists, n_geneLists_up, n_geneLists_down)",
    y = "Number of genes per list",
    color = "Statistic"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "bottom"
  )


summary_table %>%
  arrange(desc(n_genes)) %>%
  mutate(
    label_info = paste0(
      dataset, " (",
      n_publications, ", ",
      n_geneLists, ", ",
      n_geneLists_up, ", ",
      n_geneLists_down, ")"
    ),
    label_info = factor(label_info, levels = label_info)
  ) %>%
  pivot_longer(
    cols = c(
      min_n_genes_list,
      q1_n_genes_list,
      median_n_genes_list,
      mean_n_genes_list,
      q3_n_genes_list,
      max_n_genes_list,
      sd_n_genes_list
    ),
    names_to = "statistic",
    values_to = "value"
  ) %>%
  ggplot(aes(x = label_info, y = value, color = statistic, group = statistic)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.8) +
  scale_color_manual(
    values = c(
      min_n_genes_list = "gray40",
      q1_n_genes_list = "skyblue3",
      median_n_genes_list = "darkorange",
      mean_n_genes_list = "firebrick",
      q3_n_genes_list = "skyblue4",
      max_n_genes_list = "black",
      sd_n_genes_list = "purple4"
    ),
    labels = c(
      min_n_genes_list = "Min",
      q1_n_genes_list = "Q1",
      median_n_genes_list = "Median",
      mean_n_genes_list = "Mean",
      q3_n_genes_list = "Q3",
      max_n_genes_list = "Max",
      sd_n_genes_list = "SD"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(
    x = "Dataset (n_publications, n_geneLists, n_geneLists_up, n_geneLists_down)",
    y = "Number of genes per list",
    color = "Statistic"
  ) +
  theme(
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1),
    legend.position = "bottom"
  )

