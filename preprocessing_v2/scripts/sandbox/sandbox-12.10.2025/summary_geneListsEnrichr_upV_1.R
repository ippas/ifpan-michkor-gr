library(dplyr)
library(stringr)
library(ggplot2)

df_all <- gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022 %>%
  filter(grepl("NR3C1|NR3C2", Term, ignore.case = TRUE)) %>%
  mutate(
    receptor = case_when(
      str_detect(Term, "NR3C1") ~ "NR3C1 (Glucocorticoid receptor)",
      str_detect(Term, "NR3C2") ~ "NR3C2 (Mineralocorticoid receptor)"
    ),
    log10FDR = -log10(Adjusted.P.value)
  )

df_best <- df_all %>%
  group_by(receptor) %>%
  slice_min(Adjusted.P.value, n = 1, with_ties = FALSE) %>%
  mutate(
    best_log10FDR = log10FDR,
    best_n_genes = as.numeric(sub("/.*", "", Overlap))
  ) %>%
  select(receptor, best_log10FDR, best_n_genes)

df_counts <- df_all %>% count(receptor, name = "n_terms")

df_plot <- df_best %>% left_join(df_counts, by = "receptor")

ggplot(df_plot, aes(x = receptor, y = best_log10FDR)) +
  geom_col(fill = "gray80", color = "black", width = 0.6) +
  geom_point(
    data = df_all,
    aes(x = receptor, y = log10FDR),
    color = "black", size = 2, position = position_jitter(width = 0.08)
  ) +
  geom_text(aes(label = paste0(best_n_genes, " genes\n(", n_terms, " terms)")),
            vjust = -0.3, size = 4) +
  labs(
    x = NULL,
    y = expression(-log[10]("FDR")),
    title = "Top ChIP-Seq enrichment per receptor (ChEA_2022)"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))


gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs |>
  dplyr::filter(Adjusted.P.value < 0.05) |>
  dplyr::filter(grepl("dexamethasone |budesonide |methylprednisolone |hydrocortisone ", Term, ignore.case = TRUE)) |>
  (\(df) {
    total_genes <- df |>
      dplyr::pull(Genes) |>
      stringr::str_split(";") |>
      unlist() |>
      unique() |>
      length()
    
    df |>
      dplyr::mutate(
        log10FDR = -log10(Adjusted.P.value),
        Term = factor(Term, levels = rev(unique(Term)))
      ) |>
      ggplot(aes(x = log10FDR, y = Term)) +
      geom_bar(stat = "identity", fill = "grey80", color = "black", width = 0.6) +
      geom_text(aes(label = paste0("n=", n_genes)), hjust = -0.1, size = 3.5) +
      labs(
        x = expression(-log[10]("FDR")),
        y = NULL,
        title = paste0(
          "Glucocorticoid-related compounds (LINCS_L1000_Chem_Pert_Consensus_Sigs)\n",
          "Total unique genes: ", total_genes
        )
      ) +
      theme_classic(base_size = 12) +
      theme(axis.text.y = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5))
  })()


df_all <- gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>%
  filter(grepl("astrocyte|neuron|microglial|oligodendrocyte", Term, ignore.case = TRUE)) %>%
  filter(grepl("brain", Term, ignore.case = TRUE)) %>%
  mutate(
    cell_type = case_when(
      str_detect(Term, "astrocyte") ~ "Astrocyte",
      str_detect(Term, "neuron") ~ "Neuron",
      str_detect(Term, "microglial") ~ "Microglia",
      str_detect(Term, "oligodendrocyte precursor") ~ "OPC (Oligodendrocyte precursor)",
      str_detect(Term, "oligodendrocyte") ~ "Oligodendrocyte",
      TRUE ~ "Other"
    ),
    log10FDR = -log10(Adjusted.P.value)
  )

df_best <- df_all %>%
  group_by(cell_type) %>%
  slice_min(Adjusted.P.value, n = 1, with_ties = FALSE) %>%
  mutate(
    best_log10FDR = log10FDR,
    best_n_genes = as.numeric(sub("/.*", "", Overlap))
  ) %>%
  select(cell_type, best_log10FDR, best_n_genes)

df_counts <- df_all %>% count(cell_type, name = "n_terms")
df_plot <- df_best %>% left_join(df_counts, by = "cell_type")

ggplot(df_plot, aes(x = reorder(cell_type, best_log10FDR), y = best_log10FDR)) +
  geom_col(fill = "grey85", color = "black", width = 0.6) +
  geom_point(
    data = df_all,
    aes(x = cell_type, y = log10FDR),
    color = "black", fill = "black", size = 2,
    position = position_jitter(width = 0.1)
  ) +
  geom_text(aes(label = paste0(best_n_genes, " genes\n(", n_terms, " terms)")),
            vjust = -0.3, size = 3.5) +
  labs(
    x = NULL,
    y = expression(-log[10]("FDR")),
    title = "GR-upregulated genes enriched in brain cell types (CellMarker_2024)"
  ) +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))

