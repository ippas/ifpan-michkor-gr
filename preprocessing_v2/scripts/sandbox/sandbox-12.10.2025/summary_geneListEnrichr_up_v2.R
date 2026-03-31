library(dplyr)
library(ggplot2)
library(stringr)

# --- 1️⃣ ChEA (NR3C1 / NR3C2) ---
df_chea <- gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022 %>%
  filter(grepl("NR3C1|NR3C2", Term, ignore.case = TRUE)) %>%
  mutate(
    category = case_when(
      str_detect(Term, "NR3C1") ~ "NR3C1",
      str_detect(Term, "NR3C2") ~ "NR3C2",
      TRUE ~ Term
    ),
    source   = "ChEA_2022",
    log10FDR = -log10(Adjusted.P.value)
  ) %>%
  select(source, category, log10FDR, Overlap, Genes)

# --- 2️⃣ LINCS (cztery glikokortykoidy) ---
df_lincs <- gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs %>%
  filter(Adjusted.P.value < 0.05) %>%
  filter(Term %in% c(
    "Methylprednisolone Up",
    "Dexamethasone Up",
    "Budesonide Up",
    "Hydrocortisone Up"
  )) %>%
  mutate(
    category = Term,
    source   = "LINCS_L1000",
    log10FDR = -log10(Adjusted.P.value)
  ) %>%
  select(source, category, log10FDR, Overlap, Genes)

# --- 3️⃣ CellMarker (komórki mózgu) ---
df_cell <- gene_summary_list$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>%
  filter(grepl("Astrocyte|Neuron|Microglial Cell|Oligodendrocyte", Term, ignore.case = TRUE)) %>%
  filter(grepl("Brain", Term, ignore.case = TRUE)) %>%
  mutate(
    category = case_when(
      str_detect(Term, regex("Astrocyte", ignore_case = TRUE)) ~ "Astrocyte",
      str_detect(Term, regex("Neuron", ignore_case = TRUE)) ~ "Neuron",
      str_detect(Term, regex("Microglial Cell", ignore_case = TRUE)) ~ "Microglia",
      str_detect(Term, regex("Oligodendrocyte Precursor Cell", ignore_case = TRUE)) ~ "OPC",
      str_detect(Term, regex("Oligodendrocyte", ignore_case = TRUE)) ~ "Oligodendrocyte",
      TRUE ~ Term
    ),
    source   = "CellMarker_2024",
    log10FDR = -log10(Adjusted.P.value)
  ) %>%
  select(source, category, log10FDR, Overlap, Genes)

# --- 4️⃣ Połączenie danych ---
df_plot_raw <- bind_rows(df_chea, df_lincs, df_cell)

# --- 5️⃣ Podsumowanie: liczba unikalnych genów i terminów na kategorię ---
df_plot <- df_plot_raw %>%
  group_by(source, category) %>%
  summarise(
    log10FDR = max(log10FDR, na.rm = TRUE),     # najsilniejszy sygnał
    n_terms = n(),
    n_genes = length(unique(unlist(str_split(Genes, ";"))))
  ) %>%
  ungroup()

# --- 6️⃣ Kolejność kategorii ---
cat_levels <- c(
  "NR3C1", "NR3C2",
  "Methylprednisolone Up", "Dexamethasone Up", "Budesonide Up", "Hydrocortisone Up",
  "Astrocyte", "Neuron", "Microglia", "Oligodendrocyte", "OPC"
)
df_plot$category <- factor(df_plot$category, levels = cat_levels)

# --- 7️⃣ Paleta kolorów ---
palette_sources <- c(
  "ChEA_2022" = "#4b0b0b",       # ciemny burgund
  "LINCS_L1000" = "#9c5757",     # ceglasty
  "CellMarker_2024" = "#f0c1c1"  # jasny różowo-kremowy
)

# --- 8️⃣ Wykres ---
p_all <- ggplot(df_plot, aes(x = category, y = log10FDR, fill = source)) +
  geom_col(width = 0.6, position = position_dodge(width = 0.7)) +  # bez czarnych kresek
  geom_text(
    aes(label = paste0(n_genes, " genes\n(", n_terms, " terms)")),
    vjust = -0.5, size = 3.5,
    position = position_dodge(width = 0.7)
  ) +
  scale_fill_manual(values = palette_sources) +
  labs(
    x = NULL,
    y = expression(-log[10]("FDR")),
    title = "Integrated enrichment summary (ChEA + LINCS + CellMarker)",
    subtitle = paste0(
      length(unique(unlist(str_split(df_plot_raw$Genes, ";")))),
      " unique genes across all datasets"
    ),
    fill = "Database"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(face = "bold", angle = 15, hjust = 1),
    legend.position = "top"
  ) +
  coord_cartesian(clip = "off", ylim = c(0, 60))

# --- 9️⃣ Wyświetlenie ---
p_all
