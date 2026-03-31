library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)

# --- ChEA_2022 ---
df_chee <- gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$ChEA_2022 %>%
  filter(Adjusted.P.value < 0.05, n_genes >= 3) %>%
  filter(grepl("NR3C1|NR3C2", Term, ignore.case = TRUE)) %>%
  mutate(
    aggregate_term = case_when(
      grepl("NR3C1", Term, TRUE) ~ "NR3C1",
      grepl("NR3C2", Term, TRUE) ~ "NR3C2",
      TRUE ~ "other"
    )
  ) %>%
  group_by(aggregate_term) %>%
  nest() %>%
  mutate(
    aggregate_genes    = map(data, ~ .x$Genes |> unlist() |> strsplit(";") |> unlist() |> unique()),
    n_aggreate_genes   = map_int(aggregate_genes, length),
    n_traits           = map_int(data, nrow),
    source             = "ChEA_2022"
  )

# --- LINCS_L1000 ---
df_lincs <- gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$LINCS_L1000_Chem_Pert_Consensus_Sigs %>%
  filter(Adjusted.P.value < 0.05, n_genes >= 3) %>%
  filter(grepl("Betamethasone|Dexamethasone|Fluocortolone|Methylprednisolone|Paramethasone|Prednisolone|Prednisone|Triamcinolone|Hydrocortisone|Cortisone|Prednylidene|Rimexolone|Deflazacort|Cloprednol|Meprednisone|Cortivazol|Vamorolone",
               Term, TRUE)) %>%
  mutate(
    aggregate_term = case_when(
      grepl("up",   Term, TRUE) ~ "GC ↑",
      grepl("down", Term, TRUE) ~ "GC ↓",
      TRUE ~ "other"
    )
  ) %>%
  group_by(aggregate_term) %>%
  nest() %>%
  mutate(
    aggregate_genes    = map(data, ~ .x$Genes |> unlist() |> strsplit(";") |> unlist() |> unique()),
    n_aggreate_genes   = map_int(aggregate_genes, length),
    n_traits           = map_int(data, nrow),
    source             = "LINCS_L1000"
  )

# --- CellMarker_2024 ---
df_cell <- gene_summaryList_blood$n10_short_time$summary_down$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>%
  filter(Adjusted.P.value < 0.05, n_genes >= 3) %>%
  filter(grepl("lung|blood|brain", Term, TRUE)) %>%
  mutate(
    aggregate_term = case_when(
      grepl("brain", Term, TRUE) ~ "Brain cells",
      grepl("blood", Term, TRUE) ~ "Blood cells",
      grepl("lung",  Term, TRUE) ~ "Lung cells",
      TRUE ~ "other"
    )
  ) %>%
  group_by(aggregate_term) %>%
  nest() %>%
  mutate(
    aggregate_genes    = map(data, ~ .x$Genes |> unlist() |> strsplit(";") |> unlist() |> unique()),
    n_aggreate_genes   = map_int(aggregate_genes, length),
    n_traits           = map_int(data, nrow),
    source             = "CellMarker_2024"
  )

# --- łączymy
combined_df <- bind_rows(df_chee, df_lincs, df_cell)

# --- stała KOLEJNOŚĆ na osi X
term_levels <- c("NR3C1","NR3C2","GC ↑","GC ↓","Brain cells","Blood cells","Lung cells")
combined_df <- combined_df %>%
  mutate(aggregate_term = factor(aggregate_term, levels = term_levels))

# --- etykiety z (n = ...) bez nazw baz
label_map <- combined_df %>%
  select(aggregate_term, n_traits) %>%
  distinct() %>%
  mutate(label = paste0(as.character(aggregate_term), " (n = ", n_traits, ")")) %>%
  { setNames(.$label, .$aggregate_term) }

# --- kolory źródeł
palette_sources <- c(
  "ChEA_2022"     = "#4b0b0b",
  "LINCS_L1000"   = "#9c5757",
  "CellMarker_2024" = "#f0c1c1"
)

# --- wykres
ggplot(combined_df, aes(x = aggregate_term, y = n_aggreate_genes, fill = source)) +
  geom_bar(stat = "identity", width = 0.6, color = "black") +
  scale_x_discrete(labels = label_map) +
  scale_fill_manual(values = palette_sources) +
  theme_classic(base_size = 16) +
  labs(
    title = "Σ GR-related aggregated enrichments (Blood, n10_short_time)",
    x = "Aggregated term (n = number of enrichments)",
    y = expression(Sigma~"unique genes per aggregate")
  ) +
  theme(
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 25, hjust = 1)
  )
