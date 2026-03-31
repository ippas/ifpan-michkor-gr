library(dplyr)
library(purrr)
library(tidyr)


papers_data_preprocessing %>% 
  filter(simple_tissue == "blood") -> AllBloodGeneDf

AllBloodGeneDf$time %>% unique

AllBloodGeneDf %>% filter(time == "main_effect_of_treatments")

# ##############################################################################
# ---- prepare datasets ----
# ##############################################################################

dataset_list_blood <- list(
  
  # 1️⃣ Pełny zbiór
  all = AllBloodGeneDf,
  
  # 2️⃣ Listy z co najmniej 10 genami
  n10 = AllBloodGeneDf %>%
    filter(n_genes >= 10),
  
  
  # 12️⃣ Listy z efektami ostrymi (bez długiego czasu)
  short_time = AllBloodGeneDf %>%
    filter(!(time %in% c("main_effect_of_treatments",
                         "9weeks", "7weeks", "168h"))),

  n10_short_time = AllBloodGeneDf %>%
    filter(n_genes >= 10) %>% 
    filter(!(time %in% c("main_effect_of_treatments",
                         "9weeks", "7weeks", "168h")))
)


gene_summaryList_blood <- purrr::imap(dataset_list_blood, ~ {
  message("🔍 Processing dataset: ", .y)
  summarize_gene_dataset_and_enrichr(
    df = .x,
    dataset_name = .y,
    gene_list_validation = marpiech_clusters,  # lub inny zbiór referencyjny
    enrichr_run = TRUE                         # uruchamia automatycznie sekcję Enrichr
  )
})


gene_summaryList_blood$n10_short_time$summary_up$freq_genes_per_paper

 
gene_summaryList_blood$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% 
  filter(n_genes >= 3) %>% 
  filter(Adjusted.P.value < 0.05) %>%
  # select(-Genes) %>%
  filter(grepl("Lung", Term, ignore.case = T))  %>% 
  .$Genes %>% 
  strsplit(";") %>% 
  unlist %>% unique
  filter(grepl("Blood|brain|macrophage", Term, ignore.case = TRUE)) 
  
  
gene_summaryList_lung$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% 
  filter(n_genes >= 3) %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(grepl("Brain", Term, ignore.case = T)) 
  .$Genes %>% 
  strsplit(";") %>% 
  unlist %>% unique



summ$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% 
  filter(n_genes >= 3) %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(grepl("Brain", Term, ignore.case = T)) %>% 
  .$Genes %>% 
  strsplit(";") %>% 
  unlist %>% unique


gene_summaryList_brain$n10_short_time$summary_up$enrichr_summary_source$enrichr_3pub$CellMarker_2024 %>% 
  filter(n_genes >= 3) %>% 
  filter(Adjusted.P.value < 0.05) %>% 
  filter(grepl("", Term, ignore.case = T)) 
  