library(dplyr)
library(purrr)
library(tidyr)

# ##############################################################################
# ---- prepare datasets ----
# ##############################################################################

dataset_list_lung <- list(
  
  # 1️⃣ Pełny zbiór
  all = AllLungGeneDf,
  
  # 2️⃣ Listy z co najmniej 10 genami
  n10 = AllLungGeneDf %>%
    filter(n_genes >= 10),
  

  # 12️⃣ Listy z efektami ostrymi (bez długiego czasu)
  short_time = AllLungGeneDf %>%
    filter(!(time %in% c("14weeks"))),
  
  n10_short_time = AllLungGeneDf %>%
  filter(n_genes >= 10) %>% 
    filter(!(time %in% c("14weeks")))
)


gene_summaryList_lung <- purrr::imap(dataset_list_lung, ~ {
  message("🔍 Processing dataset: ", .y)
  summarize_gene_dataset_and_enrichr(
    df = .x,
    dataset_name = .y,
    gene_list_validation = marpiech_clusters,  # lub inny zbiór referencyjny
    enrichr_run = TRUE                         # uruchamia automatycznie sekcję Enrichr
  )
})
