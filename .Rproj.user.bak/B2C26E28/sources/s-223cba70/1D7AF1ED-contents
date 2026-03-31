chea_df_raw %>%
  mutate(
    TF = stringr::str_split(Term, " ", simplify = TRUE)[, 1],
    logP = -log10(Adjusted.P.value)
  ) %>%
  group_by(cluster) %>%
  arrange(Adjusted.P.value) %>%
  # slice_head(n = 10) %>%
  mutate(TF_rank = paste0("top", row_number())) %>%
  ungroup() %>% 
  mutate(
    nuclear_receptors = TF %in% nuclear_receptors$hgnc_symbol,
    minuslog10P = ifelse(logP > 10, 10, logP)
  ) %>% 
  filter(P.value < 0.05) %>% 
  # filter(Adjusted.P.value < 0.2) %>% 
  select(-c(Old.P.value, Old.Adjusted.P.value, logP)) -> data_supplement_preprocessing



# 
mapper_cluster2tissue <- read.delim("data/review-dextis/mapper_cluster2tissue-08.04.2025.tsv", 
           sep = "\t")

mapper_TFchea2mouse <- read.delim("data/review-dextis/all-signif-TF-mapping-chea.tsv", 
           sep = "\t")

all_signif_TF_chea_mean_sal <- read.delim("data/review-dextis/all-signif-TF-chea-mean-sal.tsv", 
           sep = "\t")

# read data
gtex_tissue_gene_expression <- read.delim(
  "data/databases/GTEx/gene-expression-selected-tissues-threshold1.tsv",
  sep = "\t",
  header = TRUE
)

gtex_tissue_gene_expression %>% 
  select(c(label, gene_name, mean_expression)) %>% 
  mutate(tissue = case_when(
    label == "adrenal_gland" ~ "ADR",
    label == "adipose_visceral_omentum" ~ "FAT",
    label == "liver" ~ "LIV",
    label == "muscle_skeletal" ~ "MUS",
    label == "lung" ~ "LUN",
    label == "kidney_cortex" ~ "KID",  # jeśli masz też nerkę, dodaj właściwą etykietę
    label == "pituitary" ~ "PIT",
    label == "spleen" ~ "SPL",
    label == "brain_hypothalamus" ~ "HTH",
    TRUE ~ NA_character_
  )) %>% 
  # head %>% 
  mutate(mean_cpm = mean_expression*10) %>% 
  select(-c(mean_expression, label)) -> gtex_tissue_gene_expression
# filter(mean_cpm > 50) 

gtex_tissue_gene_expression %>% 
  filter(gene_name %in% mapper_TFchea2mouse$mapped_name) -> gtex_tissue_expression_preprocessing



all_signif_TF_chea_mean_sal %>% 
  group_by(chea_name_TF, tissue) %>% 
  nest() %>% 
  mutate(
    data = data %>% map(~ 
                          .x %>% 
                          slice_max(mean, n = 1, with_ties = FALSE)
    )
  ) %>% 
  unnest(data) %>% 
  select(-gene_name) %>% 
  as.data.frame() %>%
  # filter(mean > 8) %>% 
  rename(mean = "mean_log2_microarrays") -> all_signif_TF_chea_mean_sal_preprocessing

all_signif_TF_chea_mean_sal_preprocessing %>% head
data_supplement_preprocessing %>% 
  left_join(., mapper_cluster2tissue, by = "cluster") %>% 
  left_join(., mapper_TFchea2mouse, by = c("TF" = "original_name")) %>% 
  left_join(., gtex_tissue_expression_preprocessing, by = c("mapped_name" = "gene_name", "tissue" = "tissue")) %>% 
  unique() %>% 
  left_join(., all_signif_TF_chea_mean_sal_preprocessing, by = c("tissue" = "tissue", "TF" = "chea_name_TF")) %>% 
  rename(mean_cpm = "mean_cpm_gtex") %>% 
  mutate(tissue_expression_gtex = ifelse(mean_cpm_gtex <= 50 | is.na(mean_cpm_gtex), F, T )) %>% 
  mutate(tissue_expression_microarrays = ifelse(mean_log2_microarrays <= 8 | is.na(mean_log2_microarrays), F, T)) -> data_supplement_preprocessing 

readr::write_tsv(
  data_supplement_preprocessing,
  file = "data/review-dextis/supplementary_table_fig3_review.tsv"
)


data_supplement_preprocessing %>% 
  filter(P.value < 0.05) %>% 
  filter(cluster == "cluster_C") %>% 
  as.data.frame()
  

# Oczekiwana kolejność klastrów
cluster_order <- c(paste0("cluster_", LETTERS[1:16]), "DOWN", "UP")

# Filtrowanie danych
filtered_df <- data_supplement_preprocessing %>%
  filter(P.value < 0.05)


filtered_df %>% colnames

filtered_df %>%
  # select(-c(tissue_expression_gtex, tissue_expression_microarrays)) %>%
  mutate(
    tissue_cpm_pair = paste0(tissue, ":", mean_cpm_gtex),
    tissue_microarray_pair = paste0(tissue, ":", mean_log2_microarrays)
  ) %>%
  group_by(TF) %>% colnames

filtered_df %>% 
  as.data.frame() %>%
  mutate(
    dplyr::across(
      c(mean_cpm_gtex, mean_log2_microarrays, Odds.Ratio, Combined.Score),
      ~ round(., 3)
    )
    # dplyr::across(
    #   c(P.value, Adjusted.P.value),
    #   ~ formatC(., format = "e", digits = 2)
    # )
  ) %>% 
  # select(-c(tissue_expression_gtex, tissue_expression_microarrays)) %>%
  mutate(
    tissue_cpm_pair = paste0(tissue, ":", mean_cpm_gtex),
    tissue_microarray_pair = paste0(tissue, ":", mean_log2_microarrays)
  ) %>%
  group_by(TF) %>% 
  nest() %>% 
  # nest(data = c(tissue, mapped_name, tissue_cpm_pair, tissue_microarray_pair)) %>%
  mutate(
    tissue_cpm_gtex = purrr::map_chr(data, ~ paste(unique(.x$tissue_cpm_pair), collapse = "|")),
    tissue_log2_microarray = purrr::map_chr(data, ~ paste(unique(.x$tissue_microarray_pair), collapse = "|")),
    pattern_tissue_colocalization_gtex = purrr::map_lgl(data, ~ sum(.x$tissue_expression_gtex, na.rm = TRUE) > 0),
    pattern_tissue_colocalization_microarray = purrr::map_lgl(data, ~ sum(.x$tissue_expression_microarrays, na.rm = TRUE) > 0)
  ) %>% 
  unnest(data) %>%
  ungroup %>% 
  select(-c(tissue_expression_gtex, tissue_expression_microarrays)) %>%
  select(!c(tissue, mean_cpm_gtex, mean_log2_microarrays, tissue_cpm_pair, tissue_microarray_pair)) %>% unique -> filtered_df

filtered_df %>% 
  select(-minuslog10P) %>% 
  mutate(P.value = as.numeric(P.value)) %>% 
  mutate(Adjusted.P.value = as.numeric(Adjusted.P.value)) -> filtered_df

# Tworzenie listy – tylko dostępne klastry
cluster_split <- split(filtered_df, filtered_df$cluster)

cluster_split$cluster_A

# Dodanie brakujących klastrów jako pustych data.frame
cluster_list <- lapply(cluster_order, function(clust_name) {
  if (clust_name %in% names(cluster_split)) {
    cluster_split[[clust_name]]
  } else {
    data.frame()
  }
})

cluster_list %>% str

# Nadanie nazw elementom listy
names(cluster_list) <- cluster_order

# Usunięcie przedrostka "cluster_"
names(cluster_list) <- sub("^cluster_", "", names(cluster_list))

library(openxlsx)

# Utworzenie nowego workbooka
wb <- createWorkbook()

# Dodanie każdego data.frame jako osobnego arkusza
for (i in seq_along(cluster_list)) {
  sheet_name <- names(cluster_list)[i]
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet = sheet_name, cluster_list[[i]])
}


saveWorkbook(wb, file = "data/review-dextis/supplementary_table3-new.xlsx", overwrite = TRUE)

# Styl naukowy
sci_style <- createStyle(numFmt = "0.00E+00")

# Iteracja po wszystkich arkuszach
for (sheet_name in names(wb)) {
  # Przeczytaj dane z arkusza
  df <- readWorkbook(wb, sheet = sheet_name)
  
  # Kolumny do sformatowania
  target_cols <- c("P.value", "Adjusted.P.value")
  
  for (col_name in target_cols) {
    col_index <- which(names(df) == col_name)
    
    if (length(col_index) == 1) {
      addStyle(wb, sheet = sheet_name, style = sci_style,
               rows = 2:(nrow(df) + 1), cols = col_index, gridExpand = TRUE)
    }
  }
}

# Zapisz workbook
saveWorkbook(wb, file = "data/review-dextis/supplementary_table3-new.xlsx", overwrite = TRUE)


cluster_list$cluster_A
  
data_supplement_preprocessing %>%
  filter(cluster == "cluster_A") %>% 
  filter(TF == "STAT3")
  filter(P.value < 0.05) %>% dim
  dim
  filter(TF_rank %in% c("top1", "top2", "top3", "top4", "top5", "top6", "top7", "top8", "top9", "top10")) %>% 
  .$P.value %>% max

data_supplement_preprocessing %>% 
  filter(cluster %in% c("cluster_H", "cluster_I", "cluster_L")) %>% 
  filter(TF == "LXR") %>% 
  filter(mean_log2_microarrays > 8) %>% .$tissue %>% unique
  filter(mean_cpm_gtex > 50)

  
  
# check 0.01
filtered_df %>% 
  filter(P.value < 0.01)
