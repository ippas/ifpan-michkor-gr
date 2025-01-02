source("preprocessing/functions/R/install-load-packages.R")

data_path <- "data/supplement-genes-papers/lung-33376135/146R1_Supplemental_Table_1.xlsx"

# Read sheet names and load each sheet into a list of data frames
sheets <- excel_sheets(data_path)
data <- lapply(sheets, function(sheet) read_excel(data_path, sheet = sheet, skip = 2))

process_data <- function(data_item, treatment_filter, peak_type, log2ratio_threshold = NULL) {
  processed <- data_item %>%
    bind_rows() %>%
    rename_with(~ str_replace_all(., " ", "_")) %>%
    unique() %>%
    pivot_longer(
      cols = c(ends_with("ANOVA"), ends_with("fold")),
      names_to = c("treatment", "time", "Type"),
      names_pattern = "(.*?)(\\d+h)?_(ANOVA|fold)"
    ) %>%
    pivot_wider(
      names_from = Type,
      values_from = value
    ) %>%
    mutate(treatment = str_replace_all(treatment, "_", "")) %>%
    mutate(fold = exp(fold)) %>%
    mutate(log2ratio = log2(fold)) %>%
    filter(ANOVA < 0.05)
  
  # Apply conditional filtering based on peak_type
  if (peak_type == "induction") {
    processed <- processed %>% filter(fold > 2)
  } else if (peak_type == "repression") {
    processed <- processed %>% filter(fold < 0.5)
  }
  
  # Additional log2ratio threshold filtering if applicable
  if (!is.null(log2ratio_threshold)) {
    processed <- processed %>% filter(log2ratio < log2ratio_threshold)
  }
  
  processed %>%
    select(Gene, treatment, time, ANOVA, log2ratio, starts_with("Peak")) %>%
    mutate(treatment = case_when(
      treatment == "Form" ~ "formoterol",
      treatment == "Bud" ~ "budesonide",
      TRUE ~ "budesonide_and_formoterol"
    )) %>%
    filter(treatment == treatment_filter) %>%
    mutate(peak_event = ifelse(peak_type == "induction", Peak_induction, Peak_repression)) %>%
    select(-starts_with("Peak")) %>%
    set_colnames(c("gene_name", "treatment", "time", "fdr", "log2ratio", "peak_event"))
}

# Define the treatments and peak types
treatments <- c("formoterol", "budesonide", "budesonide_and_formoterol", "formoterol", "budesonide", "budesonide_and_formoterol")
peak_types <- c("induction", "induction", "induction", "repression", "repression", "repression")

# Process each data frame and store in a list
processed_data <- lapply(seq_along(data), function(i) {
  process_data(data[[i]], treatments[i], peak_types[i])
})

processed_data 
# Combine all processed data frames into one
df_beas2b <- do.call(rbind, processed_data)  %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  mutate(cells = "BEAS-2B")
  


data_path <- "data/supplement-genes-papers/lung-33376135/146R1_Supplemental_Table_5.xlsx"

# Read sheet names and load each sheet into a list of data frames
sheets <- excel_sheets(data_path)
data <- lapply(sheets, function(sheet) read_excel(data_path, sheet = sheet, skip = 2))  
  
peak_types <- c("induction", "repression")

processed_data <- lapply(seq_along(data), function(i) {
  process_data(data[[i]], treatments[i], peak_types[i])
})


data[[1]] %>% 
  bind_rows() %>%
  rename_with(~ str_replace_all(., " ", "_")) %>%
  unique() %>%
  pivot_longer(
    cols = c(ends_with("pval"), ends_with("fold")),
    names_to = c("treatment", "time", "Type"),
    names_pattern = "(.*?)(\\d+h)?_(pval|fold)"
  ) %>% 
  pivot_wider(
    names_from = Type,
    values_from = value
  ) %>% 
  mutate(time = "6h") %>% 
  mutate(treatment = str_replace_all(treatment, "_", "")) %>%
  # mutate(fold = exp(fold)) %>%
  mutate(log2ratio = log2(fold)) %>%
  filter(pval < 0.05) %>%
  filter(fold > 2) %>%
  filter(treatment %in% c("Bud", "From", "BF")) %>% 
  select(-c(Additivity, `B+F-2`, `BF-1`, `Induced_by`, Additivity_diff, fold)) %>% 
  mutate(regulation = "up") %>% 
  set_colnames(c("gene_name", "treatment", "time", "fdr", "log2ratio", "regulation")) -> df_phbec_up
  
data[[2]] %>% 
  bind_rows() %>%
  rename_with(~ str_replace_all(., " ", "_")) %>%
  unique() %>%
  pivot_longer(
    cols = c(ends_with("pval"), ends_with("fold")),
    names_to = c("treatment", "time", "Type"),
    names_pattern = "(.*?)(\\d+h)?_(pval|fold)"
  ) %>% 
  pivot_wider(
    names_from = Type,
    values_from = value
  ) %>% 
  mutate(time = "6h") %>% 
  mutate(treatment = str_replace_all(treatment, "_", "")) %>%
  mutate(log2ratio = log2(fold)) %>%
  filter(pval < 0.05) %>% 
  filter(fold < 0.5) %>% 
  filter(treatment %in% c("Bud", "From", "BF")) %>% 
  select(-c(Repressed_by, Enhanced_repression, `BF-Bud_diff`, `BF-Form_diff`, fold)) %>% 
  mutate(regulation = "down") %>% 
  set_colnames(c("gene_name", "treatment", "time", "fdr", "log2ratio", "regulation")) -> df_phbec_down

rbind(df_phbec_up, df_phbec_down) %>% 
  mutate(treatment = case_when(
    treatment == "Form" ~ "formoterol",
    treatment == "Bud" ~ "budesonide",
    TRUE ~ "budesonide_and_formoterol"
  )) %>%
  mutate(cells = "pHBECs") -> df_phbec

rbind(df_beas2b, df_phbec) %>% 
  mutate(dose = ifelse(treatment == "formoterol", "10nM", 
                            ifelse(treatment == "budesonide", "100nM", "100nM_and_10nM"))) -> df_preprocessing

data.frame(
  article_source = "https://pubmed.ncbi.nlm.nih.gov/33376135/",
  pmid = "33376135",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "lung",
  cell = df_preprocessing$cells,
  environment = "in-vitro",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = df_preprocessing$time,  # Replace NA with actual data if available
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = 0.05,
  abs_log2ratio_threshold = 1,
  method = "microarray",
  statistical_method = "eBayes_ANOVA",
  treatment_type = "acute"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-33376135.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )

rm(
  treatments,
  peak_types,
  processed_data,
  process_data,
  df_phbec_up,
  df_phbec_down,
  df_beas2b,
  df_phbec
)





