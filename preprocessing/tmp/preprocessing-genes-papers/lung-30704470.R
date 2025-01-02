# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")


sheets1 <-  excel_sheets("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM3_ESM.xlsx")
data1 <- lapply(sheets1, function(sheet) read_excel("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM3_ESM.xlsx", sheet = sheet))


data1[[1]] %>%
  select(-c(`Entrez ID`, peak))  %>% 
  rename_with(~ gsub(" ", "_", .)) %>%
  pivot_longer(
    cols = -Gene,
    names_to = c("Time", "Metric"),
    names_pattern = "(\\d+_h)_(.*)",
    values_to = "Value"
  ) %>%
  arrange(Gene, Time, Metric) %>% 
  spread(key = Metric, value = Value) %>%  
  set_colnames(c("gene_name", "time", "fold_change", "pvalue")) %>% 
  mutate(time = str_replace(time, "_", "")) %>% 
  as.data.frame() %>% 
  mutate(treatment = "budesonide",
         dose = "300nM",
         cell = "BEAS-2B",
         log2ratio = log2(fold_change),
         regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  select(c(gene_name, log2ratio, regulation, cell, treatment, dose, time, pvalue)) -> df1
  
sheets2 <-  excel_sheets("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM4_ESM.xlsx")
data2 <- lapply(sheets2, function(sheet) read_excel("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM4_ESM.xlsx", sheet = sheet))

data2[[1]] %>% 
  select(-`Entrez ID`) %>% 
  rename_with(~ gsub(" ", "_", .)) %>%  # Replace spaces with underscores in column names
  pivot_longer(
    cols = -Gene,
    names_to = c("Hormone", "Metric"),
    names_pattern = "([A-Za-z]+)_(.*)",  # Adjust the regex pattern to capture hormone names
    values_to = "Value"
  ) %>%
  arrange(Gene, Hormone, Metric) %>% 
  spread(key = Metric, value = Value) %>%  
  setNames(c("gene_name", "treatment", "fold_change", "pvalue")) %>%  # Rename columns
  as.data.frame() %>% 
  mutate(time = "6h",
         treatment = ifelse(treatment == "Bud", "budesonide", "dexamethasone"),
         dose = ifelse(treatment == "budesonide", "300nM", "1μM"),
         cell = "A549", 
         log2ratio = log2(fold_change),
         regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  select(c(gene_name, log2ratio, regulation, cell, treatment, dose, time, pvalue)) -> df2

sheets3 <-  excel_sheets("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM5_ESM.xlsx")
data3 <- lapply(sheets3, function(sheet) read_excel("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM5_ESM.xlsx", sheet = sheet))

data3[[1]] %>% 
  select(-c(`Entrez ID`, `Group (≥2 fold, P ≤ 0.05)`, `Group (≥1.25 fold)`)) %>% 
  rename_with(~ gsub(" ", "_", .)) %>%  # Replace spaces with underscores in column names
  pivot_longer(
    cols = -Gene,
    names_to = c("Cell_Type", "Metric"),
    names_pattern = "([A-Za-z0-9]+)_(.*)",  # Adjust the regex pattern to capture cell type names
    values_to = "Value"
  ) %>%
  arrange(Gene, Cell_Type, Metric) %>% 
  spread(key = Metric, value = Value) %>%  
  setNames(c("gene_name", "cell", "fold_change", "pvalue")) %>%  # Rename columns
  as.data.frame() %>% 
  mutate(time = "6h",
         treatment = "budesonide",
         dose = "300nM",
         log2ratio = log2(fold_change),
         regulation = ifelse(log2ratio > 0, "up", "down")) %>%
  select(c(gene_name, log2ratio, regulation, cell, treatment, dose, time, pvalue)) -> df3


sheets4 <-  excel_sheets("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM9_ESM.xlsx")
data4 <- lapply(sheets4, function(sheet) read_excel("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM9_ESM.xlsx", sheet = sheet))

data4[[1]] %>% 
  select(-c(`Entrez ID`, `Group (≤0.5 fold, P ≤ 0.05)`, `Group (≤0.8 fold)`)) %>% 
  rename_with(~ gsub(" ", "_", .)) %>%  # Replace spaces with underscores in column names
  pivot_longer(
    cols = -Gene,
    names_to = c("Cell_Type", "Metric"),
    names_pattern = "([A-Za-z0-9]+)_(.*)",  # Adjust the regex pattern to capture cell type names
    values_to = "Value"
  ) %>%
  arrange(Gene, Cell_Type, Metric) %>% 
  spread(key = Metric, value = Value) %>%  
  setNames(c("gene_name", "cell", "fold_change", "pvalue")) %>%  # Rename columns
  as.data.frame() %>% 
  mutate(time = "6h",
         treatment = "budesonide",
         dose = "300nM",
         log2ratio = log2(fold_change),
         regulation = ifelse(log2ratio > 0, "up", "down")) %>%
  select(c(gene_name, log2ratio, regulation, cell, treatment, dose, time, pvalue)) -> df4


sheets5 <-  excel_sheets("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM13_ESM.xlsx")
data5 <- lapply(sheets5, function(sheet) read_excel("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM13_ESM.xlsx", sheet = sheet))

data5[[1]] %>% 
  select(-c(`Entrez ID`, `Group (≥2 fold, P ≤ 0.05)`, `Group (≥1.25 fold)`)) %>% 
  rename_with(~ gsub(" ", "_", .)) %>%  # Replace spaces with underscores in column names
  pivot_longer(
    cols = -Gene,
    names_to = c("Cell_Type", "Metric"),
    names_pattern = "([A-Za-z0-9]+)_(.*)",  # Adjust the regex pattern to capture cell type names
    values_to = "Value"
  ) %>%
  arrange(Gene, Cell_Type, Metric) %>% 
  spread(key = Metric, value = Value) %>%  
  setNames(c("gene_name", "cell", "fold_change", "pvalue")) %>%  # Rename columns
  as.data.frame() %>% 
  mutate(time = "6h",
         treatment = "budesonide",
         dose = "300nM",
         log2ratio = log2(fold_change),
         regulation = ifelse(log2ratio > 0, "up", "down")) %>%
  select(c(gene_name, log2ratio, regulation, cell, treatment, dose, time, pvalue)) -> df5


sheets6 <-  excel_sheets("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM14_ESM.xlsx")
data6 <- lapply(sheets6, function(sheet) read_excel("data/supplement-genes-papers/lung-30704470/12920_2018_467_MOESM14_ESM.xlsx", sheet = sheet))

data6[[1]] %>% 
  select(-c(`Entrez ID`, `Group (≤0.5 fold, P ≤ 0.05)`, `Group (≤0.8 fold)`)) %>% 
  rename_with(~ gsub(" ", "_", .)) %>%  # Replace spaces with underscores in column names
  pivot_longer(
    cols = -Gene,
    names_to = c("Cell_Type", "Metric"),
    names_pattern = "([A-Za-z0-9]+)_(.*)",  # Adjust the regex pattern to capture cell type names
    values_to = "Value"
  ) %>%
  arrange(Gene, Cell_Type, Metric) %>% 
  spread(key = Metric, value = Value) %>%  
  setNames(c("gene_name", "cell", "fold_change", "pvalue")) %>%  # Rename columns
  as.data.frame() %>% 
  mutate(time = "6h",
         treatment = "budesonide",
         dose = "300nM",
         log2ratio = log2(fold_change),
         regulation = ifelse(log2ratio > 0, "up", "down")) %>%
  select(c(gene_name, log2ratio, regulation, cell, treatment, dose, time, pvalue)) ->df6


rbind(df1, df2, df3, df4, df5, df6) %>% dim

rbind(df1, df2, df3, df4, df5, df6) %>% unique() %>% dim

