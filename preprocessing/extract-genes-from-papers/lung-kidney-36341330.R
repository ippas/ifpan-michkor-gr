source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

data_path <- "data/supplement-genes-papers/lung-kidney-36341330/DataSheet_1.xlsx"

# Load the data, skipping the first two rows and using the third row as headers
read_excel(data_path, skip = 2, sheet = "Lung")

# Load and process the second and third rows for column names
header_part1 <- read_excel(data_path, range = "B2:Q2", col_names = FALSE) %>%
  as.character() %>%
  purrr::keep(~ .x != "NA") %>%
  rep(each = 4) %>%
  str_replace_all(., " ", "_")

header_part2 <- read_excel(data_path, range = "B3:Q3", col_names = FALSE) %>%
  as.character() %>%
  str_replace_all(., " ", "_")

# Combine and set column names
data_colnames <- c("gene_name", paste(header_part1, header_part2, sep = "|"))

# Transform the data and rename columns
read_excel(data_path, skip = 2, sheet = "Lung") %>%
  setNames(data_colnames) %>%
  gather(key = "columns", value = "value", -gene_name) %>%
  separate(columns, into = c("comparison", "result"), sep = "\\|") %>%
  mutate(comparison = str_replace_all(comparison, "\\.", "")) %>%
  filter(result != "-log10_p_val") %>%
  spread(key = result, value = value) %>%
  setNames(c("gene_name", "comparison", "fdr", "log2ratio", "pvalue")) %>% 
  mutate(time = "14weeks") %>% 
  mutate(tissue = "lung-cSiO2") -> df_lung

read_excel(data_path, skip = 2, sheet = "Kidney") %>%
  setNames(data_colnames) %>%
  gather(key = "columns", value = "value", -gene_name) %>%
  separate(columns, into = c("comparison", "result"), sep = "\\|") %>%
  mutate(comparison = str_replace_all(comparison, "\\.", "")) %>%
  filter(result != "-log10_p_val") %>%
  spread(key = result, value = value) %>%
  setNames(c("gene_name", "comparison", "fdr", "log2ratio", "pvalue")) %>% 
  mutate(time = "14weeks") %>% 
  mutate(tissue = "kidney-cSiO2") -> df_kidney

# prepare data for blood
# Load and process the second and third rows for column names
header_part1 <- read_excel(data_path, range = "B2:BM2", col_names = FALSE, sheet = "Whole blood") %>%
  as.character() %>%
  purrr::keep(~ .x != "NA") %>%
  rep(each = 16) %>%
  str_replace_all(., " ", "_")

header_part2 <- read_excel(data_path, range = "B3:BM3", col_names = FALSE, sheet = "Whole blood") %>%
  as.character() %>%
  purrr::keep(~ .x != "NA") %>% 
  rep(each = 4) %>%
  str_replace_all(., " ", "_")

header_part3 <- read_excel(data_path, range = "B4:BM4", col_names = FALSE, sheet = "Whole blood") %>%
  as.character() %>%
  str_replace_all(., " ", "_")


# Combine and set column names
data_colnames <- c("gene_name", paste(header_part1, header_part2, header_part3, sep = "|"))


read_excel(data_path, skip = 4, sheet = "Whole blood", col_names = FALSE) %>%
  setNames(data_colnames) %>% 
  gather(key = "columns", value = "value", -gene_name) %>% 
  separate(columns, into = c("time", "comparison", "result"), sep = "\\|") %>%
  mutate(comparison = str_replace_all(comparison, "\\.", "")) %>% 
  filter(result != "-log10_p-val") %>%
  spread(key = result, value = value) %>% 
  setNames(c("gene_name", "time", "comparison", "fdr", "log2ratio", "pvalue")) %>% 
  mutate(time = ifelse(time == "Week_11", "11weeks", 
                       ifelse(time == "Week_7", "7weeks",
                              ifelse(time == "Week_9", "9weeks", time)))) %>% 
  select(c(gene_name, comparison, fdr, log2ratio, pvalue, time)) %>% 
  mutate(tissue = "blood-cSiO2") -> df_blood



rbind(df_kidney, df_lung, df_blood) %>% 
  filter(fdr < 0.05) %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  filter(abs_log2ratio > log2(1.5)) %>% 
  mutate(regulation = ifelse(log2ratio  > 0, "up", "down")) %>% 
  mutate(time = tolower(time)) %>% 
  mutate(treatment = ifelse(comparison == "cSiO2_P0_vs_VEH_P0", NA, "prednisone")) %>% 
  mutate(dose = ifelse(comparison == "cSiO2_P0_vs_VEH_P0", "0mg/kg", 
                       ifelse(comparison == "cSiO2_PM_vs_VEH_P0", "15mg/kg",
                              ifelse(comparison == "cSiO2_PL_vs_VEH_P0", "5mg/kg", "5mg/kg_vs_15mg/kg")))) -> df_preprocessing


data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9627297/",
  pmid = "36341330",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "mouse",
  tissue = df_preprocessing$tissue,
  cell = NA,
  environment = "in-vivo",
  treatment = df_preprocessing$treatment,
  dose = df_preprocessing$dose,  # Replace NA with actual data if available
  time = df_preprocessing$time,  # Replace NA with actual data if available
  pvalue = df_preprocessing$pvalue,  # Replace with actual column name
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = "0.05",
  abs_log2ratio_threshold = 0.58,
  method = "NanoString",
  strain = "NZBWF1",  # Replace with actual column name
  comparison = df_preprocessing$comparison,
  statistical_method = "nSolver_Advanced_Analysis_Module,",
  treatment_type = "chronic"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/lung-kidney-36341330.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
