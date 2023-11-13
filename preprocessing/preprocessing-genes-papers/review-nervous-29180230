# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# read data
summary_papers <- read.table("data/supplement-genes-papers/review-nervous-29180230//summary-table-papers-Arkusz1.tsv", header = T, sep = "\t") %>% 
  mutate(cleaned_year = gsub("[^0-9]", "", year)) %>% 
  mutate(label = paste(references, cleaned_year, sep = "_")) %>% 
  mutate(label = ifelse(references == "Morsink", paste0(label, "_", tissue_cells), label)) %>% 
  select(-c(species, drug, dose, cleaned_year, year, references))

read_excel("data/supplement-genes-papers/review-nervous-29180230//1-s2.0-S0278584617304633-mmc1.xlsx", sheet = 1) %>% 
  as.data.frame() %>% 
  select(-matches("^\\.\\.\\.[0-9]+$")) %>% 
  rename_all(~gsub("[^a-zA-Z0-9]", "_", .)) %>% 
  rename_all(~gsub("_+", "_", .)) %>% 
  rename_all(tolower) %>% 
  mutate(part_of_brain = gsub(" ", "-", part_of_brain)) %>% 
  mutate(subregion_cell_type = gsub(" ", "-", subregion_cell_type)) %>% 
  mutate(subregion_cell_type = tolower(subregion_cell_type)) %>% 
  mutate(first_author = ifelse(first_author == "sato", "Sato", first_author)) %>% 
  mutate(first_author = ifelse(first_author == "SLEZAK", "Slezak", first_author)) %>% 
  mutate(label = paste(first_author, year, sep = "_")) %>% 
  mutate(label = ifelse(first_author == "Morsink", paste0(label, "_", part_of_brain), label))  -> review_nervous



inner_join(summary_papers, review_nervous, by = "label") %>%
  rename(article_source = url) %>%
  select(-c(tissue_cells, label, first_author, year)) %>%
  rename(tissue = part_of_brain) %>%
  mutate(dose = sub(" ", "", dose)) %>%
  rename(time = lattency_of_tissue_collection) %>%
  mutate(time = sub(" ", "", time)) %>%
  rename(regulation = response_to_gcs) %>%
  select(-gene_name) %>%
  rename(gene_name = gene_symbol) %>%
  mutate(ensembl_id = NA,
         log2ratio = NA) %>%
  rename(cell = subregion_cell_type) %>%
  select(
    c(
      article_source,
      pmid,
      gene_name,
      ensembl_id,
      log2ratio,
      regulation,
      species,
      tissue,
      cell,
      treatment,
      dose,
      time,
      fdr_threshold,
      comparision,
      input_for_standarisation,
      input_type,
      method
    )
  ) %>%
  mutate(treatment = tolower(treatment)) %>% 
  mutate(tissue = tolower(tissue)) %>% 
  mutate(species = tolower(species)) %>% 
  mutate(species = ifelse(species == "rats", "rat", species)) %>% 
  filter(pmid != 24147833) %>%   
  filter(pmid != 24777604) %>% 
  mutate(species = ifelse(species == "mice", "mouse", species)) -> nervous_cells_database


write.table(
  nervous_cells_database,
  file = "data/supplement-genes-papers/review-nervous-29180230/nervous-cells-review-29180230.tsv",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)


