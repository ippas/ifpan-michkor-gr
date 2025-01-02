# Load necessary packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Load the first three sheets from the specified Excel file into a list of data frames
sheets <- excel_sheets("data/supplement-genes-papers/NPSCs-KO-Cav1/MCB.01121-13_zmb999100488sd1.xlsx")[1:3]
data <- lapply(sheets, function(sheet) read_excel("data/supplement-genes-papers/NPSCs-KO-Cav1/MCB.01121-13_zmb999100488sd1.xlsx", sheet = sheet))

# Rename sheet names for clarity
sheets[1] <- gsub("Cav KO Dex vs Vehicle", "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", sheets[1])
sheets[2] <- gsub("WT C57 Dex vs Vehicle", "C57-Dex_vs_C57-Ethanol", sheets[2])
sheets[3] <- gsub("Cav KO vs WT C57 Vehicle", "Cav1-KO-Ethanol_vs_C57BL6-Dex", sheets[3])

# Assign cleaned sheet names to data list
names(data) <- sheets

# Process each sheet to extract relevant columns and perform various transformations
npscs_ko_cav1_df <- data %>%
  imap_dfr(function(x, sheet_name) {
    x[, c(3, 4, 7)] %>%
      setNames(c("gene_name", "pvalue", "fold_change")) %>%
      mutate(pvalue = as.numeric(gsub("<", "", pvalue))) %>%
      mutate(comparison = sheet_name)
  }) %>% 
  as.data.frame() %>%
  mutate(
    pvalue = as.numeric(pvalue),
    strain = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "C57BL/6", 
                    ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "C57BL/6-KO-Cav1", "C57BL/6|C57BL/6-KO-Cav1")),
    log2ratio = log2(fold_change),
    regulation = ifelse(log2ratio > 0, "up", "down"),
    cell = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "NPSCs", 
                  ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "NPSCs-KO-Cav1", "NPSCs|NPSCs-KO-Cav1")),
    treatment = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "dexamethasone", 
                       ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "dexamethasone", "vehicle-ethanol"))
  )

# Create the final database with all relevant metadata and results
npscs_ko_cav1_database <- data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4097667/",
  pmid = 24777604,
  gene_name = npscs_ko_cav1_df$gene_name,
  ensembl_id = NA,
  log2ratio = npscs_ko_cav1_df$log2ratio,
  regulation = npscs_ko_cav1_df$regulation,
  species = "mouse", 
  tissue = "embryos_hypothalamic-region",
  cell = npscs_ko_cav1_df$cell,
  environment = "in-vitro",
  treatment = npscs_ko_cav1_df$treatment,
  dose = "100nM",
  time = "4h",
  pvalue = npscs_ko_cav1_df$pvalue,
  fdr_threshold = 0.1,
  method = "microarray",
  comparison = npscs_ko_cav1_df$comparison,
  strain = npscs_ko_cav1_df$strain,
  statistical_method = "BRB-ArrayTools",
  treatment_type = "acute"
)


npscs_ko_cav1_database %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/NPSCs-KO-Cav1.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )
