# Load required packages and functions
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read the names of the Excel sheets from the given file
sheets <- excel_sheets("data/supplement-genes-papers/ASM/pone.0099625.s014.xlsx")

# Load all sheets from the Excel file into a list of data frames
data <- lapply(sheets, function(sheet) {
  read_excel("data/supplement-genes-papers/ASM/pone.0099625.s014.xlsx", sheet = sheet)
})

# # Perform initial data transformation on the loaded data
# # Here we are excluding some columns and converting Ln.Fold.Change to log2 scale
# data %>% 
#   as.data.frame() %>% 
#   select(-c(Dex.FPKM, Untreated.FPKM, Test.Statistic, Locus)) %>% 
#   arrange(Ln.Fold.Change.) %>% 
#   mutate(log2ratio = Ln.Fold.Change./log(2))

# I suppose that in the paper authors make mistake between log2 and ln
# Further process the data to include only needed columns
# Also, add a new 'regulation' column based on 'log2ratio'
data %>% 
  as.data.frame() %>% 
  select(-c(Dex.FPKM, Untreated.FPKM, Test.Statistic, Locus))  %>% 
  set_names(c("gene_name", "log2ratio", "pvalue", "fdr")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> asm_df

# Create the final data frame for the database
# This includes metadata like article source, PMIDs, and experimental conditions
data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4057123/",
  pmid = 24926665,
  gene_name = asm_df$gene_name,
  ensembl_id = NA,
  log2ratio = asm_df$log2ratio,
  regulation = asm_df$regulation,
  species = "human",
  tissue = "lung",
  cell = "ASM", 
  environment = "in-vitro",
  treatment = "dexamethasone",
  dose = "1μM",
  time = "18h",
  pvalue = asm_df$pvalue,
  fdr = asm_df$fdr,
  fdr_threshold = 0.05,
  method = "RNA-seq",
  treatment_type = "acute"
) -> asm_database  # Store the final data frame in 'asm_database'
