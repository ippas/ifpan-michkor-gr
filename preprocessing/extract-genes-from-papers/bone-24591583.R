source("preprocessing/functions/R/install-load-packages.R")


################################################################################
data_path <- "data/supplement-genes-papers/bone-24591583/1400522111_sd01.xlsx"

read_excel_sheets(file_path = data_path, skip_rows = 1) -> data

names(data) <- c("dexamethasone", "dexamethasone_Hic5-modulated", "dexamethasone_Hic5-independent", "dexamethasone_Hic5-blocked")

lapply(data, dim)

# Load the data, skipping the first two rows and using the third row as headers
read_excel(data_path, skip = 1, sheet = "Sheet1") -> data

data[[1]] %>%
  as.data.frame() %>% 
  .[, c(1,2,5,6)] %>% 
  set_colnames(c("probe_id", "gene_name", "log2ratio", "fdr")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
  filter(gene_name != "NA")


lapply(data, function(x){
  x %>% 
    as.data.frame() %>% 
    .[, c(1,2,5,6)] %>% 
    set_colnames(c("probe_id", "gene_name", "log2ratio", "fdr")) %>% 
    mutate(regulation = ifelse(log2ratio > 0, "up", "down")) %>% 
    filter(gene_name != "NA")
}) %>% bind_rows(.id = "treatment") -> df_preprocessing

data.frame(
  article_source = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3964041/",
  pmid = "24591583",
  gene_name = df_preprocessing$gene_name,  # Replace with actual column name
  ensembl_id = NA,  # Replace NA with actual data if available
  log2ratio = df_preprocessing$log2ratio,  # Replace with actual column name
  regulation = df_preprocessing$regulation,  # Replace with actual column name
  species = "human",
  tissue = "bone",
  cell = "U2OS",
  environment = "in-vitro",
  treatment = df_preprocessing$treatment,
  dose = "100nM",  # Replace NA with actual data if available
  time = "4h",  # Replace NA with actual data if available
  fdr = df_preprocessing$fdr,  # Replace with actual column name
  fdr_threshold = 0.02,
  abs_log2ratio_threshold = 1,
  method = "microarray",
  statistical_method = "limma",
  treatment_type = "acute",
  geo = "GSE46448"
) %>% 
  write.table(
    .,
    file = "data/supplement-genes-papers/preprocessing-genes-papers/bone-24591583.tsv",
    quote = FALSE,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE
  )


# # Install BiocManager if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# # Update all installed packages
# update.packages(ask = FALSE, checkBuilt = TRUE)
# 
# # Install required Bioconductor packages
# BiocManager::install(c("limma", "illuminaio", "lumi"))
# 
# # Install the igraph package from CRAN
# install.packages("igraph")
# 
# # Load the required libraries
# library(limma)
# library(GEOquery)
# library(illuminaio)
# library(lumi)
# 
# # Load the GEO dataset
# gse <- getGEO("GSE46448", GSEMatrix = TRUE)
# 
# 
# # Background correction
# bgCorrected <- neqc(exprsData)
# 
# exprsData <- exprs(gse[[1]])
# bgCorrected <- exprsData
# 
# 
# # Quantile normalization and log2 transformation
# normData <- normalizeBetweenArrays(bgCorrected, method = "quantile")
# logData <- log2(normData)
# 
# logData[,1] %>% hist
# 
# # Extract expression data
# exprsData <- exprs(gse[[1]])
# 
# exprsData[,1] %>% hist
# 
# exprsData[,1] %>% log2() %>% hist
# 
# 
# 
# 
# # Apply background correction using the neqc function from limma
# bgCorrected <- neqc(exprsData)
# 
# # Display dimensions of the expression data
# dim(exprsData)
# 
# # Display dimensions of the first few rows of the expression data
# head(exprsData) %>% dim
# 
# # Extract and process phenotype data
# # This part filters and selects specific samples based on conditions
# selectedSamples <- phenoData(gse[[1]]) %>% 
#   pData %>% 
#   # Uncomment the next line if you need to filter by a specific condition
#   # filter(`transfected with:ch1` != "none (non-infected control)") %>% 
#   select(c(title, geo_accession)) %>% 
#   filter(grepl("_0h_|_4h_", title)) %>% 
#   filter(grepl("siNS", title)) %>% 
#   arrange(title) %>% 
#   .$geo_accession
# 
# # Create a design matrix for the linear model
# # Adjust the factor levels as per your experimental design
# design <- model.matrix(~ factor(c("control", "control", "control", "control", "treated", "treated", "treated", "treated")))
# 
# # Fit the linear model using lmFit from limma
# # The columns are selected based on the filtered phenotype data
# fit <- lmFit(exprsData[, selectedSamples], design)
# 
# # Apply empirical Bayes smoothing
# fit <- eBayes(fit)
# 
# # Extract the results, applying adjustments for multiple testing
# results <- topTable(fit, adjust="BH", number=Inf)
# 
# results %>% 
#   filter(adj.P.Val < 0.02) %>% 
#   filter(abs(logFC) > log(2)) %>% 
#   rownames() %>% unique -> probe_id
# # Filter results based on adjusted p-value and log fold change
# sigResults <- subset(results, adj.P.Val < 0.02 & abs(logFC) > log2(2))
# 
# sigResults %>% dim
#   # rownames_to_column(., var = "id") %>% 
#   # filter(id == "ILMN_2215639")