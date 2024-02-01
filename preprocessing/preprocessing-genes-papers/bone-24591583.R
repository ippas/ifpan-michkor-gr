# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Update all installed packages
update.packages(ask = FALSE, checkBuilt = TRUE)

# Install required Bioconductor packages
BiocManager::install(c("limma", "illuminaio", "lumi"))

# Install the igraph package from CRAN
install.packages("igraph")

# Load the required libraries
library(limma)
library(GEOquery)
library(illuminaio)
library(lumi)

# Load the GEO dataset
gse <- getGEO("GSE46448", GSEMatrix = TRUE)


# Background correction
bgCorrected <- neqc(exprsData)

exprsData <- exprs(gse[[1]])
bgCorrected <- exprsData


# Quantile normalization and log2 transformation
normData <- normalizeBetweenArrays(bgCorrected, method = "quantile")
logData <- log2(normData)

logData[,1] %>% hist

# Extract expression data
exprsData <- exprs(gse[[1]])

exprsData[,1] %>% hist

exprsData[,1] %>% log2() %>% hist




# Apply background correction using the neqc function from limma
bgCorrected <- neqc(exprsData)

# Display dimensions of the expression data
dim(exprsData)

# Display dimensions of the first few rows of the expression data
head(exprsData) %>% dim

# Extract and process phenotype data
# This part filters and selects specific samples based on conditions
selectedSamples <- phenoData(gse[[1]]) %>% 
  pData %>% 
  # Uncomment the next line if you need to filter by a specific condition
  # filter(`transfected with:ch1` != "none (non-infected control)") %>% 
  select(c(title, geo_accession)) %>% 
  filter(grepl("_0h_|_4h_", title)) %>% 
  filter(grepl("siNS", title)) %>% 
  arrange(title) %>% 
  .$geo_accession

# Create a design matrix for the linear model
# Adjust the factor levels as per your experimental design
design <- model.matrix(~ factor(c("control", "control", "control", "control", "treated", "treated", "treated", "treated")))

# Fit the linear model using lmFit from limma
# The columns are selected based on the filtered phenotype data
fit <- lmFit(exprsData[, selectedSamples], design)

# Apply empirical Bayes smoothing
fit <- eBayes(fit)

# Extract the results, applying adjustments for multiple testing
results <- topTable(fit, adjust="BH", number=Inf)

results %>% 
  filter(adj.P.Val < 0.02) %>% 
  filter(abs(logFC) > log(2)) %>% dim

# Filter results based on adjusted p-value and log fold change
sigResults <- subset(results, adj.P.Val < 0.02 & abs(logFC) > log2(2))

sigResults %>% dim
  # rownames_to_column(., var = "id") %>% 
  # filter(id == "ILMN_2215639")

################################################################################
data_path <- "/home/mateusz/Downloads/1400522111_sd01.xlsx"

# Load the data, skipping the first two rows and using the third row as headers
read_excel(data_path, skip = 1, sheet = "Sheet1") -> data

data %>% as.data.frame() %>% 
  filter(`siNS Q.Value` < 0.05) %>% 
  mutate(abs_logFC = abs(`log FC by Dex in siNS`)) %>% 
  filter(abs_logFC > 1) %>%  
  select(c(PROBE_ID, SYMBOL, `log FC by Dex in siNS`, `siNS Q.Value`)) %>% 
  set_colnames(c("probe_ID", "gene_name", "log2ratio", "fdr")) %>% 
  mutate(regulation = ifelse(log2ratio > 0, "up", "down")) -> df_preprocessing
z

