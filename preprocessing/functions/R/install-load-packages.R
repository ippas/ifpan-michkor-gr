# Function to install and load packages
install_and_load_packages <- function(packages) {
  new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)
  sapply(packages, require, character.only = TRUE)
}

# List of packages to be installed and loaded
required_packages <- c("readxl", "writexl", "magrittr", "tidyverse", "parallel", "future.apply", "future", "yaml", "furrr", "purrr")

# Install and load the packages
install_and_load_packages(required_packages)

# Load functions

# preprocessing gene papers functions
source("preprocessing/functions/R/preprocessing-gene-papers/download-and-read-geo-file.R")
source("preprocessing/functions/R/preprocessing-gene-papers/select_columns_by_pattern.R")
source("preprocessing/functions/R/preprocessing-gene-papers/perform_DESeq2_analysis.R")

# permutation functions
source("preprocessing/functions/R/overlapping-genes-analysis/gene-overlap-permutation-analysis/analyze-random-gene-sets.R")
source("preprocessing/functions/R/overlapping-genes-analysis/gene-overlap-permutation-analysis/generate-random-gene-sets.R")
source("preprocessing/functions/R/overlapping-genes-analysis/gene-overlap-permutation-analysis/perform-overlap-permutation-analysis-parallel.R")
source("preprocessing/functions/R/overlapping-genes-analysis/gene-overlap-permutation-analysis/perform-overlap-permutation-analysis.R")

# overlapping gene analysis
source("preprocessing/functions/R/overlapping-genes-analysis/read_genes_from_phenotype_models.R")
source("preprocessing/functions/R/overlapping-genes-analysis/filter_phenotypes_by_category.R")
source("preprocessing/functions/R/overlapping-genes-analysis/analyze-gene-list-overlap.R")


# Load in custom functions related to preprocessing
source("preprocessing/functions/R/gr-database-functions.R")
source("preprocessing/functions/R/draw-custom-heatmap.R")
source("preprocessing/functions/R/processing-overlap-results.R")
source("preprocessing/functions/R/manual-filter-overlap-results.R")
source("preprocessing/functions/R/gene-paper-preprocessing-functions.R")
source("preprocessing/functions/R/perform-chi2-tests.R")

