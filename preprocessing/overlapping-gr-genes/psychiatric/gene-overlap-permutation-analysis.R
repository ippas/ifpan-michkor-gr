source("preprocessing/functions/R/install-load-packages.R")

source("preprocessing/functions/R/gene-overlap-permutation-analysis/generate-random-gene-sets.R")
source("preprocessing/functions/R/gene-overlap-permutation-analysis/analyze-random-gene-sets.R")
source("preprocessing/functions/R/gene-overlap-permutation-analysis/perform-overlap-permutation-analysis-parallel.R")


Sys.time() -> start_time
permutation_results3 <- perform_overlap_permutation_analysis_parallel(
  permutations = 4,  # Number of permutations
  seed = 123,          # Initial seed for reproducibility
  num_cores = 4,       # Number of cores to use
  # Additional arguments for analyze_random_gene_sets
  reference_hgnc_vector = hgnc_symbols_vector_v110,
  size_reference_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  comparison_gene_list = phenotypes_biobank_genes_list,
  overlap_threshold = 3,
  fdr_threshold = 0.01
)

Sys.time() -> end_time

start_time - end_time
