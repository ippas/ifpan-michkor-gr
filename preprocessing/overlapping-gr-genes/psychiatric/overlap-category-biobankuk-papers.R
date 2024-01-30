source("preprocessing/functions/R/install-load-packages.R")



read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>%
  filter(category_type == "320_Origin_Categories") %>%
  group_by(category_name) %>%
  nest() %>%
  mutate(n = map_int(data, nrow)) %>%
  filter(n >= 10, n <= 300) %>% .$category_name

filter_phenotypes_by_category(
  genes_list = genes_phenotypes_PanUkBiobank,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  category_name = "Depression"
) -> depression_category_genes_list

analyze_gene_list_overlap(row_lists = depression_category_genes_list,
                          col_lists =  papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
                          reference_hgnc_vector = hgnc_symbols_vector_v110,
                          fdr_threshold = 0.01,
                          overlap_threshold = 3,
                          keep_original_data = FALSE) -> tmp




read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>%
  filter(category_type == "320_Origin_Categories") %>%
  group_by(category_name) %>%
  nest() %>%
  mutate(n = map_int(data, nrow)) %>%
  filter(n >= 10, n <= 300) %>% .$category_name 


###############################################################################

# Set the number of cores for parallel processing
num_cores <- 6  # Adjust this number based on your machine's capability
plan(multisession, workers = num_cores)

# Read and preprocess the data
categories_biobankuk <- read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>%
  filter(category_type == "320_Origin_Categories") %>%
  group_by(category_name) %>%
  nest() %>%
  mutate(n = map_int(data, nrow)) %>%
  filter(n >= 10, n <= 300) %>%
  pull(category_name) %>% .[1:6]


# Parallel processing
Sys.time() -> start_time
overlap_results <- future_map(set_names(categories_biobankuk), function(category) {
  filter_phenotypes_by_category(
    genes_list = genes_phenotypes_PanUkBiobank,
    path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
    category_name = category
  ) %>%
    analyze_gene_list_overlap(
      row_lists = .,
      col_lists = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
      reference_hgnc_vector = hgnc_symbols_vector_v110,
      fdr_threshold = 0.01,
      overlap_threshold = 3,
      keep_original_data = FALSE
    )
})
Sys.time() -> end_time

end_time - start_time

plan(NULL)
