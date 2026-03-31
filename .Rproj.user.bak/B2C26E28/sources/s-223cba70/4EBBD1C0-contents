



# # Read data
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/geneBase_060723.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_100923.tsv", sep = "\t") 
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_v5_181023.tsv", sep = "\t")
# gr_gene_database_raw <- read.csv("data/overlapping-gr-genes/publikacje_gr_v6_181023.tsv", sep = "\t")
# 
# # Extract specific keys and values from the "Info" column
# # gr_gene_database_raw %>% 
# #   extract_keys_values("Info", c("tissue", "cell","environment", "treatment", "dose", "time", "log2ratio")) -> gr_gene_database
# 
# gr_gene_database <- gr_gene_database_raw %>% 
#   extract_keys_values("info", c("tissue", "cell", "species", "environment", "treatment", "dose", "time", "log2ratio", "fdr", "statistical_method", "treatment_type", "regulation"))
# 
# # Filter rows where source is "michkor-cells", statistical_method is "ttest", and fdr is less than 0.05
# tmp_michkor <- gr_gene_database %>% 
#   filter(
#     source == "michkor-cells", 
#     statistical_method == "ttest", 
#     fdr < 0.05
#   )
# 
# # Create a new dataframe combining the relevant rows and selecting the relevant columns
# gr_gene_database <- gr_gene_database %>% 
#   filter(source != "michkor-cells") %>% 
#   bind_rows(tmp_michkor) 
# 
# gr_gene_database <- gr_gene_database %>% 
#   filter(gene_list_number != 800)

# gr_gene_database %>% 
#   filter(grepl("marpiech", gene_list_index)) %>% 
#   filter(gene_list_number != 800) %>% 
#   filter(gene_list_index != "marpiech_tissues_dex_up") %>% 
#   filter(gene_list_index != "marpiech_tissues_dex_down") %>% 
#   select(-index) %>%
#   unique() %>% 
#   select(c(gene_list_number, gene_name, hgnc_symbol)) %>% 
#   .$hgnc_symbol %>% 
#   table %>% 
#   as.data.frame() %>% 
#   arrange(Freq) %>%
#   set_colnames(c("hgnc_symbol", "freq")) %>% 
#   filter(freq > 1) %>% 
#   as.data.frame() %>% .$hgnc_symbol %>% as.character() -> hgnc_to_remove

# ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# # Retrieve the HGNC symbols for protein-coding genes
# hgnc_symbols <- biomaRt::getBM(
#   attributes = c("hgnc_symbol", "gene_biotype"),
#   filters = "biotype",
#   values = "protein_coding",
#   mart = ensembl
# )
# 
# hgnc_symbols %>% .$hgnc_symbol -> hgnc_symbols_vector

# Define the directory path
dir_path <- "data/prs-models-pan-biobank-uk/"

# List all files in the specified directory and store them in a data frame
file_list <- list.files(path = dir_path) %>% as.data.frame() %>% set_colnames("file")

# Initialize a list to store genes information from each file
genes_info_list <- list()

# Loop over each file in the directory
for (file in file_list$file) {
  # Construct the full file path
  file_path <- file.path(dir_path, file)
  
  # Load the YAML file and extract the genes information
  yaml_content <- yaml.load_file(file_path)
  genes_info <- yaml_content$description$genes
  
  # Add the extracted genes information to the list, using the file name as the list element name
  genes_info_list[[file]] <- genes_info
}

# Filter the list elements based on the specified genes
filtered_genes_info_list <- lapply(genes_info_list, function(vector) vector[vector %in% {genes_list %>% unlist() %>%  as.vector() %>% unique()}])

# Further filter the list to keep only elements where the length of the vector is greater than 2
filtered_genes_info_list <- filtered_genes_info_list[sapply(filtered_genes_info_list, function(vector) length(vector) > 2)]

filtered_genes_info_list <- filtered_genes_info_list[grepl("1e-08", names(filtered_genes_info_list))]

# Extract the names of the filtered phenotypes
filtered_phenotypes <- names(filtered_genes_info_list)

# Extract the genes information corresponding to the filtered phenotypes
filtered_phenotypes <- genes_info_list[names(genes_info_list) %in% filtered_phenotypes]

# Display the filtered phenotypes
print(filtered_phenotypes)

# Preprocess the gene database
gr_gene_database_preproccesing <- gr_gene_database %>%
  filter(!grepl("omicspred_metabolon", source)) %>%
  filter(!grepl("omicspred_nithingale_", source)) %>% 
  filter(!hgnc_symbol %in% hgnc_to_remove) %>% 
  mutate(gene_list_index = gsub("marpiech_.*", "marpiech_tissues_dex", gene_list_index)) %>% 
  mutate(source = ifelse(gene_list_index == "marpiech_tissues_dex", "marpiech_tissues_dex", source)) %>% 
  # Clean up and modify columns
  mutate(
    source = str_replace_all(source, "PMID: ", ""),
    gene_name = tolower(gene_name),
    label = ifelse(
      source == "marpiech_tissues_dex",
      paste0(source, "_", gene_list_number),
      paste0(source, "_", tissue, "_", cell)
    )
  ) %>% 
  # Filter rows based on multiple conditions
  filter(
    hgnc_symbol %in% combined_genes$external_gene_name,
    # gene_name %in% {
    #   filtered_biomart_by_go %>%
    #     pull(gene_name) %>%
    #     tolower() %>%
    #     unique()
    # },
    label != "34362910_NA_NA"
  ) %>% 
  mutate(label = paste0(source, "_", tissue, "_", cell)) %>%
  mutate(label = ifelse(source == "marpiech_tissues_dex", paste0(source, "_", gene_list_number), label))


# Overlap of genes 
gr_gene_database_preproccesing %>% filter(source == "marpiech_tissues_dex")
  
split(gr_gene_database_preproccesing$hgnc_symbol, gr_gene_database_preproccesing$label) -> gr_gene_list

results <- perform_chi2_tests(c(gr_gene_list, genes_list), hgnc_symbols_vector)

results$chi2_value_matrix

pheatmap(
  log2(results$chi2_value_matrix + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE
)

pheatmap(
  log2(results$chi2_value_matrix[19:ncol(results$chi2_value_matrix), 1:18] + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE
)

#
results <- perform_chi2_tests(c(genes_list, filtered_phenotypes), hgnc_symbols_vector)

results$p_value_matrix %>% .[1:11, 12:ncol(results$p_value_matrix)] %>% reshape2::melt() %>% 
  arrange(value) %>%
  # filter(Var1 %in% c("nervous_gene_separate_29180230", "nervous_gene_weak_29180230")) %>% 
  mutate(fdr = p.adjust(value, "fdr")) %>% 
  filter(fdr < 0.01)



# results$p_value_matrix %>% .[1:10, 1:10]

# Create heatmaps
# Visualize the p-value matrix using heatmaps
pheatmap(
  results$p_value_matrix[1:5, 6:ncol(results$p_value_matrix)],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE
)

##########################
# Your filtering code
filtered_data <- results$p_value_matrix %>% 
  .[1:11, 12:ncol(results$p_value_matrix)] %>% 
  reshape2::melt() %>% 
  arrange(value) %>%
  mutate(fdr = p.adjust(value, "fdr")) %>% 
  filter(fdr < 0.01)


unique_rows <- unique(filtered_data$Var1)
unique_cols <- unique(filtered_data$Var2)

unique_cols <- setdiff(unique_cols, unique_rows)

filtered_matrix_subset <- results$chi2_value_matrix[ unique_cols, unique_rows]



# Convert the filtered data back to a matrix
filtered_matrix <- reshape2::acast(filtered_data, Var1 ~ Var2, value.var = "fdr")


pheatmap(
  log2(filtered_matrix_subset + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE
)

#########################33
# Visualize the chi-square value matrix using heatmaps
pheatmap(
  log2(results$chi2_value_matrix[1:9, 11:ncol(results$chi2_value_matrix)] + 1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE
)

# Visualize the number overlap matrix using heatmaps
pheatmap(
  results$number_overlap_matrix[1:39, 1:39],
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE
)


# Perform chi-square tests on combined lists and a vector of symbols
results <- perform_chi2_tests(
  c(gene_list, gr_gene_list, filtered_phenotypes), 
  hgnc_symbols_vector
)

# Visualize the p-value matrix using a heatmap
pheatmap(results$p_value_matrix[1:39, 1:39], 
           cluster_rows = TRUE, 
           cluster_cols = TRUE, 
           display_numbers = TRUE,
           main = "P-value Matrix"
)

# Visualize the log-transformed chi-square value matrix using a heatmap
pheatmap(
  log2(results$chi2_value_matrix[1:39, 1:39] + 1), 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = TRUE,
  main = "Log-transformed Chi-square Value Matrix"
)

# Visualize the number overlap matrix using a heatmap
pheatmap(
  results$number_overlap_matrix[1:39, 1:39], 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = TRUE,
  main = "Number Overlap Matrix"
)

pheatmap(
  log2(results$chi2_value_matrix + 1), 
  cluster_rows = TRUE, 
  cluster_cols = TRUE, 
  display_numbers = TRUE,
  main = "Log-transformed Chi-square Value Matrix"
)
