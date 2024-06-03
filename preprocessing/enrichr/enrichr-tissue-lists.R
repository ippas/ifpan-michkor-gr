
databases2 <- c(databases, "CellMarker_2024", "GO_Cellular_Component_2023", "TRRUST_Transcription_Factors_2019", "Transcription_Factor_PPIs", "GO_Molecular_Function_2023", "Jensen_DISEASES", "Jensen_COMPARTMENTS",
                "Jensen_TISSUES", "Human_Phenotype_Ontology", "UK_Biobank_GWAS_v1", "GWAS_Catalog_2023")

multiple_tissue_enrichment_analysis <- function(data, databases, 
                                                fdr_threshold = 0.05, 
                                                overlap_threshold = 3, 
                                                regulation = "up", 
                                                rank_criterion = "log2ratio",
                                                size_list = 50){
  data %>% 
    filter(regulation == {{regulation}},
           rank_criterion == {{rank_criterion}},
           size_list == {{size_list}}) -> preprocessing_data
  
  preprocessing_data$tissue %>% unique() -> tissue_vector
  
  print(tissue_vector)
  
  lapply(tissue_vector, function(tissue){
    preprocessing_data %>% 
      filter(tissue == {{tissue}}) %>% 
      .$hgnc_symbol %>% 
      perform_enrichment_analysis(gene_list = .,
                                  databases = {{databases}},
                                  fdr_threshold = {{fdr_threshold}},
                                  overlap_threshold = {{overlap_threshold}})
  }) -> results
  
  names(results) <- tissue_vector
  
  return(results)
}

multiple_tissue_enrichment_analysis(data = detailed_tissue_master_df, databases = databases2) -> enrichr_detailed_tissue_up

multiple_tissue_enrichment_analysis(data = detailed_tissue_master_df, databases = databases2, regulation = "down") -> enrichr_detailed_tissue_down


# Process the data
processed_data <- enrichr_detailed_tissue_up %>% 
  lapply(function(df) {
    df %>% 
      arrange(Adjusted.P.value) %>% 
      filter(enrichr_database == "DisGeNET")
  }) %>% 
  lapply(head, 10) %>%  
  bind_rows(.id = "tissue") %>% 
  filter(tissue %in% c("adipose")) %>% 
  mutate(Term = factor(Term, levels = Term[order(Adjusted.P.value)])) %>% 
  mutate(n_genes = as.numeric(n_genes))

# Create the plot with reversed y-axis for each facet
ggplot(processed_data, aes(x = n_genes, y = reorder(Term, Adjusted.P.value))) +
  geom_col() +
  facet_wrap(~ tissue, ncol = 1, scales = "free_y") +
  labs(y = "Term", x = "Number of Genes") +
  theme_minimal() +
  scale_y_discrete(limits = function(x) rev(x))


# Process the data for all tissues
processed_data <- enrichr_detailed_tissue_up %>% 
  lapply(function(df) {
    df %>% 
      arrange(Adjusted.P.value) %>% 
      filter(enrichr_database == "DisGeNET")
  }) %>% 
  lapply(head, 10) %>%  
  bind_rows(.id = "tissue") %>% 
  mutate(Term = paste(tissue, Term, sep = "_")) %>% 
  mutate(Term = factor(Term, levels = Term[order(Adjusted.P.value)])) %>% 
  mutate(n_genes = as.numeric(n_genes))

processed_data %>% 
  write_tsv_xlsx(data = ., tsv_file = "results/google-drive/enrichr/enrichr-tissue-up-top10.tsv")


# Create the plot with reversed y-axis for each facet
ggplot(processed_data, aes(x = n_genes, y = reorder(Term, Adjusted.P.value))) +
  geom_col() +
  facet_wrap(~ tissue, ncol = 1, scales = "free_y") +
  labs(y = "Term", x = "Number of Genes") +
  theme_minimal() +
  scale_y_discrete(limits = function(x) rev(x))


processed_data <- enrichr_detailed_tissue_down %>% 
  lapply(function(df) {
    df %>% 
      arrange(Adjusted.P.value) %>% 
      filter(enrichr_database == "DisGeNET")
  }) %>% 
  lapply(head, 10) %>%  
  bind_rows(.id = "tissue") %>% 
  mutate(Term = paste(tissue, Term, sep = "_")) %>% 
  mutate(Term = factor(Term, levels = Term[order(Adjusted.P.value)])) %>% 
  mutate(n_genes = as.numeric(n_genes))

# Create the plot with reversed y-axis for each facet
ggplot(processed_data, aes(x = n_genes, y = reorder(Term, Adjusted.P.value))) +
  geom_col() +
  facet_wrap(~ tissue, ncol = 1, scales = "free_y") +
  labs(y = "Term", x = "Number of Genes") +
  theme_minimal() +
  scale_y_discrete(limits = function(x) rev(x))
  
read.csv("results/google-drive/enrichr/enrichr-term-manualy-filter.tsv", sep = "\t", header = TRUE) %>%
  group_by(tissue) %>%
  arrange(Adjusted.P.value) %>%
  slice(1) %>%
  ggplot(aes(x = n_genes, y = reorder(Term, Adjusted.P.value))) +
  geom_col() +
  # facet_wrap(~ tissue, ncol = 1, scales = "free_y") +
  labs(y = "Term", x = "Number of Genes") +
  theme_minimal() +
  scale_y_discrete(limits = function(x) rev(x)) +
  labs(title = "Enrichr, DisGeNET")
