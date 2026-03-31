require(rentrez)

get_paper_title <- function(pmid) {
  paper_info <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml", parsed = TRUE)
  title_path <- "//ArticleTitle"
  title <- xpathSApply(paper_info, title_path, xmlValue)
  if (length(title) > 0) {
    return(title)
  } else {
    return(NA)
  }
}


clusters_phenotypes_data$significant_data$df


clusters_papers_data$significant_data$df %>% 
  mutate(Var2 = recode(Var2, !!!clusters_mapping)) %>%
  select(c(Var1, Var2, chi2, p_value, fdr, number_overlap, overlap_genes)) %>% 
  mutate(Var1 = str_replace_all(Var1, "michkor-cells", "pmid:28381250")) %>% 
  mutate(Var1 = str_replace_all(Var1, "_NA", "")) %>% 
  set_colnames(c("GR_dependent_list", "transcriptional_pattern", "chi2_value", "p_value", "fdr", "number_of_genes", "overlapping_genes")) %>% 
  write.table("results/tables/Supplementary_table2.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


clusters_phenotypes_data$significant_data$df %>% 
  mutate(Var2 = recode(Var2, !!!clusters_mapping)) %>%
  select(c(Var1, Var2, chi2, p_value, fdr, number_overlap, overlap_genes)) %>% 
  # mutate(Var1 = str_replace_all(Var1, "michkor-cells", "pmid:28381250")) %>% 
  set_colnames(c("phenotype", "transcriptional_pattern", "chi2_value", "p_value", "fdr", "number_of_genes", "overlapping_genes")) %>% 
  write.table("results/tables/Supplementary_table3.tsv", quote = F, col.names = T, row.names = F, sep = "\t")


clusters_metabolism_data$significant_data$df %>% 
  mutate(Var2 = recode(Var2, !!!clusters_mapping)) %>%
  select(c(Var1, Var2, chi2, p_value, fdr, number_overlap, overlap_genes)) %>%
  # mutate(Var1 = str_replace_all(Var1, "michkor-cells", "pmid:28381250")) %>% 
  set_colnames(c("metabolite_trait", "transcriptional_pattern", "chi2_value", "p_value", "fdr", "number_of_genes", "overlapping_genes")) %>% 
  write.table("results/tables/Supplementary_table4.tsv", quote = F, col.names = T, row.names = F, sep = "\t")




# Define the function to fetch paper title
get_paper_title <- function(pmid) {
  paper_info <- entrez_fetch(db = "pubmed", id = pmid, rettype = "xml", parsed = TRUE)
  title_path <- "//ArticleTitle"
  titles <- xpathSApply(paper_info, title_path, xmlValue)
  if (length(titles) > 0) {
    return(titles[1])
  } else {
    return(NA)
  }
}

# Your existing data preprocessing pipeline
papers_data_preprocessed <- papers_data_preprocessing %>%
  select(source, tissue, cell, treatment, dose, time, species) %>% 
  mutate(tissue_cell = paste0(tissue, ", ", cell)) %>%
  mutate(tissue_cell = str_replace_all(tissue_cell, ", NA", "")) %>% 
  mutate(tissue_cell = str_replace_all(tissue_cell, "NA, ", "")) %>% 
  filter(!grepl("marpiech", source)) %>% 
  filter(source != "pmid:NA") %>% 
  unique() %>%
  mutate(source = str_replace_all(source, "michkor-cells", "pmid:28381250"))

unique_pmids <- distinct(papers_data_preprocessed, source) %>%
  filter(str_detect(source, "^pmid:")) %>%
  mutate(pmid = str_replace(source, "pmid:", "")) %>%
  pull(pmid)


# Fetch titles for these PMIDs
titles <- sapply(unique_pmids, get_paper_title)

# Create a lookup table for PMIDs and their titles
pmid_title_lookup <- data.frame(pmid = unique_pmids, title = titles)

# Merge the titles back into the original dataset
papers_data_final <- papers_data_preprocessed %>%
  mutate(pmid = str_replace(source, "pmid:", "")) %>%
  left_join(pmid_title_lookup, by = "pmid") %>%
  select(-pmid)

# View the final dataset
papers_data_final %>%
  mutate(source = str_replace_all(source, "pmid", "PMID")) %>% 
  select(-c(cell, tissue)) %>% 
  select(c(title, source, tissue_cell, species, treatment, dose, time)) %>% 
  as.data.frame() %>% 
  select(-time) %>% 
  unique() %>% 
   write.table("results/tables/supplementary_table1a.tsv", col.names = T, quote = F, row.names = F, sep = "\t")
