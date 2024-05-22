################################################################################
# create empty list to store differen stores blocked genes lists
gr_database_blocked_gene_lists <- list()

# prepare block with marpiech_tissues_dex genes list
gr_database_blocked_gene_lists$marpiech_cluster_dex <- marpiech_data_preprocessing %>% 
  filter((source %in% c("marpiech_clusters_dex"))) %>% 
  mutate(label = paste0("cluster_", gene_list_number)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})

# prepare block with marpiech_tissues_dex genes list
gr_database_blocked_gene_lists$marpiech_tissue_dex <- papers_data_preprocessing %>% 
  filter(gene_list_index == "all_significant_genes_marpiech") %>%
  select(c(hgnc_symbol, tissue)) %>% 
  split(.$tissue, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol}) 

gr_database_blocked_gene_lists$non_GR_dependent_gene_lists <- papers_data_preprocessing %>% 
  filter((treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})

################################################################################
gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_time_dose_regulation$metadata <- 
  papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>%
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, time, dose, regulation, sep = "_")) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_time_dose_regulation$gene_lists <- papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, time, dose, regulation, sep = "_")) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})

################################################################################
gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$metadata <- 
  papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$gene_lists <- papers_data_preprocessing %>% 
  # filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  # filter(!(source %in% c("marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  # filter(source != "marpiech_tissues") %>% 
  mutate(label = paste(source, tissue, cell, treatment, treatment_type, regulation, sep = "_")) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})



# prepare simple_tissue column

gr_database_blocked_gene_lists$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$metadata %>% 
  # select(c(tissue, cell)) %>% unique %>% as.data.frame() %>% 
  mutate(tissue = ifelse(is.na(tissue), cell, tissue)) %>% 
  # select(tissue) %>%
  categorize_tissue(tissue_column = "tissue") %>%
  select(c(source, tissue, cell, system, simple_tissue)) %>% 
  select(source, simple_tissue) %>% 
  unique() %>%
  .$simple_tissue %>% table() %>% 
  as.data.frame() %>% 
  set_colnames(c("tissue", "freq")) %>% 
  arrange(desc(freq)) %>% 
  ggplot(., aes(x = reorder(tissue, -freq), y = freq)) + 
  geom_bar(stat = "identity") +
  geom_hline(yintercept = 3, color = "red", linetype = "dashed") +
  theme(axis.text.y = element_text(angle = 90, vjust = 0.5, hjust=1))



gr_database_blocked_genes_list$GR_dependent_gene_lists_pmid_tissue_cell_treatment_type_regulation$all_list %>% 
  lapply(., function(x){length(x)}) %>% unname() %>% unlist %>% .[{. < 100}] %>% hist



papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("pmid:NA", "marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio)) 



gr_database_blocked_genes_list %>% 
  .$GR_dependent_gene_lists_pmid_tissue_cell_treatment_time_dose_regulation %>% 
  .$metadata 


papers_data_preprocessing %>% 
  filter(!(treatment %in% c("TNF", "LPS", "vehicle-ethanol", "TNFalpha"))) %>% 
  filter(!(source %in% c("pmid:NA", "marpiech_tissues_dex",  "pmid:37217509"))) %>% 
  select(-c(index, gene_name, ensembl_gene_id, ensembl_transcript_id, refseq_mrna_id, hgnc_symbol, alias, info, fdr, log2ratio))

gr_database_blocked_genes_list$GR_dependent_gene_lists_pmid_tissue_cell_treatment_regulation

papers_data_preprocessing %>% 

  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "frequency")) %>% 
  arrange(desc(frequency)) %>% 
  filter(frequency == 10) %>% 
  .$hgnc_symbol %>% 
  as.character() -> unique_genes 

  

papers_data_preprocessing %>% 
  # select(c(source, hgnc_symbol)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})
  
  
papers_data_preprocessing %>% 
  filter(!(source %in% c("pmid:NA", "marpiech_tissues_dex"))) %>% 
  filter(hgnc_symbol %in% unique_genes)  ->  tmp

split(tmp$hgnc_symbol, tmp$label) %>% lapply(., unique) %>% lapply(., length) %>% 
  as_tibble() %>% as.data.frame() %>% t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "list_name") %>% 
  rename(frequency = V1) %>% #filter(frequency > 200)
  .$frequency %>% as.numeric() %>% hist
 
