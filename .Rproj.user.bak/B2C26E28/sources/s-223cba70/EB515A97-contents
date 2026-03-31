# ================================================
# Define file paths
# ================================================
mdd2025_gwas_file <- "data/pgc-genes/pgc_mdd2025_no23andMe_eur.tsv"
gr_signatures_intersection_file <- "data/pgc-genes/intesection_results_all_gr_signatures.tsv"

# ================================================
# Load data
# ================================================
mdd2025_no23andMe_eur_gwas <- read.table(mdd2025_gwas_file, header = TRUE)
gr_signatures_intersection_data <- read.table(gr_signatures_intersection_file, header = TRUE)
protein_coding_intersection_data_GRCh37 <- read.table("data/pgc-genes/intesection_results_all_protein_coding.tsv", header = TRUE)


protein_coding_intersection_data_GRCh37 <- read.table("data/pgc-genes/intersection_results_all_protein_coding_10000.tsv", header = TRUE)


# ================================================
# Prepare preprocessed dataset with GR metadata
# ================================================
mdd2025_no23andMe_eur_gwas_preprocessing <- mdd2025_no23andMe_eur_gwas %>%
  mutate(
    gene_symbol = "None",
    signature_derivation = "depression_PMID:39814019",
    signature_name = "depression_PMID:39814019",
    range_plus_minus = "None"
  ) %>%
  select(
    rsID, chromosome, position, pvalue,
    gene_symbol, signature_name, signature_derivation, range_plus_minus
  )

mdd2025_no23andMe_eur_gwas_preprocessing %>% head
gr_signatures_intersection_data %>% head
protein_coding_intersection_data_GRCh37 %>% head


# ================================================
# Prepare depression gene lists from EnrichR
# ================================================

# 1. Run enrichment for selected databases
genes_DisGeNET <- enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "DisGeNET"
)

genes_GWAS_Catalog_2023 <- enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "GWAS_Catalog_2023"
)

genes_PhenGenI_Association_2021 <- enrichr_download_gene_lists(
  gene_list = hgnc_symbols_vector_v110,
  database = "PhenGenI_Association_2021"
)

# 2. Filter and reshape enrichment results for depression-related terms
rbind(
  {genes_DisGeNET %>% 
      filter(n_all_genes > 100) %>% 
      filter(grepl("depression", term, ignore.case = TRUE))},
  {genes_GWAS_Catalog_2023 %>% 
      filter(n_all_genes > 100) %>% 
      filter(grepl("depression", term, ignore.case = TRUE))}) %>% 
  mutate(label = paste(term, overlap, database, sep = "_")) %>% 
  mutate(label = str_replace_all(label, " ", "_")) %>%  
  ungroup %>% 
  as.data.frame() %>% 
  select(label, combined_genes) %>% 
  mutate(gene_symbol = strsplit(combined_genes, ";\\s*")) %>%
  select(-combined_genes) %>%
  unnest(gene_symbol) -> enrichr_depression_genes


enrichr_depression_genes %>% 
  left_join(., protein_coding_intersection_data_GRCh37, by = "gene_symbol") %>% 
  mutate(
    signature_derivation = label,
    signature_name = label,
    range_plus_minus = "100000"
  ) %>%
  select(
    rsID, chromosome, position, pvalue,
    gene_symbol, signature_name, signature_derivation, range_plus_minus
  ) %>% 
  as.data.frame() -> enrichr_depression_genes_preprocessing

enrichr_depression_genes_preprocessing %>% head


# ================================================
# Prepare Slezak gene signatures
# ================================================
slezak_signatures_list %>%
  tibble::enframe(name = "signature_name", value = "gene_symbol") %>%
  tidyr::unnest(gene_symbol) %>% 
  left_join(., protein_coding_intersection_data_GRCh37, by = "gene_symbol") %>% 
  filter(is.na(rsID)) %>% 
  as.data.frame() %>% select(signature_name, gene_symbol) %>% 
  .$gene_symbol %>% unique %>% length()

slezak_signatures_list %>%
  tibble::enframe(name = "signature_name", value = "gene_symbol") %>%
  tidyr::unnest(gene_symbol) %>% 
  group_by(signature_name) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~nrow(.x))) %>% 
  unnest(c(data, n_genes)) %>% 
  ungroup() %>% 
  mutate(signature_name = paste0(signature_name, "_n", n_genes)) %>% 
  left_join(., protein_coding_intersection_data_GRCh37, by = "gene_symbol") %>% 
  as.data.frame() %>% 
  filter(!is.na(rsID)) %>% 
  mutate(
    signature_derivation = "slezak_signature",
    range_plus_minus = "100000"
  ) %>% 
  select(
    rsID, chromosome, position, pvalue,
    gene_symbol, signature_name, signature_derivation, range_plus_minus
  ) %>% 
  as.data.frame() -> slezak_gene_signatures_preprocessing



rm(
  mdd2025_gwas_file,
  gr_signatures_intersection_file,
  mdd2025_no23andMe_eur_gwas,
  gr_signatures_intersection_data,
  protein_coding_intersection_data_GRCh37,
  mdd2025_no23andMe_eur_gwas_preprocessing,
  genes_DisGeNET,
  genes_GWAS_Catalog_2023,
  genes_PhenGenI_Association_2021,
  enrichr_depression_genes,
  enrichr_depression_genes_preprocessing,
  slezak_gene_signatures_preprocessing
)
