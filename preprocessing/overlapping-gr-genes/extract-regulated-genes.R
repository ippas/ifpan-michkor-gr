# Load required libraries
install.packages(c("wordcloud", "tm", "VennDiagram"))
library(dplyr)
library(wordcloud)
library(VennDiagram)
library(tm)

# Define functions to compute unique and intersect elements for a list of vectors
# Function to find the intersection of genes by group
# Extracting genes with different regulation levels
extract_regulated_genes <- function(database, column_name, regulation_level = NULL, frequency) {
  
  # Start the pipeline
  result <- database %>%
    select(c("source", "hgnc_symbol", column_name)) %>%
    na.omit() %>%
    rename(gene_name = hgnc_symbol)
  
  # Apply the filter if regulation_level is not NULL
  if (!is.null(regulation_level)) {
    result <- result %>%
      filter(get(column_name) %in% c(regulation_level))
  }
  
  # Continue the pipeline
  result %>%
    distinct(source, gene_name, .keep_all = TRUE) %>%
    pull(gene_name) %>%
    table() %>%
    as.data.frame() %>%
    arrange(desc(Freq)) %>%
    set_colnames(c("gene_name", "freq")) %>%
    filter(freq >= frequency) %>%
    pull(gene_name) %>%
    as.vector()
}


# Extract and filter gene data
gr_gene_database %>%
  select(c("source", "hgnc_symbol", "cell")) %>%
  na.omit() %>%
  rename(gene_name = hgnc_symbol) %>%
  filter(cell %in% c("NPSCs", "glia", "mglia", "astrocyte", "neuron", "astrocyte", "OPC", "oligodendrocytes")) %>%
  distinct(source, gene_name, .keep_all = TRUE) %>%
  pull(gene_name) %>%
  table() %>%
  as.data.frame() %>%
  arrange(desc(Freq)) %>% 
  filter(Freq >= 3) %>% 
  set_colnames(c("gene_name", "freq")) %>% 
  pull(gene_name) %>%
  as.vector() -> nervous_gene


# Extract genes based on their level of upregulation from the `gr_gene_database`.
# Each extraction looks at a specified `regulation_level`, i.e., "up" or "down", 
# and checks the `frequency` to determine the strength of the regulation.

# Strong upregulated genes with a frequency of 7
upregulation_genes_strong <- extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = "up",
  frequency = 7
)

# Moderately upregulated genes with a frequency of 6
upregulation_genes_middle <-  extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = "up",
  frequency = 6
)

# Weakly upregulated genes with a frequency of 5
upregulation_genes_weak <- extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = "up",
  frequency = 5
)

# Genes with downregulation at a frequency of 4
downregulation_genes <-  extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = "down",
  frequency = 4
)

# It seems you have some redundant or potentially unused lines below, 
# so they are commented with potential recommendations.

# Extraction with a NULL regulation_level; might represent a general extraction at a frequency of 7
# master_gr_genes <- extract_regulated_genes(gr_gene_database, "regulation", NULL, 7)

# Another extraction, possibly redundant
# extract_regulated_genes(gr_gene_database, "regulation", "up", 7)

# Strong master genes with a frequency of 8
master_gr_genes_strong <- extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = NULL,
  frequency = 8
)

# Moderate master genes with a frequency of 7
master_gr_genes_middle <- extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = NULL,
  frequency = 7
)

# Weak master genes with a frequency of 6
master_gr_genes_weak <- extract_regulated_genes(
  gr_gene_database,
  column_name = "regulation",
  regulation_level = NULL,
  frequency = 6
)

# Creating a list to group various sets of genes. 
# It appears `nervous_gene` was not defined in the provided code, so ensure you have defined it elsewhere.
genes_list <- list(
  master_downregulation = downregulation_genes,
  master_upregulation_weak = upregulation_genes_weak,
  master_upregulation_middle = upregulation_genes_middle,
  master_upregulation_strong = upregulation_genes_strong,
  master_nervous = nervous_gene, 
  master_gr_strong = master_gr_genes_strong,
  master_gr_middle = master_gr_genes_middle,
  master_gr_weak = master_gr_genes_weak
)


genes_list %>%
  enframe(name = "source", value = "hgnc_symbol") %>%
  unnest(cols = hgnc_symbol) %>% 
  as.data.frame() %>% 
  write.table(file = "data/supplement-genes-papers/master-gr-gene-lists.tsv",
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)


# data fraom review
# CTFG changed to CCN2 
# problem z LHFP wiele oznaczeń -> może jako Lhfp6
nervous_gene_down_29180230 <- c("Cxxc5", "Litaf", "Calm2", "Ccnd1", "Cyp7b1", "Gab1", "Gap43", "Irf1", "Jun", "Nr3c1", "Osbpl3", "Plscr1", "Sall2", "Scamp2", "Sox2", "Sox4","Sox9", "Tgfbr1", "Tle4", "Tmem109", "Vps37b", "Wnt7a", "Zfp36l1", "Cklf", "Sema6d")

nervous_gene_strong_29180230 <- c("Ddit4", "Errfi1", "Klf9", "Bcl6", "Fkbp5", "Nfkbia", "Pdk4", "Mt2A", "Adcy9", "Cxxc5", "Dusp1", "Eva1a", "Litaf", "Nedd9", "Rhob", "Sgk1", "Sult1a1", "Tiparp", "Gpd1", "Ccn2", "Plekhf1")

nervous_gene_weak_29180230 <- c(
  "Ddit4", "Errfi1", "Klf9", "Bcl6", "Fkbp5", "Nfkbia", "Pdk4", "Mt2A", "Adcy9", "Cxxc5", 
  "Dusp1", "Eva1a", "Litaf", "Nedd9", "Rhob", "Sgk1", "Sult1a1", "Tiparp", "Gpd1", "Ccn2", "Plekhf1",
  "Aldoc", "Arhgef3", "Arl4d", "Bcl6b", "Cables1", "Calm2", "Ccnd1", "Cdkn1a", "Cdo1",
  "Chst1", "Cyp7b1", "Ehd3", "Fzd1", "Gab1", "Gap43", "Gjb6", "Hepacam", "Id1", "Irf1",
  "Jun", "Klf15", "Lhfp", "Lyve1", "Mertk", "Mgst1", "Mical2", "Myh2", "Ndrg2", "Npy1r",
  "Nr3c1", "Nudt9", "Osbpl3", "Pim3", "Plscr1", "Prr5", "Rasl11b", "Rdx", "Rhou", "Sall2",
  "Scamp2", "Sdc4", "Sesn1", "Slc25a33", "Sox2", "Sox4", "Sox9", "Spsb1", "Svil", "Tgfbr1",
  "Thra", "Tle4", "Tmem109", "Tob2", "Tsc22d3", "Vps37b", "Wipf3", "Wnt16", "Wnt7a",
  "Il6R", "Dgkz", "Mtmr2", "Zfp36l1", "Azin1", "Cklf", "Ppp5c", "Sema6d", "Tle3")

# ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
# mapping_data <- biomaRt::getBM(attributes=c('hgnc_symbol', 'external_gene_name'), mart=ensembl)

# nervous_gene_weak_29180230 %>% 
#   as.data.frame() %>% 
#   set_colnames("gene_name") %>% 
#   mutate(regulation = ifelse(gene_name %in% nervous_gene_down_29180230, "down", "up")) %>% 
#   mutate(article_source = "https://pubmed.ncbi.nlm.nih.gov/29180230/") %>% 
#   mutate(pmid = "nervous_gene_weak_29180230") %>%  
#   left_join(mapping_data, by = c("gene_name" = "external_gene_name"))
# 
nervous_gene_weak_29180230 %>%
  as.data.frame() %>%
  set_colnames("gene_name") %>%
  mutate(gene_name = toupper(gene_name)) %>%  # Convert to uppercase
  mutate(regulation = ifelse(gene_name %in% nervous_gene_down_29180230, "down", "up")) %>%
  mutate(article_source = "https://pubmed.ncbi.nlm.nih.gov/29180230/") %>%
  mutate(pmid = "nervous_gene_weak_29180230") %>%
  left_join(mapping_data, by = c("gene_name" = "external_gene_name"))

nervous_gene_strong_29180230 %>%
  as.data.frame() %>%
  set_colnames("gene_name") %>%
  mutate(regulation = ifelse(gene_name %in% nervous_gene_down_29180230, "down", "up")) %>%
  mutate(article_source = "https://pubmed.ncbi.nlm.nih.gov/29180230/") %>%
  mutate(pmid = "nervous_gene_strong_29180230") %>%
  mutate(ensembl_id = NA, log2ratio = NA) %>% 
  select(c(article_source, pmid, gene_name, ensembl_id, log2ratio, regulation)) -> nervous_gene_strong_29180230_df

nervous_gene_strong_29180230 %>%
  as.data.frame() %>%
  set_colnames("gene_name") %>%
  mutate(regulation = ifelse(gene_name %in% nervous_gene_down_29180230, "down", "up")) %>%
  mutate(article_source = "https://pubmed.ncbi.nlm.nih.gov/29180230/") %>%
  mutate(pmid = "nervous_gene_strong_29180230") %>%
  mutate(ensembl_id = NA, log2ratio = NA) %>% 
  select(c(article_source, pmid, gene_name, ensembl_id, log2ratio, regulation)) -> nervous_gene_weak_29180230_df


rbind(nervous_gene_strong_29180230_df, nervous_gene_weak_29180230_df) %>% 
  write.table(., file = "data/supplement-genes-papers/review-nervous-29180230/nervous-gene-2918230.tsv",
              quote = FALSE,
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
# add nervous list
genes_list$nervous_gene_weak_29180230 <- nervous_gene_weak_29180230 %>%
  as.data.frame() %>%
  set_colnames("gene_name") %>%
  mutate(gene_name = toupper(gene_name)) %>%  # Convert to uppercase
  mutate(regulation = ifelse(gene_name %in% nervous_gene_down_29180230, "down", "up")) %>%
  mutate(article_source = "https://pubmed.ncbi.nlm.nih.gov/29180230/") %>%
  mutate(pmid = "nervous_gene_weak_29180230") %>%
  left_join(mapping_data, by = c("gene_name" = "external_gene_name")) %>% .$hgnc_symbol

genes_list$nervous_gene_strong_29180230 <- nervous_gene_strong_29180230 %>%
  as.data.frame() %>%
  set_colnames("gene_name") %>%
  mutate(gene_name = toupper(gene_name)) %>%  # Convert to uppercase
  mutate(regulation = ifelse(gene_name %in% nervous_gene_down_29180230, "down", "up")) %>%
  mutate(article_source = "https://pubmed.ncbi.nlm.nih.gov/29180230/") %>%
  mutate(pmid = "nervous_gene_strong_29180230") %>%
  left_join(mapping_data, by = c("gene_name" = "external_gene_name")) %>% .$hgnc_symbol

genes_list$nervous_gene_separate_29180230

genes_list$nervous_gene_separate_29180230 <- genes_list$nervous_gene_weak_29180230[!genes_list$nervous_gene_weak_29180230 %in% genes_list$master_gr_weak] 




