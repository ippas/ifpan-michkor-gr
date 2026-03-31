source("preprocessing/functions/R/function-gr-list-processing.R")

################################################################################
# prepare categorized_gene_lists
################################################################################
categorized_gene_lists <- list()


# gene lists clustered
categorized_gene_lists$marpiech_cluster_dex <- marpiech_data_preprocessing %>% 
  filter((source %in% c("marpiech_clusters_dex"))) %>% 
  mutate(label = paste0("cluster_", gene_list_number)) %>% 
  select(c(hgnc_symbol, label)) %>% 
  split(.$label, .$hgnc_symbol) %>% 
  lapply(., unique) %>% 
  lapply(., function(x){x$hgnc_symbol})

################################################################################
# gene lists phenotypes
################################################################################
list.files("data/prs-models-pan-biobank-uk/", full.names = TRUE) %>%
  # Set the names of the list elements to the file names without the extension
  set_names(map(., ~ basename(.x) %>% tools::file_path_sans_ext())) %>%
  # Read YAML content from each file
  map(read_yaml) %>%
  # Extract genes information and name the list elements after the files
  map(~ .x$description$genes) %>% lapply(., unique) -> phenotypes_PanUkBiobank

phenotypes_PanUkBiobank %>%
  lapply(length) %>%
  enframe(name = "label", value = "gene_count") %>% 
  unnest(gene_count) %>% 
  separate(col = label, into = c("label", "GWAS_pvalue"), sep = "(?<=-EUR)-", extra = "merge", fill = "right") -> metadata

names(phenotypes_PanUkBiobank) <- gsub("-1e-05|-1e-06|-1e-07|-1e-08", "", names(phenotypes_PanUkBiobank))

categorized_gene_lists$phenotypes_PanUkBiobank$metadata <- metadata
categorized_gene_lists$phenotypes_PanUkBiobank$gene_lists <- phenotypes_PanUkBiobank

rm(metadata, phenotypes_PanUkBiobank)

categorized_gene_lists$phenotypes_PanUkBiobank_1e08 <- categorized_gene_lists$phenotypes_PanUkBiobank %>% filter_gene_list(GWAS_pvalue == "1e-08")


################################################################################
# gene lists for metabolism 
################################################################################

gr_gene_database_preproccesing %>%
  filter(grepl("omicspred_metabolon", source)) %>% 
  extract_keys_values("info", c("Biochemical Name", "Metabolon ID", "Super Pathway", "Sub Pathway")) %>% 
  mutate(OminispreadID = str_remove(gene_list_index, "Metabolon_")) %>% 
  rename(Biochemical_Name = "Biochemical Name") %>% 
  rename(Metabolon_ID = "Metabolon ID") %>% 
  rename(Super_Pathway = "Super Pathway") %>% 
  rename(Sub_Pathway = "Sub Pathway") %>% 
  mutate(OminispreadID = str_remove(OminispreadID, "_model")) -> metabolome_data_preprocessing

gr_gene_database_preproccesing %>%
  filter(grepl("omicspred_nithingale", source)) %>% 
  extract_keys_values("info", c("Biomarker Name", "Trait ID", "Group", "Subgroup")) %>% 
  rename(Biomarker_Name = "Biomarker Name") %>% 
  rename(Trait_ID = "Trait ID") %>% 
  mutate(OminispreadID = str_remove(gene_list_index, "Nithingale_")) %>% 
  mutate(OminispreadID = str_remove(OminispreadID, "_model")) -> nightingale_data_preprocessing

refine_and_label_gene_lists(
  metabolome_data_preprocessing,
  columns = c("Biochemical_Name"),
  keep_column = c("OminispreadID", "Metabolon_ID", "Super_Pathway", "Sub_Pathway"),
  size_threshold = 3
) -> categorized_gene_lists$metabolome_gene_lists

refine_and_label_gene_lists(
  nightingale_data_preprocessing,
  columns = c("Biomarker_Name"),
  keep_column = c("OminispreadID", "Trait_ID", "Group", "Subgroup"),
  size_threshold = 3
) -> categorized_gene_lists$nightingale_gene_lists


################################################################################
# gene list for KEGG pathways
################################################################################
load("results/kegg-pathway-df.RData")


refine_and_label_gene_lists(
  kegg_pathway_df,
  columns = c("name"),
  keep_column = c("kegg_pathway_id", "class"),
  genes_column = "gene_symbol",
  size_threshold = 2
) -> categorized_gene_lists$kegg_pathways_gene_lists




################################################################################
metabolome_gene_counts = categorized_gene_lists$metabolome_gene_lists$metadata$gene_count


list(metabolome_gene_counts = categorized_gene_lists$metabolome_gene_lists$metadata$gene_count,
           nightingale_gene_counts = categorized_gene_lists$nightingale_gene_lists$metadata$gene_count) -> tmp_list

lapply(tmp_list, function(x){summary(x)})


lapply(tmp_list, function(x){hist(x)}) %>% facet_grid()

categorized_gene_lists$metabolome_gene_lists$metadata$SuperPathway %>% table

categorized_gene_lists$nightingale_gene_lists$metadata$Subgroup %>% table


categorized_gene_lists$metabolome_gene_lists$gene_lists %>% unname %>% unlist %>% unique() %>% length()

categorized_gene_lists$nightingale_gene_lists$gene_lists %>% unname %>% unlist %>% unique() %>% length()
