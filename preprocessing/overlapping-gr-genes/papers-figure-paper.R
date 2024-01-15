# gr_gene_database_preproccesing %>% 
#   filter(source == "pmid:33203447") %>% 
#   mutate(tmp = fdr) %>% 
#   mutate(fdr = regulation) %>% 
#   mutate(regulation = tmp) %>% 
#   select(c("index", "gene_name", "gene_list_index", "gene_list_number", "source", 
#            "ensembl_gene_id", "ensembl_transcript_id", "refseq_mrna_id", "hgnc_symbol", 
#            "alias", "info", "tissue", "cell", "species", "environment", "treatment", 
#            "dose", "time", "log2ratio", "fdr", "statistical_method", "treatment_type", 
#            "regulation")) -> tmp_chondrocytes

# gr_gene_database_preproccesing %>% 
#   filter(source != "pmid:33203447") %>% 
#   rbind(., tmp_chondrocytes) -> gr_gene_database_preproccesing 

# gr_gene_database_preproccesing %>% 
#   filter(!grepl("omicspred_metabolon", source)) %>%
#   filter(!grepl("omicspred_nithingale_", source)) %>%
#   .$treatment %>% unique()
#   mutate(label = paste(source, tissue, cell, regulation, treatment, sep = "_")) %>%
#   .$label %>% unique

# tmp code to repair data for PMID:24777604
gr_gene_database_preproccesing %>% 
  filter(source == "pmid:24777604") %>%
  mutate(comparison = ifelse(comparison == "Cav1-KO-DEX_vs_C57-Ethanol", "Cav1-KO-Dex_vs_Cav1-KO-Ethanol",
                             ifelse(comparison == "C57-DEX_vs_C57-Ethanol", "C57-Dex_vs_C57-Ethanol", "Cav1-KO-Ethanol_vs_C57BL6-Dex"))) %>% 
  mutate(cell = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "NPSCs", 
                           ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "NPSCs-KO-Cav1", "NPSCs|NPSCs-KO-Cav1")),
         treatment = ifelse(comparison == "C57-Dex_vs_C57-Ethanol", "dexamethasone", 
                                ifelse(comparison == "Cav1-KO-Dex_vs_Cav1-KO-Ethanol", "dexamethasone", "vehicle-ethanol"))) -> tmp_npscs_kocav1
  

gr_gene_database_preproccesing %>%
  filter(source != "pmid:24777604") %>%
  rbind(., tmp_npscs_kocav1) -> gr_gene_database_preproccesing

gr_gene_database_preproccesing %>%
  filter(!grepl("omicspred_metabolon", source)) %>%
  filter(!grepl("omicspred_nithingale_", source)) %>%
  mutate(label = ifelse(
    source == "marpiech_tissues_dex",
    paste0(source, "_", gene_list_number),
    paste0(source, "_", tissue, "_", cell)
  )) %>%
  filter(!(treatment %in% c("vitamin-d3", "vehicle-DMSO"))) %>% 
  # mutate(label = paste0(source, "_", tissue, "_", cell)) %>%
  mutate(label = paste(source, tissue,cell, regulation,sep = "_")) %>% 
  mutate(label = ifelse(
    source == "marpiech_tissues_dex",
    paste0(source, "_", gene_list_number),
    label
  )) %>% 
  group_by(label) %>%
  mutate(gene_count = n()) %>%
  ungroup() %>% 
  filter(gene_count > 5) %>% 
  # mutate(label2 = paste0(label, "_", gene_count)) %>% 
  select(-gene_count) -> papers_data_preprocessing

papers_data_preprocessing %>% 
  mutate(fdr = as.numeric(fdr),
         log2ratio = as.numeric(log2ratio)) %>% 
  mutate(abs_log2ratio = abs(log2ratio)) %>% 
  # filter(tissue != "kidney" | (tissue == "kidney" & fdr < 0.001)) %>%
  filter(source == "pmid:23593176", fdr < 0.0001, abs_log2ratio > 1) %>% 
  select(-abs_log2ratio) -> tmp_kidney

papers_data_preprocessing %>% 
  mutate(fdr = as.numeric(fdr),
         log2ratio = as.numeric(log2ratio)) %>% 
  filter(source != "pmid:23593176") %>% 
  rbind(., tmp_kidney) -> papers_data_preprocessing
 



split(papers_data_preprocessing$hgnc_symbol, papers_data_preprocessing$label) -> papers_gene_list

# calculate chi2 tests
lapply(papers_gene_list, unique) -> papers_gene_list

chi2_results_papers <- perform_chi2_tests(c(papers_gene_list), hgnc_symbols_vector_v110)


papers_data_preprocessing %>% select(label) %>% filter(grepl("pmid:NA", label)) %>% unique() %>% as.vector() %>% .$label -> remove_list

processing_overlap_results(data = chi2_results_papers ,
                           rows_to_filter = !rownames(chi2_results_papers$p_value_matrix) %in% c(remove_list, tissues_clusters),
                           # rows_to_filter = remove_list,
                           cols_to_filter = tissues_clusters[10:27],
                           overlap_threshold = 3,
                           fdr_threshold = 0.01,
                           genes_list = c(papers_gene_list)) -> clusters_papers_data


clusters_mapping <- c(
  "marpiech_tissues_dex_1" = "cluster A",  
  "marpiech_tissues_dex_10" = "cluster A", 
  "marpiech_tissues_dex_11" = "cluster K", 
  "marpiech_tissues_dex_12" = "cluster L", 
  "marpiech_tissues_dex_13" = "cluster M", 
  "marpiech_tissues_dex_14" = "cluster N",
  "marpiech_tissues_dex_15" = "cluster O", 
  "marpiech_tissues_dex_16" = "cluster P", 
  "marpiech_tissues_dex_17" = "DOWN", 
  "marpiech_tissues_dex_18" = "UP", 
  "marpiech_tissues_dex_2" = "cluster B",  
  "marpiech_tissues_dex_3" = "cluster C", 
  "marpiech_tissues_dex_4" = "cluster D",  
  "marpiech_tissues_dex_5" = "cluster E", 
  "marpiech_tissues_dex_6" = "cluster F",  
  "marpiech_tissues_dex_7" = "cluster G",  
  "marpiech_tissues_dex_8" = "cluster H",  
  "marpiech_tissues_dex_9" = "cluster I"
)

clusters_mapping %>%
  setNames(names(.)) %>%
  purrr::map_chr(~ paste("      ", ., sep = "")) -> clusters_mapping


papers_mapping <- c("michkor-cells_NA_astrocyte_down" = "astrocytes, down, PMID:28381250", 
                    "michkor-cells_NA_astrocyte_up" = "astrocytes, up, PMID:28381250", 
                    "pmid:18489715_hypothalamus_NA_up" = "hypothalamus, up, PMID:18489715",
                    "pmid:19801529_lung_A549_down" = "lung, A549, down, PMID:19801529", 
                    "pmid:19801529_lung_A549_up" = "lung, A549, up, PMID:19801529", 
                    "pmid:21187916_adipocyte-tissue_3T3-L1_down" = "adipocyte tissue, 3T3-L1, down, PMID:21187916",
                    "pmid:21187916_adipocyte-tissue_3T3-L1_up" = "adipocyte tissue, 3T3-L1, up, PMID:21187916", 
                    "pmid:26606517_embryos_hypothalamic-region_NPSCs_up" = "embryos, hypothalamic region, up, PMID:26606517", 
                    "pmid:21846803_hippocampus_dentate-gyrus_up" = "hippocampus, dentate gyrus, up, PMID:21846803",
                    "pmid:22108209_adipose_NA_up" = "adipose, up, PMID:22108209", 
                    "pmid:22108209_skeletal-muscle_NA_up" = "skeletal muscle, up, PMID:22108209", 
                    "pmid:23110767_cortex_astrocytes_down" = "cortex, astrocytes, down, PMID:23110767",
                    "pmid:23110767_cortex_astrocytes_up" = "cortex, astrocytes, up, PMID:23110767", 
                    "pmid:23339081_striatum_astrocytes_up" = "striatum, astrocytes, up, PMID:23339081", 
                    "pmid:23339081_striatum_neurones_down" = "striatum, neurones, down, PMID:23339081",
                    "pmid:23593176_kidney_HPCs-differentiated_up" = "kidney, HPCs differentiated, up, PMID:23593176", 
                    "pmid:23633533_hippocampus_dentate-gyrus_up" = "hippocampus, dentate gyrus, up, PMID:23633533", 
                    "pmid:23736296_hippocampus_ca1_down" = "hippocampus, CA1, down, PMID:23736296",
                    "pmid:23736296_hippocampus_ca1_up" = "hippocampus, CA1, up, PMID:23736296", 
                    "pmid:24147833_brain_astrocytes_up" = "brain, astrocytes, up, PMID:24147833", 
                    "pmid:24342991_hippocampus_NA_down" = "hippocampus, down, PMID:24342991",
                    "pmid:24342991_hippocampus_NA_up" = "hippocampus, up, PMID:24342991", 
                    "pmid:24395918_blood_macrophages_up" = "blood, macrophages, up, PMID:24395918", 
                    "pmid:24777604_embryos_hypothalamic-region_NPSCs_up" = "embryos, hypothalamic region, NPSCs, up, PMID:24777604",
                    "pmid:24926665_lung_ASM_down" = "lung, ASM, down, PMID:24926665", 
                    "pmid:28432191_skeletal-muscle_myotubes_down" = "skeletal muscle, myotubes, down, PMID:28432191", 
                    "pmid:28432191_skeletal-muscle_myotubes_up" = "skeletal muscle, myotubes, up, PMID:28432191",
                    "pmid:28500512_prefrontal-cortex_NA_down" = "prefrontal cortex, down, PMID:28500512", 
                    "pmid:28500512_prefrontal-cortex_NA_up" = "prefrontal cortex, up, PMID:28500512", 
                    "pmid:33203447_knee-cartilage_osteoarthritis-chon_down" = "knee cartilage, osteoarthritis chondrocytes, down, PMID:33203447",
                    "pmid:33203447_knee-cartilage_osteoarthritis-chon_up" = "knee cartilage, osteoarthritis chondrocytes, up, PMID:33203447", 
                    "pmid:33923915_lung_A549_up" = "lung, A549, up, PMID:33923915", 
                    "pmid:34362910_hippocampus_NA_down" = "hippocampus, down, PMID:34362910",
                    "pmid:34362910_hippocampus_NA_up" = "hippocampus, up, PMID:34362910", 
                    "pmid:34731602_liver_NA_1" = "liver,  cluster 1, PMID:34731602", 
                    "pmid:34731602_liver_NA_2" = "liver,  cluster 2, PMID:34731602",
                    "pmid:34731602_liver_NA_3" = "liver,  cluster 3, PMID:34731602", 
                    "pmid:34731602_liver_NA_4" = "liver,  cluster 4, PMID:34731602", 
                    "pmid:35904237_adrenal-gland_Y-1_down" = "adrenal gland, Y-1, down, PMID:35904237",
                    "pmid:35904237_adrenal-gland_Y-1_up" = "adrenal gland, Y-1, up, PMID:35904237", 
                    "pmid:36284713_blood_THP-1_down" = "blood, THP-1, down, PMID:36284713", 
                    "pmid:36284713_blood_THP-1_up" = "blood, THP-1, up, PMID:36284713",
                    "pmid:34735568_kidney_NA_up" = "kidney, up, PMID:34735568",
                    "pmid:34735568_kidney_NA_down" = "kidney, down, PMID:34735568")

# Replacing " ID" with "                               ID"
papers_mapping <- gsub(" ID", "                               ID", papers_mapping)

# Printing the modified vector
print(papers_mapping)


papers_mapping <- c(
  "michkor-cells_NA_astrocyte_down" = "astrocytes,   \u2212, PMID:28381250", 
  "michkor-cells_NA_astrocyte_up" = "astrocytes, \u002B, PMID:28381250", 
  "pmid:18489715_hypothalamus_NA_up" = "hypothalamus, \u002B, PMID:18489715",
  "pmid:19801529_lung_A549_down" = "lung, A549,   \u2212, PMID:19801529", 
  "pmid:19801529_lung_A549_up" = "lung, A549, \u002B, PMID:19801529", 
  "pmid:21187916_adipocyte-tissue_3T3-L1_down" = "adipocyte tissue, 3T3-L1,  \u2212, PMID:21187916",
  "pmid:21187916_adipocyte-tissue_3T3-L1_up" = "adipocyte tissue, 3T3-L1, \u002B, PMID:21187916", 
  "pmid:26606517_embryos_hypothalamic-region_NPSCs_up" = "embryos, hypothalamic region, \u002B, PMID:26606517", 
  "pmid:21846803_hippocampus_dentate-gyrus_up" = "hippocampus, dentate gyrus, \u002B, PMID:21846803",
  "pmid:22108209_adipose_NA_up" = "adipose, \u002B, PMID:22108209", 
  "pmid:22108209_skeletal-muscle_NA_up" = "skeletal muscle, \u002B, PMID:22108209", 
  "pmid:23110767_cortex_astrocytes_down" = "cortex, astrocytes,  \u2212, PMID:23110767",
  "pmid:23110767_cortex_astrocytes_up" = "cortex, astrocytes, \u002B, PMID:23110767", 
  "pmid:23339081_striatum_astrocytes_up" = "striatum, astrocytes, \u002B, PMID:23339081", 
  "pmid:23339081_striatum_neurones_down" = "striatum, neurones,  \u2212, PMID:23339081",
  "pmid:23593176_kidney_HPCs-differentiated_up" = "kidney, HPCs differentiated, \u002B, PMID:23593176", 
  "pmid:23633533_hippocampus_dentate-gyrus_up" = "hippocampus, dentate gyrus, \u002B, PMID:23633533", 
  "pmid:23736296_hippocampus_ca1_down" = "hippocampus, CA1,  \u2212, PMID:23736296",
  "pmid:23736296_hippocampus_ca1_up" = "hippocampus, CA1, \u002B, PMID:23736296", 
  "pmid:24147833_brain_astrocytes_up" = "brain, astrocytes, \u002B, PMID:24147833", 
  "pmid:24342991_hippocampus_NA_down" = "hippocampus,  \u2212, PMID:24342991",
  "pmid:24342991_hippocampus_NA_up" = "hippocampus, \u002B, PMID:24342991", 
  "pmid:24395918_blood_macrophages_up" = "blood, macrophages, \u002B, PMID:24395918", 
  "pmid:24777604_embryos_hypothalamic-region_NPSCs_up" = "embryos, hypothalamic region, NPSCs, \u002B, PMID:24777604",
  "pmid:24926665_lung_ASM_down" = "lung, ASM,  \u2212, PMID:24926665", 
  "pmid:28432191_skeletal-muscle_myotubes_down" = "skeletal muscle, myotubes,  \u2212, PMID:28432191", 
  "pmid:28432191_skeletal-muscle_myotubes_up" = "skeletal muscle, myotubes, \u002B, PMID:28432191",
  "pmid:28500512_prefrontal-cortex_NA_down" = "prefrontal cortex,  \u2212, PMID:28500512", 
  "pmid:28500512_prefrontal-cortex_NA_up" = "prefrontal cortex, \u002B, PMID:28500512", 
  "pmid:33203447_knee-cartilage_osteoarthritis-chon_down" = "knee cartilage, osteoarthritis chondrocytes,  \u2212, PMID:33203447",
  "pmid:33203447_knee-cartilage_osteoarthritis-chon_up" = "knee cartilage, osteoarthritis chondrocytes, \u002B, PMID:33203447", 
  "pmid:33923915_lung_A549_up" = "lung, A549, \u002B, PMID:33923915", 
  "pmid:34362910_hippocampus_NA_down" = "hippocampus,  \u2212, PMID:34362910",
  "pmid:34362910_hippocampus_NA_up" = "hippocampus, \u002B, PMID:34362910", 
  "pmid:34731602_liver_NA_1" = "liver GR knockout, cluster 1, PMID:34731602", 
  "pmid:34731602_liver_NA_2" = "liver GR knockout, cluster 2, PMID:34731602",
  "pmid:34731602_liver_NA_3" = "liver GR knockout, cluster 3, PMID:34731602", 
  "pmid:34731602_liver_NA_4" = "liver GR knockout, cluster 4, PMID:34731602", 
  "pmid:35904237_adrenal-gland_Y-1_down" = "adrenal gland, Y-1,  \u2212, PMID:35904237",
  "pmid:35904237_adrenal-gland_Y-1_up" = "adrenal gland, Y-1, \u002B, PMID:35904237", 
  "pmid:36284713_blood_THP-1_down" = "blood, THP-1,  \u2212, PMID:36284713", 
  "pmid:36284713_blood_THP-1_up" = "blood, THP-1, \u002B, PMID:36284713",
  "pmid:34735568_kidney_NA_up" = "kidney, \u002B, PMID:34735568",
  "pmid:34735568_kidney_NA_down" = "kidney,  \u2212, PMID:34735568"
)

papers_mapping %>%
  setNames(names(.)) %>%
  sapply(function(x) str_replace(x, ",(?=[^,]*$)", "    ")) %>%
  setNames(names(papers_mapping)) %>%   sapply(function(x) str_replace(x, ",(?=[^,]*$)", "    ")) %>%
  setNames(names(papers_mapping)) -> papers_mapping 
  
  # purrr::map_chr(~ paste("    ", ., sep = "")) -> papers_mapping

# papers_mapping <- paste("   ", papers_mapping, sep = "")


# papers_mapping <- c(
#   "michkor-cells_NA_astrocyte_down" = "astrocytes, PMID:28381250",
#   "michkor-cells_NA_astrocyte_up" = "astrocytes, PMID:28381250",
#   "pmid:18489715_hypothalamus_NA_up" = "hypothalamus, PMID:18489715",
#   "pmid:19801529_lung_A549_down" = "lung, A549, PMID:19801529",
#   "pmid:19801529_lung_A549_up" = "lung, A549, PMID:19801529",
#   "pmid:21187916_adipocyte-tissue_3T3-L1_down" = "adipocyte tissue, 3T3-L1, PMID:21187916",
#   "pmid:21187916_adipocyte-tissue_3T3-L1_up" = "adipocyte tissue, 3T3-L1, PMID:21187916",
#   "pmid:26606517_embryos_hypothalamic-region_NPSCs_up" = "embryos, hypothalamic region, PMID:26606517",
#   "pmid:21846803_hippocampus_dentate-gyrus_up" = "hippocampus, dentate gyrus, PMID:21846803",
#   "pmid:22108209_adipose_NA_up" = "adipose, PMID:22108209",
#   "pmid:22108209_skeletal-muscle_NA_up" = "skeletal muscle, PMID:22108209",
#   "pmid:23110767_cortex_astrocytes_down" = "cortex, astrocytes, PMID:23110767",
#   "pmid:23110767_cortex_astrocytes_up" = "cortex, astrocytes, PMID:23110767",
#   "pmid:23339081_striatum_astrocytes_up" = "striatum, astrocytes, PMID:23339081",
#   "pmid:23339081_striatum_neurones_down" = "striatum, neurones, PMID:23339081",
#   "pmid:23593176_kidney_HPCs-differentiated_up" = "kidney, HPCs differentiated, PMID:23593176",
#   "pmid:23633533_hippocampus_dentate-gyrus_up" = "hippocampus, dentate gyrus, PMID:23633533",
#   "pmid:23736296_hippocampus_ca1_down" = "hippocampus, CA1, PMID:23736296",
#   "pmid:23736296_hippocampus_ca1_up" = "hippocampus, CA1, PMID:23736296",
#   "pmid:24147833_brain_astrocytes_up" = "brain, astrocytes, PMID:24147833",
#   "pmid:24342991_hippocampus_NA_down" = "hippocampus, PMID:24342991",
#   "pmid:24342991_hippocampus_NA_up" = "hippocampus, PMID:24342991",
#   "pmid:24395918_blood_macrophages_up" = "blood, macrophages, PMID:24395918",
#   "pmid:24777604_embryos_hypothalamic-region_NPSCs_up" = "embryos, hypothalamic region, NPSCs, PMID:24777604",
#   "pmid:24926665_lung_ASM_down" = "lung, ASM, PMID:24926665",
#   "pmid:28432191_skeletal-muscle_myotubes_down" = "skeletal muscle, myotubes, PMID:28432191",
#   "pmid:28432191_skeletal-muscle_myotubes_up" = "skeletal muscle, myotubes, PMID:28432191",
#   "pmid:28500512_prefrontal-cortex_NA_down" = "prefrontal cortex, PMID:28500512",
#   "pmid:28500512_prefrontal-cortex_NA_up" = "prefrontal cortex, PMID:28500512",
#   "pmid:33203447_knee-cartilage_osteoarthritis-chon_down" = "knee cartilage, osteoarthritis chondrocytes, PMID:33203447",
#   "pmid:33203447_knee-cartilage_osteoarthritis-chon_up" = "knee cartilage, osteoarthritis chondrocytes, PMID:33203447",
#   "pmid:33923915_lung_A549_up" = "lung, A549, PMID:33923915",
#   "pmid:34362910_hippocampus_NA_down" = "hippocampus, PMID:34362910",
#   "pmid:34362910_hippocampus_NA_up" = "hippocampus, PMID:34362910",
#   "pmid:34731602_liver_NA_1" = "liver, cluster 1, PMID:34731602",
#   "pmid:34731602_liver_NA_2" = "liver, cluster 2, PMID:34731602",
#   "pmid:34731602_liver_NA_3" = "liver, cluster 3, PMID:34731602",
#   "pmid:34731602_liver_NA_4" = "liver, cluster 4, PMID:34731602",
#   "pmid:35904237_adrenal-gland_Y-1_down" = "adrenal gland, Y-1, PMID:35904237",
#   "pmid:35904237_adrenal-gland_Y-1_up" = "adrenal gland, Y-1, PMID:35904237",
#   "pmid:36284713_blood_THP-1_down" = "blood, THP-1, PMID:36284713",
#   "pmid:36284713_blood_THP-1_up" = "blood, THP-1, PMID:36284713",
#   "pmid:34735568_kidney_NA_up" = "kidney, PMID:34735568",
#   "pmid:34735568_kidney_NA_down" = "kidney, PMID:34735568"
# )
# 
# sapply(names(papers_mapping), function(name) {
#   if (grepl("_down", name)) {
#     return("darkblue")
#   } else if (grepl("_up", name)) {
#     return("darkred")
#   } else {
#     return("black")
#   }
# }) -> row_colors
# 
# sapply(names(clusters_mapping ), function(name) {
#   if (grepl("_17", name)) {
#     return("darkblue")
#   } else if (grepl("_18", name)) {
#     return("darkred")
#   } else {
#     return("black")
#   }
# }) -> col_colors
# 
# col_colors <- c("marpiech_tissues_dex_12" = "darkred", "marpiech_tissues_dex_15" = "darkred", 
#                 "marpiech_tissues_dex_17" = "darkblue", "marpiech_tissues_dex_18" = "darkred", 
#                 "marpiech_tissues_dex_4" = "darkblue", "marpiech_tissues_dex_6" = "darkblue", 
#                 "marpiech_tissues_dex_11" = "darkred", "marpiech_tissues_dex_13" = "darkred", 
#                 "marpiech_tissues_dex_16" = "darkred", "marpiech_tissues_dex_9" = "darkred", 
#                 "marpiech_tissues_dex_3" = "darkblue", "marpiech_tissues_dex_10" = "darkred", 
#                 "marpiech_tissues_dex_7" = "darkred", "marpiech_tissues_dex_8" = "darkred", 
#                 "marpiech_tissues_dex_2" = "darkblue")

manual_filter_overlap_results(data = clusters_papers_data, 
                              rows_to_remove = c("pmid:24777604_embryos_hypothalamic-region_NPSCs-KO-Cav1_up",
                                                 "pmid:26606517_embryos_hypothalamic-region_NPSCs_up")) -> clusters_papers_data

clusters_papers_data$significant_uniq_data

  
draw_custom_heatmap(
  clusters_papers_data,
  data_type = "manual_filter_overlap_results",
  # data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  papers_mapping,
  # column_names_side = "top",
  fdr_threshold = 0.01,
  fdr_thresholds = c(0.01, 0.0001),
  # color_rects =  c("green", "#FF00FF"),
  color_rects =  c("#4C8D05", "#66023C"),
  color_rect = "green",
  lwd_rect = 2,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T,
  row_dend_width = unit(4, "cm"),  # Adjust row dendrogram width
  column_dend_height = unit(3, "cm"),  # Adjust column dendrogram height
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16),
  column_names_rot = 45
)



getwd()
dev.off()

svg("results/figures/cluster-papers-tissue-cell2.svg", width = 19, height = 21)

draw_custom_heatmap(
  clusters_papers_data,
  data_type = "manual_filter_overlap_results",
  # data_type = "significant_uniq_data",
  palette = c(
    "pastel_blue"       = "white",
    "pastel_light_blue" = "#f8dedd",
    "white"        = "#f1bcbb",
    "pastel_orange"= "#edacab",
    "pastel_red"   = "#e68a89"
  ),
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  papers_mapping,
  fdr_threshold = 0.01,
  fdr_thresholds = c(0.01, 0.0001),
  # color_rects =  c("green", "#FF00FF"),
  color_rects =  c("#4C8D05", "#66023C"),
  color_rect = "green",
  lwd_rect = 2,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T,
  row_dend_width = unit(4, "cm"),  # Adjust row dendrogram width
  column_dend_height = unit(3, "cm"),  # Adjust column dendrogram height
  row_names_gp = gpar(fontsize = 16),
  column_names_gp = gpar(fontsize = 16),
  column_names_rot = 45,
  column_dend_side = "bottom",
  column_names_side = "top"
)

dev.off()




# Extracting keys and values
keys <- names(clusters_mapping)
values <- clusters_mapping

# Splitting the values into two columns: 'cluster' and 'tissue'
split_values <- strsplit(values, "\n")
clusters <- sapply(split_values, `[`, 1)
tissues <- sapply(split_values, `[`, 2)

# Handling missing tissue values
tissues[is.na(tissues)] <- ""

# Creating a dataframe
df <- data.frame(marpiech_cluster = keys, cluster = clusters, tissue = tissues, stringsAsFactors = FALSE)

# View the dataframe
print(df)

write.csv(df, "results/tables/cluster-description.tsv", col.names = T, row.names = F, sep = "\t")



clusters_papers_data2$significant_uniq_data$df %>% .$overlap_genes


