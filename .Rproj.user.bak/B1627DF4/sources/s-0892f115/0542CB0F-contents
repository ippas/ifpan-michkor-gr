install.packages("rentrez")
library(rentrez)
install.packages("XML")
library(XML)

gr_gene_database_preproccesing %>%
  filter(!grepl("omicspred_metabolon", source)) %>%
  filter(!grepl("omicspred_nithingale_", source)) %>%
  mutate(label = ifelse(
    source == "marpiech_tissues_dex",
    paste0(source, "_", gene_list_number),
    paste0(source, "_", tissue, "_", cell)
  )) %>%
  mutate(label = paste0(source, "_", tissue, "_", cell)) %>%
  mutate(label = ifelse(
    source == "marpiech_tissues_dex",
    paste0(source, "_", gene_list_number),
    label
  )) %>% 
  group_by(label) %>%
  mutate(gene_count = n()) %>%
  ungroup() %>% 
  filter(gene_count > 5) %>% 
  mutate(label2 = paste0(label, "_", gene_count)) %>% 
  select(-gene_count) -> papers_data_preprocessing


split(papers_data_preprocessing$hgnc_symbol, papers_data_preprocessing$label) -> papers_gene_list

# calculate chi2 tests
chi2_results_papers <- perform_chi2_tests(c(papers_gene_list), hgnc_symbols_vector_v110)

# prepare mapping vectors to visualisation
# Manually create the mapping for rows
papers_mapping <- c(
  "michkor-cells_NA_astrocyte" = "michkor, astrocyte",
  # "michkor-cells_NA_neuron" = "michkor, neuron",
  "pmid:19801529_lung_A549" = "19801529, lung: A549",
  "pmid:24147833_brain_astrocytes" = "24147833, astrocytes",
  "pmid:24147833_brain_mglia" = "24147833, microglia",
  "pmid:24147833_brain_oligodendrocytes" = "24147833, oligodendrocytes",
  # "pmid:24147833_brain_OPC" = "24147833, OPC",
  "pmid:26606517_embryos_hypothalamic-region_NPSCs" = "26606517, embryos: hypothalamic region, NPSCs",
  "pmid:24926665_lung_ASM" = "24926665, lung: ASM",
  "pmid:24777604_embryos_hypothalamic-region_NPSCs" = "24777604, embryos: hypothalamic region, NPSCs",
  "pmid:28500512_prefrontal-cortex_NA" = "28500512, prefrontal cortex",
  "pmid:31176677_ARC_glia" = "31176677, ARC: glia",
  "pmid:33981007_placenta_CTB" = "33981007, placenta: CTB",
  "pmid:34362910_hippocampus_NA" = "34362910, hippocampus",
  "pmid:35904237_adrenal-gland_Y-1" = "35904237, adrenal-gland: Y1",
  "pmid:35948369_embryos_epi" = "35948369, embryos: epi",
  "pmid:35948369_embryos_mural" = "35948369, embryos: mural",
  "pmid:35948369_embryos_pe" = "35948369, embryos: pe",
  "pmid:35948369_embryos_polar" = "35948369, embryos: polar",
  "master_gr_weak" = "master glucocorticoids genes"
)

clusters_mapping <- c(
  "marpiech_tissues_dex_1" = "cluster 1\nLIV",  
  "marpiech_tissues_dex_10" = "cluster 10\nPIT", 
  "marpiech_tissues_dex_11" = "cluster 11\nnon-specific", 
  "marpiech_tissues_dex_12" = "cluster 12\nKID", 
  "marpiech_tissues_dex_13" = "cluster 13\nHTH,LUN", 
  "marpiech_tissues_dex_14" = "cluster 14\nMUS",
  "marpiech_tissues_dex_15" = "cluster 15\nFAT", 
  "marpiech_tissues_dex_16" = "cluster 16\nFAT", 
  "marpiech_tissues_dex_17" = "cluster DOWN", 
  "marpiech_tissues_dex_18" = "cluster UP", 
  "marpiech_tissues_dex_2" = "cluster 2\nLUN",  
  "marpiech_tissues_dex_3" = "cluster 3\nFAT", 
  "marpiech_tissues_dex_4" = "cluster 4\nnon-specific",  
  "marpiech_tissues_dex_5" = "cluster 5\nKID", 
  "marpiech_tissues_dex_6" = "cluster 6\nPIT",  
  "marpiech_tissues_dex_7" = "cluster 7\nSPL",  
  "marpiech_tissues_dex_8" = "cluster 8\nLIV",  
  "marpiech_tissues_dex_9" = "cluster 9\nKID,PIT"
)


clusters_mapping <- c(
  "marpiech_tissues_dex_1" = "cluster A",  
  "marpiech_tissues_dex_10" = "cluster J", 
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


tissues_mapping <- c(
  "pmid:NA_adrenal-cortex_NA" = "adrenal cortex",
  "pmid:NA_anterior-thigh_NA" = "anterior thigh",
  "pmid:NA_hypothalamus_NA" = "hypothalamus",
  "pmid:NA_kidneys_NA" = "kidneys",
  "pmid:NA_liver_NA" = "liver",
  "pmid:NA_lung_NA" = "lung",
  "pmid:NA_perigonadal-adipose-tissue_NA" = "perigonadal\nadipose tissue",
  "pmid:NA_pituitary-gland_NA" = "pituitary gland",
  "pmid:NA_spleen_NA" = "spleen"
)

# Original list of row names to exclude
tissues_clusters <- c("pmid:NA_adrenal-cortex_NA", "pmid:NA_anterior-thigh_NA", "pmid:NA_hypothalamus_NA", 
                   "pmid:NA_kidneys_NA", "pmid:NA_liver_NA", "pmid:NA_lung_NA", 
                   "pmid:NA_perigonadal-adipose-tissue_NA", "pmid:NA_pituitary-gland_NA", "pmid:NA_spleen_NA",
                   "marpiech_tissues_dex_1", "marpiech_tissues_dex_10", "marpiech_tissues_dex_11",
                   "marpiech_tissues_dex_12", "marpiech_tissues_dex_13", "marpiech_tissues_dex_14",
                   "marpiech_tissues_dex_15", "marpiech_tissues_dex_16", "marpiech_tissues_dex_17",
                   "marpiech_tissues_dex_18", "marpiech_tissues_dex_2", "marpiech_tissues_dex_3",
                   "marpiech_tissues_dex_4", "marpiech_tissues_dex_5", "marpiech_tissues_dex_6",
                   "marpiech_tissues_dex_7", "marpiech_tissues_dex_8", "marpiech_tissues_dex_9")

##################################
# results for clusters vs papers #
##################################

###
processing_overlap_results(data = chi2_results_papers ,
                           rows_to_filter = !rownames(chi2_results_papers $p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[10:27],
                          overlap_threshold = 3,
                          fdr_threshold = 0.01,
                           genes_list = c(papers_gene_list)) -> clusters_papers_data


# tmp$gene_list_sizes <- c(papers_gene_list, genes_list["master_gr_weak"]) %>% sapply(., length)

draw_custom_heatmap(
  clusters_papers_data,
  data_type = "significant_uniq_data",
  col_mapping_vector =  clusters_mapping,
  row_mapping_vector =  setNames(papers_data_preprocessing$label2, papers_data_preprocessing$label),
  fdr_threshold = 0.01,
  fdr_thresholds = c(0.01, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "gray",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = T
)


processing_overlap_results(data = chi2_results_papers ,
                           rows_to_filter = !rownames(chi2_results_papers $p_value_matrix) %in% tissues_clusters,
                           cols_to_filter = tissues_clusters[1:9],
                           genes_list = c(papers_gene_list, genes_list["master_gr_weak"])) -> tissue_papers_data

draw_custom_heatmap(
  tissue_papers_data,
  data_type = "original_data",
  col_mapping_vector = tissues_mapping,
  # row_mapping_vector = phenotypes_mapping_vector,
  fdr_threshold = 0.1,
  fdr_thresholds = c(0.05, 0.0001),
  color_rects =  c("green", "#FF00FF"),
  color_rect = "green",
  lwd_rect = 3,
  alpha_rect = 1,
  apply_filling = F,
  color_filling = "green",
  alpha_filling = 0.6,
  size_filling = 1,
  pch_filling = 16,
  col_significant = F
)



  


  

