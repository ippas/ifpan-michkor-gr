# Wstępne filtrowanie danych
data_filtered <- rbind(
  pgc_mdd2025_no23andMe_eur_processing,
  intesection_results_all_gr_signatures,
  slezak_gene_signatures_preprocessing,
  enrichr_depression_genes_preprocessing
) %>% 
  filter(range_plus_minus %in% c("100000", "all"),
         pvalue < 1e-4)

# Stałe
background_list <- c("GWAS", "Mental_Depression_546/575_DisGeNET")
signature_names <- c(background_list,   
                     "michkor_indicate_T.S1_Cx43DGE_n67"
                     # "T.S4_Astrocytes",
                     # "T.S4_Astrocytes_n140",
                     # # "T.S5_Endothelial_n34",
                     # "T.S6_ExcitatoryNeurons_n66",
                     # "T.S7_InhibitoryNeurons_n35",
                     # "T.S8_Microglia_n28",
                     # "T.S9_Oligodendrocytes_n78",
                     # "T.S10_OPCs_n48",
                     # "T.S11_Pericytes_n35",
                     # "T.S12_VLMCs_n30",
                     # "T.S13_Unclassified_n356"
                     )

plot_signature_density(data_filtered, background_list,
                       background_colors = c("#7F7F7F", "#1B9E77"),
                       signature_names, log_y_axis = FALSE)



data_filtered %>% 
  filter(signature_name == "michkor_indicate_T.S1_Cx43DGE_n67") %>% 
  select(rsID, pvalue, gene_symbol) %>% 
  group_by(gene_symbol) %>% 
  mutate(n_rsID = n()) %>% 
  slice_min(pvalue, with_ties = FALSE)
