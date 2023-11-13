# Load necessary libraries
library(pheatmap)

# ---------------------
# 1. Perform Chi-squared Tests on Gene List
# ---------------------

# Perform chi-squared test on gr_gene_list and hgnc_symbols_vector
chi2_results <- perform_chi2_tests(c(gr_gene_list, genes_list["master_gr_weak"]), hgnc_symbols_vector)

# ---------------------
# 2. Visualize Gene Expression Profiles
# ---------------------

# Define a subset of rows for heatmap visualization
selected_rows <- c(19, 20, 30:ncol(chi2_results$chi2_value_matrix))
# Assuming you have the original column names as:
expression_data_clustered <- log2(chi2_results$chi2_value_matrix[selected_rows, 1:18] + 1)


# Manually create the mapping for rows
row_mapping <- c(
  "michkor-cells_NA_astrocyte" = "michkor, astrocyte",
  "michkor-cells_NA_neuron" = "michkor, neuron",
  "pmid:19801529_lung_A549" = "19801529, lung: A549",
  "pmid:24147833_brain_astrocytes" = "24147833, astrocytes",
  "pmid:24147833_brain_mglia" = "24147833, microglia",
  "pmid:24147833_brain_oligodendrocytes" = "24147833, oligodendrocytes",
  "pmid:24147833_brain_OPC" = "24147833, OPC",
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

cluster_vector <- c(
  "marpiech_tissues_dex_1" = "cluster 1",  
  "marpiech_tissues_dex_10" = "cluster 10", 
  "marpiech_tissues_dex_11" = "cluster 11", 
  "marpiech_tissues_dex_12" = "cluster 12", 
  "marpiech_tissues_dex_13" = "cluster 13", 
  "marpiech_tissues_dex_14" = "cluster 14",
  "marpiech_tissues_dex_15" = "cluster 15", 
  "marpiech_tissues_dex_16" = "cluster 16", 
  "marpiech_tissues_dex_17" = "cluster 17", 
  "marpiech_tissues_dex_18" = "cluster 18", 
  "marpiech_tissues_dex_2" = "cluster 2",  
  "marpiech_tissues_dex_3" = "cluster 3", 
  "marpiech_tissues_dex_4" = "cluster 4",  
  "marpiech_tissues_dex_5" = "cluster 5", 
  "marpiech_tissues_dex_6" = "cluster 6",  
  "marpiech_tissues_dex_7" = "cluster 7",  
  "marpiech_tissues_dex_8" = "cluster 8",  
  "marpiech_tissues_dex_9" = "cluster 9"
)

tissue_vector <- c(
  "pmid: NA_adrenal-cortex_NA" = "adrenal cortex",
  "pmid: NA_anterior-thigh_NA" = "anterior thigh",
  "pmid: NA_hypothalamus_NA" = "hypothalamus",
  "pmid: NA_kidneys_NA" = "kidneys",
  "pmid: NA_liver_NA" = "liver",
  "pmid: NA_lung_NA" = "lung",
  "pmid: NA_perigonadal-adipose-tissue_NA" = "perigonadal\nadipose tissue",
  "pmid: NA_pituitary-gland_NA" = "pituitary gland",
  "pmid: NA_spleen_NA" = "spleen"
)


# a. Gene expression profile with clustering
selected_rows <- c(19, 20, 30:ncol(chi2_results$chi2_value_matrix))

# Original list of row names to exclude
excluded_rows <- c("pmid:NA_adrenal-cortex_NA", "pmid:NA_anterior-thigh_NA", "pmid:NA_hypothalamus_NA", 
                   "pmid:NA_kidneys_NA", "pmid:NA_liver_NA", "pmid:NA_lung_NA", 
                   "pmid:NA_perigonadal-adipose-tissue_NA", "pmid:NA_pituitary-gland_NA", "pmid:NA_spleen_NA",
                   "marpiech_tissues_dex_1", "marpiech_tissues_dex_10", "marpiech_tissues_dex_11",
                   "marpiech_tissues_dex_12", "marpiech_tissues_dex_13", "marpiech_tissues_dex_14",
                   "marpiech_tissues_dex_15", "marpiech_tissues_dex_16", "marpiech_tissues_dex_17",
                   "marpiech_tissues_dex_18", "marpiech_tissues_dex_2", "marpiech_tissues_dex_3",
                   "marpiech_tissues_dex_4", "marpiech_tissues_dex_5", "marpiech_tissues_dex_6",
                   "marpiech_tissues_dex_7", "marpiech_tissues_dex_8", "marpiech_tissues_dex_9")


# select significant results
extract_data(chi2_results, !rownames(chi2_results$p_value_matrix) %in% excluded_rows, 1:18) %>% 
  rename(paper = Var1, cluster = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05)

# select overlapping genes
extract_data(chi2_results, !rownames(chi2_results$p_value_matrix) %in% excluded_rows, 1:18) %>% 
  rename(paper = Var1, cluster = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

  
# b. Gene expression profile without clustering for specific columns
selected_columns <- c(1:18)
expression_data_specific <- log(chi2_results$chi2_value_matrix[!rownames(chi2_results$p_value_matrix) %in% excluded_rows, selected_columns] + 1)

# heatmap
pheatmap(
  expression_data_specific,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_row = row_mapping[rownames(expression_data_specific)],
  labels_col = cluster_vector[colnames(expression_data_specific)],
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)

# ---------------------
# 3. Tissue Specific Gene Expression Visualization
# ---------------------

# filter significant comparison
extract_data(chi2_results, !rownames(chi2_results$p_value_matrix) %in% excluded_rows, excluded_rows[1:9]) %>% 
  rename(paper = Var1, tissue = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) 

# select overlapping genes
extract_data(chi2_results, !rownames(chi2_results$p_value_matrix) %in% excluded_rows, excluded_rows[1:9]) %>% 
  rename(paper = Var1, tissue = Var2) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))

extract_data(chi2_results, !rownames(chi2_results$p_value_matrix) %in% excluded_rows, excluded_rows[1:9]) 

chi2_results$p_value_matrix[selected_rows, 21:29] %>% 
  melt() %>% 
  set_colnames(c("list_genes", "tissues_marpiech", "p_value")) %>% 
  mutate(tissues_marpiech = str_replace(tissues_marpiech, "pmid: NA_", "")) %>%
  mutate(tissues_marpiech = str_replace(tissues_marpiech, "_NA$", "")) %>% 
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% 
  .$list_genes %>%  table %>% 
  as.data.frame() %>% 
  arrange(Freq)
 

# prepare data to heatmpa for tissues
expression_data_tissue <- log(chi2_results$chi2_value_matrix[!rownames(chi2_results$p_value_matrix) %in% excluded_rows, excluded_rows[1:9]] + 1)

# Visualize the tissue-specific gene expression using heatmap
pheatmap(
  expression_data_tissue,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_row = row_mapping,
  labels_col = tissue_vector,
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)



# ---------------------
# 4. Perform Chi-squared Tests for Overlap with PRS Models from Biobank UK Phenotypes
# ---------------------

# Perform chi-squared test between gr_gene_list and filtered_phenotypes
chi2_results_models <- perform_chi2_tests(c(gr_gene_list[excluded_rows], filtered_phenotypes), hgnc_symbols_vector)

# ---------------------
# 5. Select Significant Results from the p-value Matrix
# ---------------------

# Extract significant results based on adjusted p-values (FDR)

# filter significant comparison
extract_data(chi2_results_models , !rownames(chi2_results_models $p_value_matrix) %in% excluded_rows, excluded_rows[10:27]) %>% 
  rename(phenotype = Var1, cluster = Var2) %>% 
  filter(number_overlap > 2) %>%
  # .$p_value %>% p.adjust(., method = "fdr") %>% {. < 0.05}
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr)  %>% 
  filter(fdr < 0.05) %>% 
  .$overlap_genes %>% 
  strsplit(., split = ",") %>% 
  unlist %>% 
  unique() %>% 
  walk(., ~cat(.x, "\n"))


significant_results <- extract_data(chi2_results_models , !rownames(chi2_results_models $p_value_matrix) %in% excluded_rows, excluded_rows[10:27]) %>% 
  rename(phenotype = Var1, cluster = Var2) %>% 
  filter(number_overlap > 2) %>%
  # .$p_value %>% p.adjust(., method = "fdr") %>% {. < 0.05}
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% .$phenotype %>% unique() %>% as.character()

# ---------------------
# 6. Choose Phenotypes for Visualization
# ---------------------

# List of phenotypes for heatmap visualization
selected_phenotypes <- c(
  "biobankuk-30280-both_sexes--immature_reticulocyte_fraction-EUR-1e-08.yml",
  "biobankuk-3064-both_sexes--peak_expiratory_flow_pef_-EUR-1e-08.yml",
  "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08.yml",
  "biobankuk-30630-both_sexes--apolipoprotein_a-EUR-1e-08.yml",
  "biobankuk-20004-both_sexes-1479-operation_code-EUR-1e-08.yml",
  "biobankuk-sbp-both_sexes--systolic_blood_pressure_automated_reading_adjusted_by_medication-EUR-1e-08.yml",
  "biobankuk-30620-both_sexes--alanine_aminotransferase-EUR-1e-08.yml"
)

selected_phenotypes_mapping <- c(
  "biobankuk-30280-both_sexes--immature_reticulocyte_fraction-EUR-1e-08.yml" = "Immature reticulocyte fraction",
  "biobankuk-3064-both_sexes--peak_expiratory_flow_pef_-EUR-1e-08.yml" = "Peak expiratory flow (PEF)",
  "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08.yml" = "Trciglycerides",
  "biobankuk-30630-both_sexes--apolipoprotein_a-EUR-1e-08.yml" = "Apolipoprotein A",
  "biobankuk-20004-both_sexes-1479-operation_code-EUR-1e-08.yml" = "Varicose vein surgery",
  "biobankuk-sbp-both_sexes--systolic_blood_pressure_automated_reading_adjusted_by_medication-EUR-1e-08.yml" = "Systolic blood pressure",
  "biobankuk-30620-both_sexes--alanine_aminotransferase-EUR-1e-08.yml" = "Alanine aminotransferase"
)

selected_phenotypes_mapping <- c(
  "biobankuk-30280-both_sexes--immature_reticulocyte_fraction-EUR-1e-08.yml" = "immature reticulocyte fraction",
  "biobankuk-3064-both_sexes--peak_expiratory_flow_pef_-EUR-1e-08.yml" = "peak expiratory flow",
  "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08.yml" = "triglycerides",
  "biobankuk-20004-both_sexes-1479-operation_code-EUR-1e-08.yml" = "varicose vein surgery",
  "biobankuk-30290-both_sexes--high_light_scatter_reticulocyte_percentage-EUR-1e-08.yml" = "hight light scatter reticulocyte percentage",
  "biobankuk-30300-both_sexes--high_light_scatter_reticulocyte_count-EUR-1e-08.yml" = "hight light scatter reticulocyte count",
  "biobankuk-sbp-both_sexes--systolic_blood_pressure_automated_reading_adjusted_by_medication-EUR-1e-08.yml" = "systolic blood pressure (1)",
  "biobankuk-4080-both_sexes--systolic_blood_pressure_automated_reading-EUR-1e-08.yml" = "systolic blood pressure (2)",
  "biobankuk-30620-both_sexes--alanine_aminotransferase-EUR-1e-08.yml" = "alanine aminotransferase",
  "biobankuk-30120-both_sexes--lymphocyte_count-EUR-1e-08.yml" = "lymphocyte count",
  "biobankuk-30130-both_sexes--monocyte_count-EUR-1e-08.yml" = "monocyte count",
  "biobankuk-5256-both_sexes--corneal_hysteresis_right_-EUR-1e-08.yml" = "corneal hysteris right",
  "biobankuk-30240-both_sexes--reticulocyte_percentage-EUR-1e-08.yml" = "reticulocyte percentage"
) %>% unique

cluster_vector <- c(
  "marpiech_tissues_dex_1" = "cluster 1",  
  "marpiech_tissues_dex_10" = "cluster 10", 
  "marpiech_tissues_dex_11" = "cluster 11", 
  "marpiech_tissues_dex_12" = "cluster 12", 
  "marpiech_tissues_dex_13" = "cluster 13", 
  "marpiech_tissues_dex_14" = "cluster 14",
  "marpiech_tissues_dex_15" = "cluster 15", 
  "marpiech_tissues_dex_16" = "cluster 16", 
  "marpiech_tissues_dex_17" = "cluster 17", 
  "marpiech_tissues_dex_18" = "cluster 18", 
  "marpiech_tissues_dex_2" = "cluster 2",  
  "marpiech_tissues_dex_3" = "cluster 3", 
  "marpiech_tissues_dex_4" = "cluster 4",  
  "marpiech_tissues_dex_5" = "cluster 5", 
  "marpiech_tissues_dex_6" = "cluster 6",  
  "marpiech_tissues_dex_7" = "cluster 7",  
  "marpiech_tissues_dex_8" = "cluster 8",  
  "marpiech_tissues_dex_9" = "cluster 9"
)

# ---------------------
# 7. Visualize Significant Results using Heatmap
# ---------------------

# Visualize chi-squared results for the selected phenotypes
heatmap_data <- chi2_results_models$chi2_value_matrix[unique(significant_results), names(cluster_vector)]
pheatmap(
  log2(heatmap_data +1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  # main = "Overlap with PRS Models from Biobank UK Phenotypes",
  labels_row = selected_phenotypes_mapping,
  labels_col = cluster_vector,
  fontsize = 12,             
  fontsize_row = 15,          
  fontsize_col = 15
)

# Visualize chi-squared tissue specific resuts for the significatn phenotypes

significant_results  <- extract_data(chi2_results_models , !rownames(chi2_results_models$p_value_matrix) %in% excluded_rows, excluded_rows[1:9]) %>% 
  rename(phenotype = Var1, tissue = Var2) %>% 
  filter(number_overlap > 2) %>%
  # .$p_value %>% p.adjust(., method = "fdr") %>% {. < 0.05}
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05) %>% .$phenotype %>% unique() %>% as.character()
  
  
extract_data(chi2_results_models , !rownames(chi2_results_models$p_value_matrix) %in% excluded_rows, excluded_rows[1:9]) %>% 
  rename(phenotype = Var1, tissue = Var2) %>% 
  filter(number_overlap > 2) %>%
  # .$p_value %>% p.adjust(., method = "fdr") %>% {. < 0.05}
  mutate(fdr = p.adjust(p_value, method = "fdr")) %>% 
  arrange(fdr) %>% 
  filter(fdr < 0.05)
# chi2_results$p_value_matrix %>%  
#   .[21:29, 47:ncol(chi2_results$p_value_matrix)] %>% 
#   reshape2::melt() %>% 
#   arrange(value) %>%
#   mutate(fdr = p.adjust(value, "fdr")) %>% 
#   filter(fdr < 0.05) %>% 
#   .$Var2 %>% as.character() %>% unique() -> significant_results

heatmap_data <-  chi2_results_models$chi2_value_matrix[significant_results, excluded_rows[1:9]]



phenotypes_vector <- c(
  "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08.yml" = "triglycerides",
  "biobankuk-30830-both_sexes--shbg-EUR-1e-08.yml" = "SHBG",
  "biobankuk-30630-both_sexes--apolipoprotein_a-EUR-1e-08.yml" = "apolipoprotein A",
  "biobankuk-20154-both_sexes--forced_expiratory_volume_in_1_second_fev1_predicted_percentage-EUR-1e-08.yml" = "forced expiratory volume",
  "biobankuk-30880-both_sexes--urate-EUR-1e-08.yml" = "urate"
)

phenotypes_vector <- c(
  "biobankuk-5256-both_sexes--corneal_hysteresis_right_-EUR-1e-08.yml" = "Corneal hystereis right",
  "biobankuk-30870-both_sexes--triglycerides-EUR-1e-08.yml" = "triglycerides",
  "biobankuk-30620-both_sexes--alanine_aminotransferase-EUR-1e-08.yml" = "alanine aminotransferase",
  "biobankuk-30830-both_sexes--shbg-EUR-1e-08.yml" = "SHBG"
)

# Determine the breaks for the color scale
breaks <- seq(min(log2(heatmap_data + 1), na.rm = TRUE), 
              max(log2(heatmap_data + 1), na.rm = TRUE), 
              length.out = 101)

# Create custom legend labels
legend_labels <- c(min(breaks), "log2(chi2+1) value", max(breaks))

pheatmap(
  log2(heatmap_data +1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  labels_col = tissue_vector,
  labels_row = phenotypes_vector,
  fontsize = 15,             
  fontsize_row = 15,          
  fontsize_col = 15
)

# Define the breaks and colors for the legend
breaks = seq(min(log2(heatmap_data + 1)), max(log2(heatmap_data + 1)), length.out = 5)
colors = colorRampPalette(c("blue", "white", "red"))(length(breaks) - 1)

# Plot the legend
legend("topright", legend = breaks, fill = colors, title = "Your Legend Description")


pheatmap(
  log2(heatmap_data +1),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  main = "Overlap with PRS Models from Biobank UK Phenotypes",
  annotation = list(significance = significance_matrix),
  annotation_colors = list(significance = c("" = "white", "*" = "red"))
)



# Remove variables from the first corrected code
rm(chi2_results)
rm(expression_data_clustered)
rm(expression_data_specific)
rm(tissue_data)
rm(tissue_colnames)

# Remove variables from the second corrected code
rm(significant_results)
rm(selected_phenotypes)
rm(heatmap_data)


# testing
# gr_gene_database_preproccesing <- 
gr_gene_database %>%
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
  mutate(label = ifelse(source == "marpiech_tissues_dex", paste0(source, "_", gene_list_number), label)) -> gr_gene_database_preproccesing


# Overlap of genes 
gr_gene_database_preproccesing %>% filter(source == "marpiech_tissues_dex")

split(gr_gene_database_preproccesing$hgnc_symbol, gr_gene_database_preproccesing$label) -> gr_gene_list
