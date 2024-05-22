################################################################################
databases <- c("BioPlanet_2019", "KEGG_2021_Human", "WikiPathway_2023_Human", "Elsevier_Pathway_Collection", "DisGeNET", "GO_Biological_Process_2023", "ChEA_2022", "MGI_Mammalian_Phenotype_Level_4_2021")
fdr_threshold <- 0.05
overlap_threshold <- 3


regulation_norm_master_df %>% 
  filter(regulation %in% c("up")) %>% 
  filter(rank_criterion == "log2ratio", 
         size_list == 50, 
         normalize == "n_list + 2") %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold, overlap_threshold) -> results_enrichr


results_enrichr %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% nrow
results_enrichr %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% .$Genes %>% convert_genes_to_vector(split = ";") %>% entropy
results_enrichr %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% .$Genes %>% convert_genes_to_vector(split = ";") %>% unique() %>% length()


results_enrichr %>% filter(enrichr_database == "DisGeNET") %>% nrow
results_enrichr %>% filter(enrichr_database == "DisGeNET") %>% .$Genes %>% convert_genes_to_vector(split = ";") %>% entropy
results_enrichr %>% filter(enrichr_database == "DisGeNET") %>% .$Genes %>% convert_genes_to_vector(split = ";") %>% unique() %>% length()

results_enrichr %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Up-Regulated Genes by Log2(FC)")

new_enrichr_results_log2ratio_top50_down %>% 
  filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Up-Regulated Genes by Log2(FC)") -> plot2

plot1 + plot2

regulation_master_df %>% 
  filter(size_list == 50,
         rank_criterion == "log2ratio",
         regulation == "up") %>% .$hgnc_symbol
jaccard_index(set1 = new_master_upregulated_genes, 
              set2 = regulation_master_df %>% 
                filter(size_list == 50,
                       rank_criterion == "log2ratio",
                       regulation == "up") %>% .$hgnc_symbol
)