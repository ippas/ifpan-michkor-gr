install.packages("enrichR")
require(enrichR)

dbs <- listEnrichrDbs()

regulation_master_df %>% 
  filter(regulation == "up") %>% 
  filter(rank_criterion == "fdr", 
         size_list == 50) %>% .$hgnc_symbol %>% cat(sep = "\n")

enrichr_pathway_names <- c("BioPlanet_2019", "KEGG_2021_Human", "WikiPathway_2023_Human", "Elsevier_Pathway_Collection")
enrichr_diseases_drugs_names <- c("DisGeNET")
enrichr_ontologies <- c("GO_Biological_Process_2023", "MGI_Mammalian_Phenotype_Level_4_2021")
enrichr_transcription <- c("ChEA_2022")


################################################################################
databases <- c("BioPlanet_2019", "KEGG_2021_Human", "WikiPathway_2023_Human", "Elsevier_Pathway_Collection", "DisGeNET", "GO_Biological_Process_2023", "ChEA_2022", "MGI_Mammalian_Phenotype_Level_4_2021")
fdr_threshold <- 0.05
overlap_threshold <- 3


regulation_master_df %>% 
  filter(regulation == "up") %>% 
  filter(rank_criterion == "log2ratio", 
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold, overlap_threshold) -> enrichr_results_log2ratio_top50_up


write_tsv_xlsx(data = enrichr_results_log2ratio_top50_up,
               tsv_file = "results/google-drive/enrichr/enrichr-top50-log2ratio-up.tsv")

regulation_master_df %>% 
  filter(regulation == "down") %>% 
  filter(rank_criterion == "log2ratio", 
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold, overlap_threshold) -> enrichr_results_log2ratio_top50_down

write_tsv_xlsx(data = enrichr_results_log2ratio_top50_down,
               tsv_file = "results/google-drive/enrichr/enrichr-top50-log2ratio-down.tsv")



regulation_master_df %>% 
  filter(regulation == "up") %>% 
  filter(rank_criterion == "fdr", 
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold, overlap_threshold) -> enrichr_results_fdr_top50_up


write_tsv_xlsx(data = enrichr_results_fdr_top50_up,
               tsv_file = "results/google-drive/enrichr/enrichr-top50-fdr-up.tsv")

regulation_master_df %>% 
  filter(regulation == "down") %>% 
  filter(rank_criterion == "fdr", 
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold, overlap_threshold) -> enrichr_results_fdr_top50_down

write_tsv_xlsx(data = enrichr_results_fdr_top50_down,
               tsv_file = "results/google-drive/enrichr/enrichr-top50-fdr-down.tsv")

################################################################################
tissue_master_df %>% 
  filter(regulation == "down") %>% 
  filter(rank_criterion == "log2ratio", 
         tissue == "brain",
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold, overlap_threshold) -> enrichr_results_log2ratio_top50_down_brain


tissue_master_df %>% 
  filter(regulation == "up") %>% 
  filter(rank_criterion == "log2ratio", 
         tissue == "blood",
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold = 0.05, overlap_threshold = 3) -> enrichr_results_log2ratio_top50_up_blood

tissue_master_df %>% 
  filter(regulation == "up") %>% 
  filter(rank_criterion == "log2ratio", 
         tissue == "liver",
         size_list == 50) %>% .$hgnc_symbol %>% 
  perform_enrichment_analysis(gene_list = ., databases, fdr_threshold = 0.05, overlap_threshold = 3) -> enrichr_results_log2ratio_top50_up_lung



################################################################################
# visualization


enrichr_results_log2ratio_top50_up %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Up-Regulated Genes by Log2(FC)") -> plot1


enrichr_results_fdr_top50_up %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Up-Regulated Genes by FDR") -> plot2

plot1 + plot2

enrichr_results_log2ratio_top50_down %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Down-Regulated Genes by Log2(FC)") -> plot3


enrichr_results_fdr_top50_down %>% filter(enrichr_database == "MGI_Mammalian_Phenotype_Level_4_2021") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Down-Regulated Genes by FDR") -> plot4

plot1 + plot2 + plot3 + plot4


################################################################################
enrichr_results_log2ratio_top50_up %>% filter(enrichr_database == "DisGeNET") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Up-Regulated Genes by Log2(FC)", xlab = "DisGeNET") -> plot1


enrichr_results_fdr_top50_up %>% filter(enrichr_database == "DisGeNET") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Up-Regulated Genes by FDR", xlab = "DisGeNET") -> plot2

plot1 + plot2

enrichr_results_log2ratio_top50_down %>% filter(enrichr_database == "DisGeNET") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Down-Regulated Genes by Log2(FC)", xlab = "DisGeNET") -> plot3


enrichr_results_fdr_top50_down %>% filter(enrichr_database == "DisGeNET") %>% 
  plotEnrich(., showTerms = 40, numChar = 50, y = "Count", orderBy = "P.value", title = "Ranked Master Down-Regulated Genes by FDR", xlab = "DisGeNET") -> plot4

plot1 + plot2 + plot3 + plot4
