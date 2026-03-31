disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
  filter_min_vector_length(min_len = 10) %>% 
  lapply(., length) %>% 
  unname %>% unlist %>% se



pgc_geneList_10e4 %>% 
  lapply(., length) %>% 
  unname %>% unlist %>% se

disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia %>% length()

disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder %>% length()

intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia,
  pgc_geneList_10e4$pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv
)



intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder,
  pgc_geneList_10e4$pgc_jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb.txt.tsv
) %>% length()


intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder,
  pgc_geneList_10e4$pgc_jamapsy_Giannakopoulou_2021_exclude_whi_23andMe.txt.tsv
) %>% length()


intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder,
  pgc_geneList_10e4$`pgc_mdd_symptoms_2023-Clin-MDD9_death.txt.tsv`
) %>% length()


intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder,
  pgc_geneList_10e4$pgc_PGC_MDD2018_10kSNPs.tsv
) %>% length()


intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder,
  pgc_geneList_10e4$`pgc_mdd_symptoms_2023-Comm-MDD9_death.txt.tsv`
) %>% length()

intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Major_Depressive_Disorder,
  pgc_geneList_10e4$`pgc_mdd_symptoms_2023-Clin-MDD7_worthless.txt.tsv`
) %>% length()

pgc_annotation_geneCenter50kb_p1e4 %>%
  filter(pvalue < 0.0001) %>%
  select(gene_symbol, source_file) %>% 
  unique %>%
  group_by(source_file) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10)
