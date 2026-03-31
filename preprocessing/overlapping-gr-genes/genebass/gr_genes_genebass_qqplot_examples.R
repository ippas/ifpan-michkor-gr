# ##############################################################################
# ---- uses data ----
genebass_mentalHealth_skat_all %>% 
  filter(!is.na(pvalue)) -> genebass_mentalHealth_skat_all

grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up 


# ##############################################################################
# ---- brain_up examples of qqplot ----
# ##############################################################################
p1 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20518",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20518) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = T,
  gene_signature_name = "brain_up",
  lambda_gc_text_size = 6
)

p2 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20409",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20409) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = T,
  gene_signature_name = "brain_up",
  lambda_gc_text_size = 6
)

p3 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20514",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20514) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = T,
  gene_signature_name = "brain_up",
  lambda_gc_text_size = 6
)


p4 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20413",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20413) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = T,
  gene_signature_name = "brain_up",
  lambda_gc_text_size = 6
)


p1 + p2 + p3 + p4

# ---- example for brain_up with synonymous annotation ----
p5 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "synonymous"),
  phenocode = "20407",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20407) %>% 
    filter(annotation == "synonymous") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 16,
  subtitle_text_size = 16,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "brain_up",
  lambda_gc_text_size = 6
)


# ##############################################################################
# ---- examples for metasignature_down, pLoF annotation ----
# ##############################################################################
p1 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20518",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_down %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20518) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_down",
  lambda_gc_text_size = 6
)

p2 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20407",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_down %>% 
    filter(pvalue < 0.05) %>%  
    filter(phenocode == 20407) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_down",
  lambda_gc_text_size = 6
)

p3 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20408",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_down %>% 
    filter(pvalue < 0.05) %>%  
    filter(phenocode == 20408) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_down",
  lambda_gc_text_size = 6
)


p4 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20512",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_down %>% 
    filter(pvalue < 0.05) %>%  
    filter(phenocode == 20512) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_down",
  lambda_gc_text_size = 6
)

p1+p2+p3+p4

# ##############################################################################
# ---- metasignature up ----
# ##############################################################################
grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up %>% filter(pvalue < 0.05) %>% 
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  filter(fdr_0.2)

p1 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20507",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20507) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_up",
  lambda_gc_text_size = 6
)

p2 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20463",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20463) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_up",
  lambda_gc_text_size = 6
)


p3 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "synonymous"),
  phenocode = "20514",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20514) %>% 
    filter(annotation == "synonymous") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_up",
  lambda_gc_text_size = 6
)


p4 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "missense|LC"),
  phenocode = "20518",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$metasignature_up %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20518) %>% 
    filter(annotation == "missense|LC") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "metasignature_up",
  lambda_gc_text_size = 6
)

p1 + p2 + p3 + p4

# ##############################################################################
# ---- brain down ----
# ##############################################################################
grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_down %>% filter(pvalue < 0.05) %>% 
  mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
  mutate(fdr_0.2 = ifelse(pvalue < p20_pvalue, T, F)) %>% 
  head


p1 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20518",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_down %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20518) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "brain_down",
  lambda_gc_text_size = 6
)


p2 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20514",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_down %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20514) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "brain_down",
  lambda_gc_text_size = 6
)

p3 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "pLoF"),
  phenocode = "20408",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_down %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20408) %>% 
    filter(annotation == "pLoF") %>% .$gene_symbol %>% unique,
  label_box_padding = 2,
  label_point_padding = 0.8,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol", "pvalue", "annotation"),
  lambda_gc = F,
  gene_signature_name = "brain_down",
  lambda_gc_text_size = 6
)

p4 <- qqplot_genebass_by_phenocode_v2(
  genebass_data = genebass_mentalHealth_skat_all %>% filter(annotation == "synonymous"),
  phenocode = "20544",
  pvalue_column = "pvalue",
  highlight_color = "blue",
  highlight_genes = grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association$brain_down %>% filter(pvalue < 0.05) %>%  filter(phenocode == 20544) %>% 
    filter(annotation == "synonymous") %>% .$gene_symbol %>% unique,
  label_box_padding = 4,
  label_point_padding = 2,
  title_text_size = 14,
  subtitle_text_size = 14,
  axis_title_size = 16,
  axis_text_size = 16,
  label_text_size = 5,
  label_strategy = "lowest_p",
  label_fields = c("gene_symbol"),
  lambda_gc = F,
  gene_signature_name = "brain_down",
  lambda_gc_text_size = 6
)
p1 + p2 + p3 + p4
