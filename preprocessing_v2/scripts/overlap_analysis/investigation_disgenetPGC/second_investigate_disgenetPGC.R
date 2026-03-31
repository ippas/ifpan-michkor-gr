intersect(
  disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia,
  pgc_geneList_10e4$pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv
)

read.delim(file = "data/databases/GWAS_Catalog/gwas-association-downloaded_2025-12-02-MONDO_0005090-withChildTraits.tsv") %>%
  filter(P.VALUE < 1e-20) %>% 
  .$MAPPED_GENE %>% 
  # strsplit(" - ") %>% 
  # unlist %>% 
  # strsplit("-") %>% 
  # unlist %>% 
  # strsplit(", ") %>% 
  unlist -> schizophrenia_GWASCAtalog


read.delim(file = "data/databases/GWAS_Catalog/gwas-association-downloaded_2025-12-02-MONDO_0005090-withChildTraits.tsv") %>%
  filter(MAPPED_GENE %in% 
           intersect(
             disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia,
             schizophrenia_GWASCAtalog
           )) %>% 
  select(MAPPED_GENE, P.VALUE) %>% arrange(MAPPED_GENE) %>% 
  group_by(MAPPED_GENE) %>% 
  slice_min(order_by = P.VALUE, n = 1) %>% 
  as.data.frame() %>% 
  arrange(P.VALUE)


read.delim(file = "data/databases/GWAS_Catalog/gwas-association-downloaded_2025-12-02-MONDO_0005090-withChildTraits.tsv") %>%
  select(MAPPED_GENE, P.VALUE) %>% 
  arrange(P.VALUE) %>% 
  head(40)

pgc_annotation_geneCenter50kb_p1e4 %>% 
  filter(gene_symbol %in% intersect(
    disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia,
    schizophrenia_GWASCAtalog
  )) %>% 
  filter(source_file == "pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv") %>% 
  group_by(gene_symbol) %>% 
  slice_min(order_by = pvalue, n = 1) %>% as.data.frame() %>% .$pvalue %>% summary


pgc_annotation_geneCenter50kb_p1e4 %>% 
  filter(gene_symbol %in% intersect(
    disgenet_mentalDisorders$geneLists_scoreMin0.5$Schizophrenia,
    pgc_geneList_10e4$pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv
  )) %>% 
  filter(source_file == "pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv") %>% 
  group_by(gene_symbol) %>% 
  slice_min(order_by = pvalue, n = 1) %>% 
  as.data.frame() %>% 
  select(gene_symbol, pvalue)



pgc_annotation_geneCenter50kb_p1e4 %>% 
  filter(source_file == "pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv") %>% 
  group_by(gene_symbol) %>% 
  slice_min(order_by = pvalue, n = 1) %>% as.data.frame() %>% .$pvalue %>% summary


pgc_annotation_geneCenter50kb_p1e4 %>% 
  filter(source_file == "pgc_PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv.tsv")  %>% 
  group_by(gene_symbol) %>% 
  slice_min(order_by = pvalue, n = 1) %>% 
  as.data.frame() %>% 
  select(gene_symbol, pvalue) %>% 
  arrange(desc(pvalue)) %>% 
  tail(100)

  
  
