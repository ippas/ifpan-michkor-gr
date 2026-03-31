# ##############################################################################
# ---- uses data ----
# ##############################################################################

# genebass
genebass_online_mental %>% 
  filter(pvalue < 0.01) %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>% 
  select(c(gene_symbol, description_format)) %>% 
  unique() %>% 
  group_by(description_format) %>% 
  nest %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>%
  select(-n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% 
  distinct() %>%
  group_by(description_format) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe() -> genebass_geneLists

genebass_geneLists %>% 
  unname() %>% 
  unlist %>% 
  unique -> genebass_uniqueGenes
  

# DisGeNET
disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
  filter_min_vector_length(min_len = 10) -> DisGeNET_geneLists

disgenet_mentalDisorders$geneLists_scoreMin0.5 %>%
  filter_min_vector_length(min_len = 10) %>% 
  unname() %>% 
  unlist %>% 
  unlist %>% 
  unique() -> digenet_uniqueGenes

# PGC
pgc_annotation_geneCenter50kb_p1e4 %>%
  filter(pvalue < 0.0001) %>%
  select(gene_symbol, source_file) %>% 
  filter(gene_symbol %in% hgnc_symbols_vector_v110) %>% 
  unique %>%
  group_by(source_file) %>% 
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>% 
  select(-n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% 
  distinct() %>%
  group_by(source_file) %>%
  summarise(genes = list(unique(gene_symbol)), .groups = "drop") %>%
  deframe() -> PGC_geneLists

PGC_geneLists %>% 
  unname() %>% 
  unlist() %>% 
  unique() -> PGC_uniqueGenes

# GWAS Catalog
moodDisorders_GWASCatalog %>% 
  select(mappedGenes, efoTraits, pValue) %>% 
  separate_rows(mappedGenes, sep = ",") %>% 
  filter(mappedGenes %in% hgnc_symbols_vector_v110) %>% 
  group_by(efoTraits) %>% 
  nest %>% 
  mutate(n_genes = map(data, ~ .x %>% .$mappedGenes %>% length)) %>% 
  unnest(n_genes) %>% 
  filter(n_genes >= 10) %>% 
  unnest() %>% 
  select(c(efoTraits, mappedGenes)) %>% 
  { split(.$mappedGenes, .$efoTraits) } %>% 
  lapply(unique) -> GWASCatalog_geneLists

GWASCatalog_geneLists %>% 
  unname %>% 
  unlist %>% 
  unique() -> GWASCatalog_uniqueGenes

PGC_uniqueGenes
digenet_uniqueGenes
genebass_uniqueGenes
GWASCatalog_uniqueGenes

# comparison
with(
  list(B = unique(PGC_uniqueGenes),
       A = unique(digenet_uniqueGenes),
       C = unique(genebass_uniqueGenes)),
  c(
    A_only  = length(setdiff(A, union(B, C))),
    B_only  = length(setdiff(B, union(A, C))),
    C_only  = length(setdiff(C, union(A, B))),
    AB_only = length(setdiff(intersect(A, B), C)),
    AC_only = length(setdiff(intersect(A, C), B)),
    BC_only = length(setdiff(intersect(B, C), A)),
    ABC     = length(Reduce(intersect, list(A, B, C)))
  )
)


genebass_geneLists %>% length()
PGC_geneLists %>% length()
DisGeNET_geneLists %>% length()
GWASCatalog_geneLists %>% length()

genebass_geneLists %>% unname %>% unlist %>% unique %>%  length()
PGC_geneLists %>% unname %>% unlist %>% unique %>%  length()
DisGeNET_geneLists %>% unname %>% unlist %>% unique %>%  length()
GWASCatalog_uniqueGenes %>% length()


DisGeNET_geneLists %>% lapply(length) %>% unname() %>% unlist %>% mean
PGC_geneLists %>% lapply(length) %>% unname() %>% unlist %>% mean
genebass_geneLists %>% lapply(length) %>% unname() %>% unlist %>% mean
GWASCatalog_geneLists %>% lapply(length) %>% unname() %>% unlist %>% mean


DisGeNET_geneLists %>% lapply(length) %>% unname() %>% unlist %>% sd
PGC_geneLists %>% lapply(length) %>% unname() %>% unlist %>% sd
genebass_geneLists %>% lapply(length) %>% unname() %>% unlist %>% sd
GWASCatalog_geneLists %>% lapply(length) %>% unname() %>% unlist %>% sd


with(
  list(
    B = unique(PGC_uniqueGenes),
    A = unique(digenet_uniqueGenes),
    D = unique(genebass_uniqueGenes),
    C = unique(GWASCatalog_uniqueGenes)
  ),
  c(
    A_only = length(setdiff(A, union(union(B, C), D))),
    B_only = length(setdiff(B, union(union(A, C), D))),
    C_only = length(setdiff(C, union(union(A, B), D))),
    D_only = length(setdiff(D, union(union(A, B), C))),
    
    AB_only = length(setdiff(intersect(A, B), union(C, D))),
    AC_only = length(setdiff(intersect(A, C), union(B, D))),
    AD_only = length(setdiff(intersect(A, D), union(B, C))),
    BC_only = length(setdiff(intersect(B, C), union(A, D))),
    BD_only = length(setdiff(intersect(B, D), union(A, C))),
    CD_only = length(setdiff(intersect(C, D), union(A, B))),
    
    ABC_only = length(setdiff(Reduce(intersect, list(A, B, C)), D)),
    ABD_only = length(setdiff(Reduce(intersect, list(A, B, D)), C)),
    ACD_only = length(setdiff(Reduce(intersect, list(A, C, D)), B)),
    BCD_only = length(setdiff(Reduce(intersect, list(B, C, D)), A)),
    
    ABCD = length(Reduce(intersect, list(A, B, C, D)))
  )
)




# comparison
with(
  list(B = unique(GWASCatalog_uniqueGenes),
       A = unique(digenet_uniqueGenes),
       C = unique(genebass_uniqueGenes)),
  c(
    A_only  = length(setdiff(A, union(B, C))),
    B_only  = length(setdiff(B, union(A, C))),
    C_only  = length(setdiff(C, union(A, B))),
    AB_only = length(setdiff(intersect(A, B), C)),
    AC_only = length(setdiff(intersect(A, C), B)),
    BC_only = length(setdiff(intersect(B, C), A)),
    ABC     = length(Reduce(intersect, list(A, B, C)))
  )
)


# schizophrenia
intersect(GWASCatalog_geneLists$schizophrenia, DisGeNET_geneLists$Schizophrenia) 

# MDD
GWASCatalog_geneLists$`major depressive disorder` %>% length()
DisGeNET_geneLists$Major_Depressive_Disorder %>% length()
intersect(GWASCatalog_geneLists$`major depressive disorder`, DisGeNET_geneLists$Major_Depressive_Disorder) %>% length()

# bipolar disorder

GWASCatalog_geneLists$`bipolar disorder` %>% length()
DisGeNET_geneLists$Bipolar_Disorder %>% length()
intersect(GWASCatalog_geneLists$`bipolar disorder`, DisGeNET_geneLists$Bipolar_Disorder) %>% length()

GWASCatalog_geneLists$`autism spectrum disorder` %>% length()
DisGeNET_geneLists$Autistic_Disorder %>% length()
intersect(GWASCatalog_geneLists$`autism spectrum disorder`,
          DisGeNET_geneLists$Autistic_Disorder) %>% length()


GWASCatalog_geneLists$`Alzheimer disease` %>% length()
DisGeNET_geneLists$`Alzheimer's_Disease` %>% length()
intersect(GWASCatalog_geneLists$`Alzheimer disease`,
          DisGeNET_geneLists$`Alzheimer's_Disease`) %>% length()


GWASCatalog_geneLists$`anxiety disorder` %>% length()
DisGeNET_geneLists$Anxiety_Disorders %>% length()
intersect(GWASCatalog_geneLists$`anxiety disorder`,
          DisGeNET_geneLists$Anxiety_Disorders) %>% length()

GWASCatalog_geneLists$`panic disorder` %>% length()
DisGeNET_geneLists$Panic_Disorder %>% length()
intersect(GWASCatalog_geneLists$`panic disorder`,
          DisGeNET_geneLists$Panic_Disorder) %>% length()

GWASCatalog_geneLists$`attention deficit hyperactivity disorder` %>% length()
DisGeNET_geneLists$sleep %>% length()
intersect(GWASCatalog_geneLists$`panic disorder`,
          DisGeNET_geneLists$Panic_Disorder) %>% length()

# check MDD in GTEx
MDD_GWASCatalog_GTEx <- multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression(gene_symbol = GWASCatalog_geneLists$`major depressive disorder`)

MDD_GWASCatalog_GTEx %>% 
  filter(tissue == "Lung") %>% 
  select(-c(median_mean, median_median, median_sd, median_min, error_msg)) %>% 
  filter(median_max > 1 ) %>% 
  filter(query_gene_symbol %in% flat_allGrSignatures_31.10.2025$minusGlobalUpDown5TissuesDerivedCells_LungCellsDown)


MDD_GWASCatalog_GTEx %>% 
  filter(tissue == "Whole_Blood") %>% 
  select(-c(median_mean, median_median, median_sd, median_min, error_msg)) %>% 
  filter(median_max > 1 ) %>% dim

MDD_GWASCatalog_GTEx %>% 
  filter(median_max > 1) %>% 
  distinct(tissue, query_gene_symbol) %>% 
  count(tissue, name = "n_genes") %>% as.data.frame() %>% 
  arrange(n_genes)

sig_names
