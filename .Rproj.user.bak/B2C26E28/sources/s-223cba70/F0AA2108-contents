
# wczytanie
read_tsv("data/databases/pgc_gene_list.tsv", col_names = T) %>% 
  unique() %>%
  rename(file = "pgc_gwas") %>% 
  mutate(trait = case_when(
    grepl("14102594", pgc_gwas) ~ "BIP",
    grepl("14671971", pgc_gwas) ~ "ET",
    grepl("14671980", pgc_gwas) ~ "ED",
    grepl("14671989", pgc_gwas) ~ "ASD",
    grepl("14671998", pgc_gwas) ~ "BIP",
    grepl("14672019", pgc_gwas) ~ "CDG",
    grepl("14672040", pgc_gwas) ~ "CDG",
    grepl("14672085", pgc_gwas) ~ "MDD",
    grepl("14672133", pgc_gwas) ~ "PTSD",
    grepl("14672187", pgc_gwas) ~ "SUD",
    grepl("14842689", pgc_gwas) ~ "ANX",
    grepl("14842692", pgc_gwas) ~ "SUD",
    grepl("19193084", pgc_gwas) ~ "SCZ",
    grepl("19426775", pgc_gwas) ~ "SCZ",
    grepl("22564390", pgc_gwas) ~ "ADHD",
    grepl("22564402", pgc_gwas) ~ "BIP",
    grepl("22745573", pgc_gwas) ~ "MDD",
    grepl("24268882", pgc_gwas) ~ "SUD",
    grepl("26349322", pgc_gwas) ~ "PTSD",
    grepl("27061255", pgc_gwas) ~ "MDD",
    grepl("27216117", pgc_gwas) ~ "BIP",
    pgc_gwas == "raw/daner_adhd_meta_filtered_NA_iPSYCH23_PGC11_sigPCs_woSEX_2ell6sd_EUR_Neff_70.meta.gz" ~ "ADHD",
    pgc_gwas == "raw/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB.gz" ~ "MDD",
    pgc_gwas == "raw/CLOZUK_PGC2noclo.METAL.assoc.dosage.fix.gz" ~ "SCZ",
    pgc_gwas == "raw/daner_PGC_SCZ52_0513a.hq2.gz" ~ "SCZ",
    pgc_gwas == "raw/PGC.ASD.euro.all.25Mar2015.txt.gz" ~ "ASD",
    pgc_gwas == "raw/scz.swe.pgc1.2013-11b/scz.swe.pgc1.results.v3.txt.gz" ~ "SCZ"
  )) -> pgc_genes_df

pgc_genes_df %>% 
  select(c(gene_symbol, pgc_gwas, trait)) %>% 
  unique() %>% 
  mutate(pgc_gwas = str_replace_all(pgc_gwas, "raw/", "")) %>% 
  group_by(pgc_gwas) %>% 
  nest() %>% 
  mutate(gene_count = map(data, ~nrow(.x))) -> pgc_genes_df


categorized_gene_lists$phenotypes_PGC$metadata <- pgc_genes_df %>% 
  unnest() %>% 
  select(-gene_symbol) %>% 
  unique 


categorized_gene_lists$phenotypes_PGC$gene_lists <- pgc_genes_df %>% 
  unnest() %>% 
  select(pgc_gwas, gene_symbol) %>% 
  split(.$pgc_gwas, .$gene_symbol) %>% 
  lapply(., function(x){x$gene_symbol})

  
pgc_genes_df %>% 
  unnest() %>% 
  ungroup %>% 
  select(trait, gene_symbol) %>% 
  unique() %>% 
  group_by(trait) %>% 
  nest() %>% 
  mutate(gene_count = map(data, ~nrow(.x))) %>% 
  select(-data) %>% 
  unnest %>% 
  as.data.frame() 


categorized_gene_lists$diseases_PGC$metadata <- pgc_genes_df %>% 
  unnest() %>% 
  ungroup %>% 
  select(trait, gene_symbol) %>% 
  unique() %>% 
  group_by(trait) %>% 
  nest() %>% 
  mutate(gene_count = map(data, ~nrow(.x))) %>% 
  select(-data) %>% 
  unnest %>% 
  as.data.frame() 

categorized_gene_lists$diseases_PGC$gene_lists <- pgc_genes_df %>% 
  unnest() %>% 
  ungroup %>% 
  select(trait, gene_symbol) %>%
  split(.$trait, .$gene_symbol) %>% 
  lapply(., function(x){x$gene_symbol})
  



