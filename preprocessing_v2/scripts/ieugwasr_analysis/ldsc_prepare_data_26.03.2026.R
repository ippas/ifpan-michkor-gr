read_tsv(file = "data/ieu_open_gwas_project/metadata_iuegwas_26.03.2026.tsv") %>% 
  as.data.frame() %>% 
  filter(population == "European") %>% 
  filter(!is.na(sample_size)) %>%
  filter(!is.na(subcategory)) %>% 
  filter(sample_size > 10000) -> metadata_ieu_EURsampleSize10000

metadata_ieu_EURsampleSize10000


ldsc_results_433traits <- read_tsv(file = "data/ldsc_analysis/main_sampleSize10000noMissingSubcategory/parse_ldsc_rg_logs/rg_results_all_vs_all_parsed.tsv")

ldsc_results_clean <- ldsc_results_433traits %>%
  as.data.frame() %>%
  select(
    # --- NAJWAŻNIEJSZE: jakie fenotypy ---
    p1_id,
    p2_id,
    
    # --- WYNIKI LDSC ---
    summary_rg,
    summary_rg_se,
    summary_rg_z,
    summary_rg_p,
    
    # --- dodatkowe statystyki (opcjonalnie ale przydatne) ---
    gencov,
    gencov_se,
    gencov_intercept,
    gencov_intercept_se,
    
    summary_h2_obs,
    summary_h2_obs_se,
    summary_h2_int,
    summary_h2_int_se,
    
    # --- kontrola jakości / flagi ---
    status,
    finished_successfully,
    has_rg,
    has_gencov,
    has_p1_h2,
    has_p2_h2,
    
    # --- techniczne ---
    p1_read_snps,
    p2_read_snps,
    merge_sumstats_snps,
    valid_alleles_snps,
    total_time_sec,
    
    # --- NA SAMYM KOŃCU: ścieżki ---
    log_filename,
    log_file,
    out_prefix,
    p1_path,
    p2_path,
    summary_p1_path,
    summary_p2_path
  ) %>% 
  mutate(p2_id = str_remove(p2_id, ".sumstats.gz \\\\"))

ldsc_results_clean

metadata_ieu_EURsampleSize10000 %>% 
  select(id, trait, subcategory, category)

# przygotuj metadata (dla bezpieczeństwa)
metadata_clean <- metadata_ieu_EURsampleSize10000 %>%
  select(id, trait, subcategory, category)

ldsc_results_annotated <- ldsc_results_clean %>%
  
  # --- join dla p1 ---
  left_join(metadata_clean, by = c("p1_id" = "id")) %>%
  rename(
    p1_trait = trait,
    p1_subcategory = subcategory,
    p1_category = category
  ) %>%
  
  # --- join dla p2 ---
  left_join(metadata_clean, by = c("p2_id" = "id")) %>%
  rename(
    p2_trait = trait,
    p2_subcategory = subcategory,
    p2_category = category
  )

ldsc_results_annotated <- ldsc_results_annotated %>%
  select(
    # --- NAJWAŻNIEJSZE ---
    p1_id, p1_trait, p1_category, p1_subcategory,
    p2_id, p2_trait, p2_category, p2_subcategory,
    
    # --- WYNIKI ---
    summary_rg,
    summary_rg_se,
    summary_rg_z,
    summary_rg_p,
    
    # --- RESZTA ---
    everything()
  )


ldsc_results_annotated %>% 
  filter(p1_id == "ieu-a-1187") %>% 
  filter(summary_rg_p < 0.05)
