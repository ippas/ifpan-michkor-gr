
# ##############################################################################
# ---- load data ----
# ##############################################################################
read.delim(gzfile("data/PGC/mergedPGC_annotation_geneCenter50kb_p1e4.tsv"),
           header = TRUE, sep = "\t") %>% 
  filter(source_file != "") -> pgc_annotation_geneCenter50kb_p1e4

# ##############################################################################
# ---- prepare_metadata ----
# ##############################################################################
pgc_annotation_geneCenter50kb_p1e4$source_file %>% 
  unique() %>% 
  as.data.frame() %>% 
  set_colnames(c("filename")) %>% 
  mutate(
    disease = case_when(
      # ADHD = Attention Deficit Hyperactivity Disorder
      str_detect(filename, regex("adhd2019", ignore_case = TRUE)) ~ "ADHD",
      str_detect(filename, regex("adhd2022", ignore_case = TRUE)) ~ "ADHD",
      str_detect(filename, regex("adhd2018", ignore_case = TRUE)) ~ "ADHD",
      str_detect(filename, regex("adhd2010", ignore_case = TRUE)) ~ "ADHD",
      str_detect(filename, regex("adhd", ignore_case = TRUE)) ~ "ADHD",
      # ALZ = Alzheimer's Disease
      str_detect(filename, regex("alz2021", ignore_case = TRUE)) ~ "ALZ",
      str_detect(filename, regex("alz2019", ignore_case = TRUE)) ~ "ALZ",
      str_detect(filename, regex("alz", ignore_case = TRUE)) ~ "ALZ",
      # ANX = Anxiety Disorder
      str_detect(filename, regex("panic2019", ignore_case = TRUE)) ~ "ANX",
      str_detect(filename, regex("anxiety", ignore_case = TRUE)) ~ "ANX",
      # ASD = Autism Spectrum Disorder
      str_detect(filename, regex("ASD", ignore_case = TRUE)) ~ "ASD",
      str_detect(filename, regex("asd2017", ignore_case = TRUE)) ~ "ASD",
      # BIP = Bipolar Disorder
      str_detect(filename, regex("bip2024", ignore_case = TRUE)) ~ "BIP",
      str_detect(filename, regex("bip2021_noUKB", ignore_case = TRUE)) ~ "BIP",
      str_detect(filename, regex("bip2021", ignore_case = TRUE)) ~ "BIP",
      str_detect(filename, regex("bip32b_mds7a_0416a", ignore_case = TRUE)) ~ "BIP",
      str_detect(filename, regex("bip.full.2012-04.txt", ignore_case = TRUE)) ~ "BIP",
      str_detect(filename, regex("daner_bip_pgc3_nm_noukbiobank", ignore_case = TRUE)) ~ "BIP",
      # CDG = Cross-Disorder
      str_detect(filename, regex("PFactor_2025.tsv.gz|F5_SubstanceUse_2025.tsv.gz|F4_Internalizing_2025.tsv.gz|F3_Neurodevelopmental_2025.tsv.gz|F2_SchizophreniaBipolar_2025.tsv.gz|F1_CompulsiveDisorders_2025.tsv.gz", ignore_case = TRUE)) ~ "CDG",
      str_detect(filename, regex("MDD_BIP_METACARPA_BD1_INFO6_A5_NTOT_Top_10K_Clumped|MDD_BIP_METACARPA_BD2_INFO6_A5_NTOT_Top_10K_Clumped|MDD_BIP_METACARPA_SAB_INFO6_A5_NTOT_Top_10K_Clumped|MDD_MHQ_BIP_METACARPA_INFO6_A5_NTOT_Top_10K_Clumped|MDD_MHQ_BIP_METACARPA_INFO6_A5_NTOT_no23andMe_noUKBB|MDD_MHQ_METACARPA_INFO6_A5_NTOT_Top_10K_Clumped|MDD_MHQ_METACARPA_INFO6_A5_NTOT_no23andMe_noUKBB|MDD_MHQ_Recurrent_METACARPA_INFO6_A5_NTOT_Top_10K_Clumped|MDD_MHQ_Single_METACARPA_INFO6_A5_NTOT_Top_10K_Clumped|MDD_MHQ_Subthreshold_METACARPA_INFO6_A5_NTOT_Top_10K_Clumped|MHQ_Depression_WG_MAF1_INFO4_HRC_Only_Filtered_Dups_FOR_METACARPA_INFO6_A5_NTOT|MHQ_Recurrent_Depression_WG_MAF1_INFO4_HRC_Only_Filtered_Dups_FOR_METACARPA_INFO6_A5_NTOT|MHQ_Single_Depression_WG_MAF1_INFO4_HRC_Only_Filtered_Dups_FOR_METACARPA_INFO6_A5_NTOT|MHQ_Subthreshold_WG_MAF1_INFO4_HRC_Only_Filtered_Dups_FOR_METACARPA_INFO6_A5_NTOT|daner_PGC_BIP32b_mds7a_mds7a_BD1.0416a_INFO6_A5_NTOT|daner_PGC_BIP32b_mds7a_mds7a_BD2.0416a_INFO6_A5_NTOT", ignore_case = TRUE)) ~ "CDG",
      str_detect(filename, regex("cdg2019", ignore_case = TRUE)) ~ "CDG",
      str_detect(filename, regex("BDSCZvsCONT|BDvsCONT|SCZvsBD|sczvscont-sumstat", ignore_case = TRUE)) ~ "CDG",
      str_detect(filename, regex("pgc.cross.full.2013", ignore_case = TRUE)) ~ "CDG",
      # ED = Eating Disorders
      str_detect(filename, regex("pgcAN2.2019-07", ignore_case = TRUE)) ~ "ED",
      str_detect(filename, regex("pgc.ed.freeze1.summarystatistics.July2017", ignore_case = TRUE)) ~ "ED",
      # MDD = Major Depressive Disorder
      str_detect(filename, regex("ppd2023", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("mdd2025", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("mdd_symptoms_2023", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("mdd2023", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("jamapsy_Giannakopoulou_2021_exclude_whi_23andMe_ukb|jamapsy_Giannakopoulou_2021_exclude_whi_23andMe", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("PGC_UKB_23andMe_depression_10000|PGC_UKB_depression_genome-wide", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("daner_pgc_mdd_meta_w2_no23andMe_rmUKBB", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("MDD2018_ex23andMe|PGC_MDD2018_10kSNPs|daner_pgc_mdd_meta_w2_no23andMe_rmUKBB", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("pgc.mdd.clump.2012-04|pgc.mdd.full.2012-04", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("mdd2018_noUKBB-PMID29700475", ignore_case = TRUE)) ~ "MDD",
      str_detect(filename, regex("AntiDepNonRemissionEAS2021|AntiDepNonRemissionEUR2021|AntiDepPercImprovEAS2021|AntiDepPercImprovEUR2021", ignore_case = TRUE)) ~ "MDD",
      # OCD-TS = OCD & Tourette Syndrome
      str_detect(filename, regex("TS_Oct2018", ignore_case = TRUE)) ~ "OCD-TS",
      str_detect(filename, regex("ocs2024", ignore_case = TRUE)) ~ "OCD-TS",
      str_detect(filename, regex("daner_OCD_full_incl23andMe_clump_W-3000_R-0p2_P-0p05_hrcALL.clumped_10000|daner_OCD_full_wo23andMe_190522|daner_OCDmeta_wo23andMe_LOOUKBB", ignore_case = TRUE)) ~ "OCD-TS",
      str_detect(filename, regex("ocd_aug2017", ignore_case = TRUE)) ~ "OCD-TS",
      str_detect(filename, regex("hoarding2022", ignore_case = TRUE)) ~ "OCD-TS",
      # OTHER, Clozapine
      str_detect(filename, regex("daner_cia_combo_eur-qc.ch.fl", ignore_case = TRUE)) ~ "CLOZAPINE",
      # PTSD = Post Traumatic Stress Disorder
      str_detect(filename, regex("aam_ptsd_pcs_v5_jan4_2022|eur_ptsd_pcs_v4_aug3_2021|eur_ptsdcasecontrol_pcs_v4_aug3_2021|hna_ptsd_pcs_v4_aug3_2021|trans_ptsd_pcs_v4_aug3_2021", ignore_case = TRUE)) ~ "PTSD",
      str_detect(filename, regex("pts_aam_freeze2_overall|pts_all_freeze2_overall|pts_eur_freeze2_overall|pts_lat_freeze2_overall", ignore_case = TRUE)) ~ "PTSD",
      str_detect(filename, regex("SORTED_PTSD_AA7_ALL_study_specific_PCs1|SORTED_PTSD_EA9_AA7_LA1_SA2_ALL_study_specific_PCs1|SORTED_PTSD_EA9_ALL_study_specific_PCs1", ignore_case = TRUE)) ~ "PTSD",
      str_detect(filename, regex("Freeze1.5_eur_samplesize_info_frq_ordinal_AND_linear_maf5_MAF_FILTER1|Freeze2_CT1_eur_samplesize_info_frq_ordinal_AND_linear_maf5_MAF_FILTER1|Trauma_One_Age_Sex_WG_MAF1_INFO4_Filtered_Dups_maf5", ignore_case = TRUE)) ~ "PTSD",
      # SCZ = Schizophrenia
      str_detect(filename, regex("PGC3_SCZ_wave3.afram.autosome.public.v3.vcf.tsv|PGC3_SCZ_wave3.asian.autosome.public.v3.vcf.tsv|PGC3_SCZ_wave3.core.autosome.public.v3.vcf.tsv|PGC3_SCZ_wave3.core.chrX.public.v3.vcf.tsv|PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv|PGC3_SCZ_wave3.latino.autosome.public.v3.vcf.tsv|PGC3_SCZ_wave3.primary.autosome.public.v3.vcf.tsv|PGC3_SCZ_wave3.primary.chrX.public.v3.vcf.tsv|daner_PGC_SCZ_w3_14_0618a_asn_female|daner_PGC_SCZ_w3_14_0618a_asn_male|daner_PGC_SCZ_w3_75_0618a_eur_female|daner_PGC_SCZ_w3_75_0618a_eur_male|daner_PGC_SCZ_w3_81_0618a_all_female|daner_PGC_SCZ_w3_81_0618a_all_male|daner_PGC_SCZ_w3_90_0418b_ukbbdedupe|daner_scz_w3_HRC_chrX_asn_fem_run2|daner_scz_w3_HRC_chrX_asn_mal_run2|daner_scz_w3_HRC_chrX_eur_fem_deduped_0518e|daner_scz_w3_HRC_chrX_eur_mal_deduped_0518e", ignore_case = TRUE)) ~ "SCZ",
      str_detect(filename, regex("scz2019asi-PMID31740837", ignore_case = TRUE)) ~ "SCZ",
      str_detect(filename, regex("scz2018clozuk-PMID29483656.tsv", ignore_case = TRUE)) ~ "SCZ",
      str_detect(filename, regex("scz2014-PMID25056061.tsv", ignore_case = TRUE)) ~ "SCZ",
      str_detect(filename, regex("scz.swe.pgc1.results.v3.txt", ignore_case = TRUE)) ~ "SCZ",
      str_detect(filename, regex("pgc.scz.clump.2012-04", ignore_case = TRUE)) ~ "SCZ",
      # SUD = Substance Use Disorders
      str_detect(filename, regex("OD_cases_vs._opioid-exposed_controls_in_African-ancestry_cohorts|OD_cases_vs._opioid-exposed_controls_in_European-ancestry_cohorts|OD_cases_vs._opioid-exposed_controls_in_the_trans-ancestry_meta-analysis|OD_cases_vs._opioid-unexposed_controls_in_African-ancestry_cohorts|OD_cases_vs._opioid-unexposed_controls_in_European-ancestry_cohorts|OD_cases_vs._opioid-unexposed_controls_in_the_trans-ancestry_meta-analysis|opioid-exposed_vs._opioid-unexposed_controls_in_African-ancestry_cohorts|opioid-exposed_vs._opioid-unexposed_controls_in_European-ancestry_cohorts|opioid-exposed_vs._opioid-unexposed_controls_in_the_trans-ancestry_meta-analysis", ignore_case = TRUE)) ~ "SUD",
      str_detect(filename, regex("CUD_AFR_full_public_11.14.2020|CUD_EUR_casecontrol_public_11.14.2020|CUD_EUR_full_public_11.14.2020", ignore_case = TRUE)) ~ "SUD",
      str_detect(filename, regex("AUDIT_UKB_2018_AJP.txt", ignore_case = TRUE)) ~ "SUD",
      str_detect(filename, regex("pgc_alcdep.afr_discovery.aug2018_release|pgc_alcdep.afr_unrel_genotyped.aug2018_release|pgc_alcdep.afr_unrelated.aug2018_release|pgc_alcdep.discovery.aug2018_release|pgc_alcdep.eur_discovery.aug2018_release|pgc_alcdep.eur_unrel_genotyped.aug2018_release|pgc_alcdep.eur_unrelated.aug2018_release|pgc_alcdep.trans_fe_unrel_geno.aug2018_release|pgc_alcdep.trans_mantra_unrel_geno.aug2018_release|pgc_alcdep.trans_re2_unrel_geno.aug2018_release", ignore_case = TRUE)) ~ "SUD",
      str_detect(filename, regex("Hatoum2023AddictionAfrican.txt|Hatoum2023AddictionEuropean.txt", ignore_case = TRUE)) ~ "SUD",
      # SUI = Suicide attempt
      
      
      TRUE ~ "Other"
    )
  ) %>% 
  mutate(
    disease_name = case_when(
      disease == "ADHD"      ~ "Attention-Deficit/Hyperactivity Disorder",
      disease == "ALZ"       ~ "Alzheimer's Disease",
      disease == "ANX"       ~ "Anxiety Disorders",
      disease == "ASD"       ~ "Autism Spectrum Disorder",
      disease == "BIP"       ~ "Bipolar Disorder",
      disease == "CDG"       ~ "Cross-Disorder Psychiatric Traits",
      disease == "ED"        ~ "Eating Disorders",
      disease == "MDD"       ~ "Major Depressive Disorder",
      disease == "OCD-TS"    ~ "Obsessive-Compulsive and Tourette Disorders",
      disease == "CLOZAPINE" ~ "Clozapine-Induced Agranulocytosis",
      disease == "PTSD"      ~ "Post-Traumatic Stress Disorder",
      disease == "SCZ"       ~ "Schizophrenia",
      disease == "SUD"       ~ "Substance Use Disorders",
      disease == "Other"     ~ "Other / Not Classified",
      TRUE ~ disease
    )
  ) %>% 
  mutate(disease_name = str_replace_all(disease_name, " ", "_")) -> pgc_metadata_diseases


# ##############################################################################
# ---- save to tsv ----
# ##############################################################################
write.table(
  pgc_metadata_diseases,
  file = "data/PGC/pgc_metadata_diseases.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
