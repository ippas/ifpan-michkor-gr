
compute_rank_score_log2ratioPlusFDR(data = tissue_master_lists_noReg) -> tissue_master_lists_noReg_rsLog2ratioPlusFDR

tissue_master_lists_noReg_rsLog2ratioPlusFDR %>% 
  bind_rows(.id = "name") %>% 
  separate(name, into = c("tissue", "size_list", "rank_criterion"), sep = "_") %>% 
  write_tsv_xlsx(., tsv_file = "results/google-drive/master-tissue-lists/tissue-master-lits-noREG-rsLog2ratioPlusFDR.tsv")

