# ##############################################################################
# ---- prepare gene lists ----
# ##############################################################################

papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label = paste(label, regulation, sep = "_")) %>% 
  dplyr::select(label, hgnc_symbol) %>% 
  group_by(label) %>%
  nest() %>% 
  mutate(n_genes = map(data, ~ .x %>% unique %>% nrow)) %>% 
  unnest(n_genes) %>%
  # filter(n_genes >= 10) %>% 
  # filter(n_genes < 500) %>% 
  dplyr::select(-n_genes) %>% 
  unnest(data) %>% 
  group_split() %>%
  setNames(map_chr(., ~ unique(.x$label))) %>% 
  lapply(., function(x){
    x$hgnc_symbol %>% unique()
  }) %>% 
  c(., gene_lists_all) -> AllBrain2BrainSignatures2Global



papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label = paste(label, regulation, sep = "_")) %>% 
  # dplyr::select(label, hgnc_symbol) %>% 
  group_by(label) %>%
  nest() %>% 
  mutate(n_genes = map(data, ~ .x$hgnc_symbol %>% unique %>% length)) %>% 
  unnest(n_genes) %>% 
  unnest(data) %>% 
  ungroup %>% 
  as.data.frame() -> AllBrainGeneDf

AllBrainGeneLists <- AllBrain2BrainSignatures2Global[
  !names(AllBrain2BrainSignatures2Global) %in% 
    c("metasignature_up", "metasignature_down", "brain_up", "brain_down")
]

AllBrain2BrainSignatures2Global[
  setdiff(names(AllBrain2BrainSignatures2Global),
          c("metasignature_up", "metasignature_down", "brain_up", "brain_down"))
] %>% 
  unname %>% 
  unlist
  
read.delim(file = "results/gr-signatures/gr-signatures-multi-approach-24.04.2025.tsv") %>% 
  filter(signature_derivation == "marpiech_cluster") %>% 
  select(-signature_derivation) %>% 
  filter(signature_name == "cluster_P") %>% .$hgnc_symbol
  filter(hgnc_symbol == "FKBP5")
  
c(AllBrain2BrainSignatures2Global,
  read.delim(file = "results/gr-signatures/gr-signatures-multi-approach-24.04.2025.tsv") %>% 
    filter(signature_derivation == "marpiech_cluster") %>% 
    select(-signature_derivation) %>% 
    group_by(signature_name) %>% 
    summarise(genes = list(hgnc_symbol)) %>% 
    deframe()
) -> AllBrain2BrainSignatures2GlobalMarpiechClusters


# summary brain
papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  mutate(label = paste(label, regulation, sep = "_")) %>% 
  filter(hgnc_symbol == "ABCA1")



papers_data_preprocessing %>% 
  filter(simple_tissue == "brain") %>% 
  filter(source != "pmid:28500512") %>% 
  mutate(label = paste(label, regulation, sep = "_")) %>%
  select(regulation, hgnc_symbol) %>% unique %>% group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(n_reg = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_reg) %>% filter(n_reg == 2) %>% .$hgnc_symbol -> tmp


papers_data_preprocessing %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  filter(simple_tissue == "brain") %>% 
  # filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "10days", "3weeks", "3months"))) %>% 
  # filter(source != "pmid:34362910") %>% 
  select(regulation, hgnc_symbol) %>% unique %>% group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(n_reg = map(data, ~ .x %>% nrow)) %>% 
  unnest(n_reg) %>% filter(n_reg == 2) %>% .$hgnc_symbol -> tmp


papers_data_preprocessing %>% 
  mutate(label2 = paste(label, regulation, sep = "_")) %>% 
  filter(simple_tissue == "brain") %>% 
  filter(!(time %in% c("240h", "720h", "240h_vs_720h", "672h", "10days", "3weeks", "3months"))) %>% 
  filter(hgnc_symbol %in% tmp) %>% 
    select(label2, regulation, source, hgnc_symbol) %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(n_lists = map(data, ~ .x %>% nrow)) %>% 
  mutate(n_up_lists = map(data, ~ .x %>% filter(regulation == "up") %>% nrow)) %>% 
  mutate(n_down_lists = map(data, ~ .x %>% filter(regulation == "down") %>% nrow)) %>% 
  unnest(c(n_lists, n_up_lists, n_down_lists)) -> doubleReg_brainGrGenes
  

doubleReg_brainGrGenes %>% 
  filter(n_lists > 5) %>% 
  .$n_lists %>% table



doubleReg_brainGrGenes %>% 
  filter(n_lists > 4) %>% dim
  mutate(up_down_ratio = log2(n_up_lists/n_down_lists)) %>% 
  filter(up_down_ratio  > 0)


doubleReg_brainGrGenes %>% 
  filter(n_lists < 4) %>% dim
  

