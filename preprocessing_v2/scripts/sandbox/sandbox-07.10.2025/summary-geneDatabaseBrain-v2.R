
tmp <- summarize_gene_dataset(df = AllBrainGeneDf)

tmp$regulation_summary$n_genes 


# remove gene list with less than 10 genes
AllBrainGeneDf %>% filter(n_genes >= 10)
# remove gene list with less than 50 genes
AllBrainGeneDf %>% filter(n_genes >= 50)
# remove gene list with less than 100 genes
AllBrainGeneDf %>% filter(n_genes >= 100)

AllBrainGeneDf %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(presentFKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>% 
  filter(presentFKBP5) %>% unnest

AllBrainGeneDf %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>% 
  filter(presentClusterP) %>% unnest 
  
AllBrainGeneDf %>% 
  filter(is.na(log2ratio))

AllBrainGeneDf %>% 
  filter(!is.na(log2ratio))

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(is.na(log2ratio))

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(!is.na(log2ratio))

AllBrainGeneDf %>% 
  filter(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))

AllBrainGeneDf %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days")))

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days")))


AllBrainGeneDf %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))) %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(presentFKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>% 
  filter(presentFKBP5)

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))) %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(presentFKBP5 = map_lgl(data, ~ any(.x$hgnc_symbol %in% "FKBP5"))) %>% 
  filter(presentFKBP5)

AllBrainGeneDf %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))) %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>% 
  filter(presentClusterP) %>% unnest 

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))) %>% 
  group_by(label) %>% 
  nest() %>% 
  mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>% 
  filter(presentClusterP) %>% unnest 

AllBrainGeneDf %>% 
  group_by(source) %>% 
  nest() %>% 
  mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>% 
  filter(presentClusterP) %>% unnest 

AllBrainGeneDf %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))) %>% 
  group_by(source) %>% 
  nest() %>% 
  mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>% 
  filter(presentClusterP) %>% unnest 

AllBrainGeneDf %>% 
  filter(n_genes >= 10) %>% 
  filter(!(time %in% c("240h", "720h", "240_vs_720h", "672h", "3months", "3weeks", "10days"))) %>% 
  group_by(source) %>% 
  nest() %>% 
  mutate(presentClusterP = map_lgl(data, ~ any(.x$hgnc_symbol %in% AllBrain2BrainSignatures2GlobalMarpiechClusters$cluster_P))) %>% 
  filter(presentClusterP) %>% unnest 

