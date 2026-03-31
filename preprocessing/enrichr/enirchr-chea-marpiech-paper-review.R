chea_list2 <-
  enrichr_multiple_list_and_databases(data = gr_database_blocked_gene_lists$marpiech_cluster_dex_letters,
                                      databases = c("ChEA_2022"))


gtex_tissue_gene_expression <- read.delim(
  "data/databases/GTEx/gene-expression-selected-tissues-threshold1.tsv",
  sep = "\t",
  header = TRUE
)

gtex_tissue_gene_expression %>% head

ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

nuclear_receptors <- biomaRt::getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id", "description", "external_synonym"),
  filters = "go",
  values = "GO:0004879", # nuclear receptor activity
  mart = ensembl
)

plot_df <- chea_df_raw %>%
  mutate(
    TF = stringr::str_split(Term, " ", simplify = TRUE)[, 1],
    logP = -log10(Adjusted.P.value)
  ) %>%
  group_by(cluster) %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 10) %>%
  mutate(TF_rank = paste0("top", row_number())) %>%
  ungroup() %>% 
  mutate(
    nuclear_receptors = TF %in% nuclear_receptors$hgnc_symbol,
    logP = ifelse(logP > 10, 10, logP)
  ) %>% 
  filter(Adjusted.P.value < 0.2)

plot_df$tissue_agree10 <- FALSE
plot_df$tissue_agree5 <- FALSE
plot_df$tissue_agree1 <- FALSE


plot_df %>% 
  filter(Adjusted.P.value < 0.2) %>% 
  filter(cluster == "cluster_M") %>% .$TF %>% sort %>% unique

gtex_tissue_gene_expression %>% 
  filter(grepl("JUN", gene_name)) %>% 
  filter(mean_expression > 1)

gtex_tissue_gene_expression %>% 
  filter(mean_expression > 10) %>% 
  filter(gene_name %in% {
    plot_df %>% 
      filter(Adjusted.P.value < 0.2) %>% 
      filter(cluster == "cluster_M") %>% .$TF %>% c(.)}) %>% 
  filter(label %in% c("adrenal_gland", "liver", "kidney", "pituitary", "adipose_visceral_omentum",
                      "brain_hypothalamus")) %>%
  # filter(label %in% c("adrenal_gland", "adipose_visceral_omentum", "kidney", "liver", "pituitary")) %>%
  .$gene_name %>% unique() %>% sort

# # Klaster A = ADR, FAT, LIV, MUS
# klaster B  = ADR, FAT, LIV, LUN
# Klaster C = ADR, FAT, 
# klaster E = KID
# klaster F = ADR, FAT, KID, LIV, PIT
# klaster G = SPL
# klaster H = LIV
# klaster J = ADR, FAT, HTH, KID, LIV, MUS, PIT
# klaster L =  KID, LIV, PIT
# klaster M = ADR, FAT, HTH, KID, LIV, MUS, SPL, LUN
# klaster N = HTH, KID, LIV, MUS, PIT

# cluster_D
# SMRT -> NCOR2
# SMAD2/3 -> SMAD2 i SMAD3
# RACK7 -> ZMYND8
# RXR -> RXRA


plot_df %>% 
  # cluster I
  mutate(tissue_agree10 = ifelse(cluster == "cluster_D" & TF %in% c("RELA", "IRF8", "SMAD2/3", "SMRT"), TRUE, tissue_agree10)) %>% 
  mutate(tissue_agree5 = ifelse(cluster == "cluster_D" & TF %in% c("RELA", "IRF8", "SMAD2/3", "SOX2", "NFKB1", "KDM2B", "SMRT"), TRUE, tissue_agree5)) %>% 
  mutate(tissue_agree1 = ifelse(cluster == "cluster_D" & TF %in% c("SMAD2/3","SOX2", "RELA", "KDM2B", "NFKB1", "STAT4", "IRF8", "RELB", "GATA3", "SMRT"), TRUE, tissue_agree1)) %>% 
  # cluster O
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_O" &
      TF %in% c("PPARG"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_O" &
      TF %in% c("GATA2", "NR3C1", "PPARG"),
    TRUE,
    tissue_agree5
  # )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_O" &
      TF %in% c("EZH2",
                "FOXA2",
                "GATA2",
                "NR3C1",
                "PPARG",
                "RNF2",
                "RUNX2",
                "SUZ12",
                "TP53",
                "TP53"),
    TRUE,
    tissue_agree1
  )) %>%
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_P", 
    FALSE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_P" &
      TF %in% c("AR", "NR3C1"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_P" &
      TF %in% c("AR", "CLOCK", "FOXA1", "NR3C1", "TAL1"),
    TRUE,
    tissue_agree1
  )) %>%
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_K" &
      TF %in% c("SMARCA4", "RXR", "RACK7"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_K" &
      TF %in% c("KDM2B", "NR3C1", "SMARCA4", "RXR", "RACK7"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_K" &
      TF %in% c("CLOCK", "KDM2B", "NR3C1", "NR3C2", "SMARCA4", "VDR", "RXR", "RACK7"),
    TRUE,
    tissue_agree1
  )) %>%
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_I" &
      TF %in% c("FOXO1", "PPARA", "PPARG", "RXR"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_I" &
      TF %in% c("FOXO1", "NR3C1", "PPARA", "PPARG", "LXR", "RXR"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_I" &
      TF %in% c("ESR1", "FOXA1", "FOXO1", "NR3C1", "PPARA", "PPARG", "LXR", "RXR"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "UP" &
      TF %in% c("PPARA", "RXR"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "UP" &
      TF %in% c("NR1I2", "NR3C1", "PPARA", "LXR", "RXR"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "UP" &
      TF %in% c("CLOCK", "ESR1", "FOXA2", "NR1I2", "NR3C1", "PPARA", "LXR", "RXR"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "DOWN" &
      TF %in% c("FOXO3", "IRF8", "NFE2L2", "THRA", "SMRT"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "DOWN" &
      TF %in% c("FOXO3", "IRF8", "NFE2L2", "THRA", "SMRT"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "DOWN" &
      TF %in% c("FOXO3", "IRF8", "NFE2L2", "RELB", "SUZ12", "THRA", "SMRT"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_G" &
      TF %in% c("ARID1A", "STAT1"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_G" &
      TF %in%  c("ARID1A", "SETDB1", "STAT1"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_G" &
      TF %in% c("ARID1A", "CREB1", "GATA3", "NR3C1", "PPARG", "RUNX2", "SETDB1", "STAT1", "SUZ12"),
    TRUE,
    tissue_agree1
  )) %>%   
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_A" &
      TF %in% c( "PPARA", "RXR",  "STAT3"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_A" &
      TF %in%   c("LXR", "PPARA", "RXR", "STAT3", "THRA"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_A" &
      TF %in% c("KDM2B", "LMO2", "LXR", "PPARA", "PPARD", "RXR", "STAT3", "THRA"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_B" &
      TF %in% c( "HTT"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_B" &
      TF %in%   c("AR", "CNOT3", "ERG", "GATA2", "HTT", "SMC1"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_B" &
      TF %in% c("AR", "CNOT3", "DNAJC2", "ERG", "GATA2", "HTT", "SMC1", "SMC3"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_C" &
      TF %in% c( ""),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_C" &
      TF %in% c(""),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_C" &
      TF %in% c("SUZ12"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_F" &
      TF %in% c( ""),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_F" &
      TF %in% c(""),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_F" &
      TF %in% c("ESR2"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_H" &
      TF %in% c("CEBPB", "CJUN", "PPARA", "RXR"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_H" &
      TF %in% c("CEBPB", "CJUN", "NR1H3", "PPARA", "RXR", "LXR"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_H" &
      TF %in%   c("ARNT", "CEBPB", "ESR1", "FOXO1", "CJUN", "NR1H3", "PPARA", "RXR", "LXR"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_L" &
      TF %in% c("CEBPA", "CEBPB", "EGR1", "HNF4A", "PPARA", "RXR"),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_L" &
      TF %in%  c("CEBPA", "CEBPB", "EGR1", "HNF4A", "LXR", "NR1I2", "PPARA", "RXR"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_L" &
      TF %in% c("CEBPA", "CEBPB", "EGR1", "ESR1", "FOXO1", "HNF4A", "LXR", "NR1I2", "PPARA", "RXR"),
    TRUE,
    tissue_agree1
  )) %>% 
  mutate(tissue_agree10 = ifelse(
    cluster == "cluster_M" &
      TF %in% c(""),
    TRUE,
    tissue_agree10
  )) %>%
  mutate(tissue_agree5 = ifelse(
    cluster == "cluster_M" &
      TF %in%  c("SOX2"),
    TRUE,
    tissue_agree5
  )) %>%
  mutate(tissue_agree1 = ifelse(
    cluster == "cluster_M" &
      TF %in% c("SOX2"),
    TRUE,
    tissue_agree1
  ))-> plot_df
  
  



ggplot(plot_df, aes(x = cluster, y = TF_rank, fill = logP)) +
  geom_tile(color = "black", width = 0.92, height = 0.92) +
  
  # 🔽 Dodaj fioletowy kwadracik dla receptorów jądrowych
  geom_point(
    data = subset(plot_df, nuclear_receptors == TRUE),
    aes(x = cluster, y = TF_rank),
    shape = 15, size = 6, color = "#724E91",
    position = position_nudge(y = -0.2, x = -0.1)
  ) +
  geom_point(
    data = subset(plot_df, TF == "NR3C1"),
    aes(x = cluster, y = TF_rank),
    shape = 15, size = 6, color = "#451F55",
    position = position_nudge(y = -0.2, x = -0.1)
  ) +
  geom_point(
    data = subset(plot_df, tissue_agree5 == TRUE),
    aes(x = cluster, y = TF_rank),
    shape = 15, size = 6, color = "#69995D",
    position = position_nudge(y = -0.2, x = 0.12)
  ) +
  
  geom_text(aes(label = TF), size = 5, color = "black",
            position = position_nudge(y = 0.1, x = 0)
            ) +
  scale_y_discrete(limits = paste0("top", 10:1)) +
  scale_x_discrete(position = "top") +
  scale_fill_gradient(
    low = "#ffffff", high = "#b2182b",
    name = expression(-log[10](italic(p)))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))



# Uzupełnienie siatki, by zachować puste klastry
all_clusters <- unique(chea_df_raw$cluster)
complete_df <- expand.grid(
  cluster = all_clusters,
  TF_rank = paste0("top", 1:10)
)

plot_df_complete <- complete_df %>%
  left_join(plot_df, by = c("cluster", "TF_rank"))

# Finalny wykres
ggplot(plot_df_complete, aes(x = cluster, y = TF_rank, fill = logP)) +
  geom_tile(width = 0.95, height = 0.95) +
  
  geom_point(
    data = subset(plot_df_complete, nuclear_receptors == TRUE),
    aes(x = cluster, y = TF_rank),
    # shape = 15, size = 6, color = "#724E91",
    shape = 15, size = 6, color = "#ab87c9",
    position = position_nudge(y = -0.2, x = -0.1),
    na.rm = TRUE
  ) +
  geom_point(
    data = subset(plot_df_complete, TF == "NR3C1"),
    aes(x = cluster, y = TF_rank),
    shape = 15, size = 6, color = "#451F55",
    position = position_nudge(y = -0.2, x = -0.1),
    na.rm = TRUE
  ) +
  geom_point(
    data = subset(plot_df_complete, tissue_agree5 == TRUE),
    aes(x = cluster, y = TF_rank),
    shape = 15, size = 6, color = "#69995D",
    position = position_nudge(y = -0.2, x = 0.12),
    na.rm = TRUE
  ) +
  
  geom_text(aes(label = TF), size = 5, color = "black",
            position = position_nudge(y = 0.1, x = 0),
            na.rm = TRUE
  ) +
  scale_y_discrete(limits = paste0("top", 10:1)) +
  scale_x_discrete(position = "top", limits = sort(all_clusters)) +
  
  # 🔥 KLUCZOWA ZMIANA TUTAJ 🔥 #
  scale_fill_gradient(
    low = "#ffffff", high = "#b2182b",
    limits = c(0, 10), # Skala od 0 do 10
    name = expression(-log[10](italic(p))),
    na.value = "white"
  ) +
  
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
