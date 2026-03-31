get_overlap_genes <- function(signature_name) {
  
  cleanOverlapResults_grSystemicTissues_ICD10F %>% 
    filter(grepl(signature_name, grSignature)) %>% 
    filter(p_value < 0.05) %>% 
    filter(observed_overlap >= 3) %>% 
    pull(overlap_genes) %>% 
    strsplit(",") %>% 
    unlist() %>% 
    unique()
}

genes_systemic <- get_overlap_genes("systemic")
genes_neural   <- get_overlap_genes("neural")
genes_blood    <- get_overlap_genes("blood")
genes_lung     <- get_overlap_genes("lung")


prepare_hpa_df <- function(gene_list, tissue_name) {
  
  HumanProteinAtlas_rnaTissueConsensus %>% 
    filter(Organ %in% c(
      "brain",
      "bone marrow & lymphoid tissues",
      "respiratory system"
    )) %>% 
    select(Gene, gene_symbol, Organ, organ_nTPM_max) %>% 
    filter(gene_symbol %in% gene_list) %>% 
    distinct() %>% 
    rename(gene_ensembl = Gene) %>% 
    mutate(
      log_expr = log2(organ_nTPM_max + 1),
      tissue   = tissue_name
    )
}



df_plot <- bind_rows(
  prepare_hpa_df(genes_systemic, "systemic"),
  prepare_hpa_df(genes_neural,   "neural"),
  prepare_hpa_df(genes_blood,    "blood"),
  prepare_hpa_df(genes_lung,     "lung")
) %>%
  mutate(
    tissue = factor(tissue, levels = c("systemic", "neural", "blood", "lung"))
  ) %>% 
  mutate(gene_symbol = ifelse(gene_symbol == "PRODH", paste0(gene_symbol, "_", gene_ensembl), gene_symbol))



anova_results <- df_plot %>%
  group_by(tissue) %>%
  anova_test(
    dv = log_expr,
    wid = gene_symbol,
    within = Organ
  )


posthoc_results <- df_plot %>%
  group_by(tissue) %>%
  pairwise_t_test(
    log_expr ~ Organ,
    paired = TRUE,
    p.adjust.method = "BH"
  )

signature_colors <- c(
  systemic = "#335C67",
  neural   = "#E09F3E",
  blood    = "#9E2A2B",
  lung     = "#540B0E"
)

plot_box <- df_plot %>%
  
  ggplot(aes(x = Organ, y = log_expr)) +
  
  geom_boxplot(
    outlier.shape = NA,
    color = "black"
  ) +
  
  geom_jitter(
    aes(color = tissue),
    width = 0.2,
    alpha = 0.8,
    size = 2
  ) +
  
  scale_color_manual(values = signature_colors) +
  
  facet_wrap(
    ~ tissue,
    nrow = 1
  ) +
  
  theme_classic() +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  
  labs(
    x = "Organ",
    y = "log2(max nTPM + 1)",
    title = "Expression of overlap genes across organs"
  )

plot_box
