# ##############################################################################
# ---- functions ----
# ##############################################################################

run_enrichr <- function(gene_list, database) {
  if (!requireNamespace("enrichR", quietly = TRUE)) {
    install.packages("enrichR")
  }
  library(enrichR)
  
  # sprawdzenie dostępnych baz (opcjonalne)
  available_dbs <- enrichR::listEnrichrDbs()
  
  # wykonanie wzbogacenia
  enriched <- enrichR::enrichr(gene_list, databases = database)
  
  # jeśli baza jest pojedyncza, to zwracamy bez listy
  if (length(database) == 1) {
    return(enriched[[1]])
  } else {
    return(enriched)
  }
}


# ##############################################################################
# ---- analysis ----
# ##############################################################################
grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association %>% 
  lapply(., function(x){
    x %>% 
      mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
      filter(fdr_0.1)
  }) %>% 
  bind_rows() %>% 
  filter(signature_name == "brain_up") %>% .$gene_symbol %>% unique  %>% 
  run_enrichr(gene_list = ., database = "DisGeNET") -> disgenet_mentalHelathBrainUp


grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association %>% 
  lapply(., function(x){
    x %>% 
      mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
      filter(fdr_0.1)
  }) %>% 
  bind_rows() %>% 
  filter(signature_name == "metasignature_up") %>% .$gene_symbol %>% unique  %>% 
  run_enrichr(gene_list = ., database = "DisGeNET") -> disgenet_mentalHelathMetasignatureUp


grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association %>% 
  lapply(., function(x){
    x %>% 
      mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
      filter(fdr_0.1)
  }) %>% 
  bind_rows() %>% 
  filter(signature_name == "metasignature_down") %>% .$gene_symbol %>% unique  %>% 
  run_enrichr(gene_list = ., database = "DisGeNET") -> disgenet_mentalHelathMetasignatureDown


# ##############################################################################
# ---- prepare table, fdr < 0.01, min 3 genes, OR > 1  ----
# ##############################################################################
input_sizes <- c(
  brain_up = 32,
  metasignature_up = 22,
  metasignature_down = 57
)


bind_rows(
  disgenet_mentalHelathBrainUp %>% mutate(signature_name = "brain_up"),
  disgenet_mentalHelathMetasignatureUp %>% mutate(signature_name = "metasignature_up"),
  disgenet_mentalHelathMetasignatureDown %>% mutate(signature_name = "metasignature_down")
) %>% 
  select(!c(Old.P.value, Old.Adjusted.P.value)) %>% 
  set_colnames(c("term", "overlap", "pvalue", "fdr", "odds_ratio", "combinded_score", "genes", "signature_name")) %>% 
  separate(overlap, into = c("n_genes", "total_genes"), sep = "/", convert = TRUE) %>% 
  # mutate(proportion = n_genes / total_genes) %>% 
  mutate(
    input_size = input_sizes[signature_name],
    precision = n_genes / input_size,
    recall = n_genes / total_genes,
    jaccard = n_genes / (input_size + total_genes - n_genes),
    f1_score = ifelse((precision + recall) > 0, 2 * (precision * recall) / (precision + recall), 0)
  ) %>% 
  filter(fdr < 0.1) %>% 
  filter(n_genes > 2) -> grGenesMentalHealth_enrichrDisGeNET
  
grGenesMentalHealth_enrichrDisGeNET %>% 
  filter(signature_name == "brain_up") %>% 
  filter(term %in% c(
    "Heart failure",
    "Congestive heart failure",
    "Ischemic cardiomyopathy",
    "Aortic Valve Stenosis",
    "Cardiovascular Diseases"
  ))
  

grGenesMentalHealth_enrichrDisGeNET %>% 
  filter(precision > 0.4)

c(
  "RAMP2", "OGN", "PDK4", "IGF2", "RHOU", 
  "FLNC", "HOPX", "TSC22D3", "PRG4", "TNFRSF11A"
)





library(ggplot2)
library(dplyr)

disgenet_mentalHelathMetasignatureDown %>%
  head(10) %>%
  ggplot(aes(x = reorder(Term, -log10(P.value)), 
             y = -log10(P.value))) +
  geom_col(fill = "#80b1d3") +
  coord_flip() +
  labs(
    title = "Top DisGeNET terms (mental health metasignature UP)",
    x = NULL,
    y = expression(-log[10](P.value))
  ) +
  theme_minimal() +
  expand_limits(y = max(-log10(disgenet_mentalHelathMetasignatureUp$P.value[1:10])) * 1.2) +
  scale_y_continuous(limits = c(4, NA))


disgenet_mentalHelathMetasignatureUp %>%
  head(10) %>%
  ggplot(aes(y = reorder(Term, -log10(P.value)), 
             x = -log10(P.value))) +
  geom_col(fill = "firebrick", alpha = 0.25) +
  geom_text(aes(label = Term), 
            x = 2,  # << stała pozycja tekstu
            hjust = 0, 
            color = "black", 
            size = 7) +
  coord_cartesian(xlim = c(2, 12)) +
  labs(
    title = "Top DisGeNET terms (mental health metasignature UP)",
    x = expression(-log[10](P.value)),
    y = NULL
  ) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> p1



disgenet_mentalHelathMetasignatureDown %>%
  head(10) %>%
  ggplot(aes(y = reorder(Term, -log10(P.value)), 
             x = -log10(P.value))) +
  geom_col(fill = "#005b96", alpha = 0.25) +
  geom_text(aes(label = Term), 
            x = 2,  # << stała pozycja tekstu
            hjust = 0, 
            color = "black", 
            size = 7) +
  coord_cartesian(xlim = c(2, 12)) +
  labs(
    title = "Top DisGeNET terms (mental health metasignature UP)",
    x = expression(-log[10](P.value)),
    y = NULL
  ) +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) -> p2


p3 <- disgenet_mentalHelathBrainUp %>%
  head(10) %>%
  ggplot(aes(y = reorder(Term, -log10(P.value)), 
             x = -log10(P.value))) +
  geom_col(fill = "orange", alpha = 0.25) +
  geom_text(aes(label = Term), 
            x = 2,  
            hjust = 0, 
            color = "black", 
            size = 7) +
  coord_cartesian(xlim = c(2, 12)) +
  labs(
    title = "Top DisGeNET terms (mental health metasignature UP)",
    x = expression(-log[10](P.value)),
    y = NULL
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 18),
    axis.title.x = element_text(size = 20),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )



p1 <- p1 + labs(title = NULL)
p2 <- p2 + labs(title = NULL)
p3 <- p3 + labs(title = NULL)

p2 /p1/p3

(p2 +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )) /
  (p1 +
     theme(
       axis.text.x = element_blank(),
       axis.ticks.x = element_blank(),
       axis.title.x = element_blank()
     )) /
  p3


svg("data/genebass/figures/grSignatures_mentalHelath_Disgenet.svg", width = 10, height = 12)
p2/p1/p3
dev.off()

