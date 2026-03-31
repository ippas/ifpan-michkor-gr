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
  run_enrichr(gene_list = ., database = "GO_Biological_Process_2025") -> disgenet_mentalHelathBrainUp


grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association %>% 
  lapply(., function(x){
    x %>% 
      mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
      filter(fdr_0.1)
  }) %>% 
  bind_rows() %>% 
  filter(signature_name == "metasignature_up") %>% .$gene_symbol %>% unique  %>% 
  run_enrichr(gene_list = ., database = "GO_Biological_Process_2025") -> disgenet_mentalHelathMetasignatureUp


grSignature_genebass_associations_FDRMonteCarlo$original_genebass_association %>% 
  lapply(., function(x){
    x %>% 
      mutate(fdr_0.1 = ifelse(pvalue < p10_pvalue, T, F)) %>% 
      filter(fdr_0.1)
  }) %>% 
  bind_rows() %>% 
  filter(signature_name == "metasignature_down") %>% .$gene_symbol %>% unique  %>% 
  run_enrichr(gene_list = ., database = "GO_Biological_Process_2025") -> disgenet_mentalHelathMetasignatureDown







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
  coord_cartesian(xlim = c(2, 6)) +
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
  coord_cartesian(xlim = c(2, 6)) +
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
  coord_cartesian(xlim = c(2, 6)) +
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

