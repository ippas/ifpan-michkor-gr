# ##############################################################################
# ---- uses data ----
# ##############################################################################

flat_allGrSignatures_17.10.2025 %>% class

flat_allGrSignatures_17.10.2025 %>% 
  lapply(length) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("signature") %>% 
  dplyr::rename(n_genes = 2) %>% 
  dplyr::mutate(
    regulation = ifelse(grepl("Up", signature), "up", "down")
  ) %>% 
  filter(regulation == "down")



flat_allGrSignatures_17.10.2025 %>% 
  lapply(length) %>% 
  unlist() %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("signature") %>% 
  dplyr::rename(n_genes = 2) %>% 
  dplyr::mutate(
    regulation = ifelse(grepl("Up", signature), "up", "down")
  ) %>% 
  dplyr::filter(regulation == "up") %>%
  ggplot(aes(x = reorder(signature, n_genes), y = n_genes)) +
  geom_col(fill = "firebrick4") +
  coord_flip() +
  theme_minimal(base_size = 13) +
  labs(
    x = "Signature",
    y = "Number of genes",
    title = "Down-regulated GR-dependent gene signatures"
  )


