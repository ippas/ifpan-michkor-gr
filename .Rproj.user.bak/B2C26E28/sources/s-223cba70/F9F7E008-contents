# ##############################################################################
# ---- functions ----
# ##############################################################################
calculate_jaccard_similarity_matrix <- function(named_list_of_vectors, heatmap = FALSE) {
  stopifnot(is.list(named_list_of_vectors))
  stopifnot(!is.null(names(named_list_of_vectors)))
  
  ids <- names(named_list_of_vectors)
  
  # Posortuj numerycznie po nazwie (np. "factor_1", "factor_2", ...)
  extract_numeric <- function(x) as.numeric(stringr::str_extract(x, "\\d+"))
  sorted_ids <- ids[order(extract_numeric(ids))]
  
  named_list_of_vectors <- named_list_of_vectors[sorted_ids]
  
  jaccard_index <- function(x, y) {
    length(intersect(x, y)) / length(union(x, y))
  }
  
  mat <- outer(
    sorted_ids, sorted_ids,
    Vectorize(function(i, j) jaccard_index(named_list_of_vectors[[i]], named_list_of_vectors[[j]]))
  )
  
  diag(mat) <- NA  # Zamień 1 na NA na przekątnej
  
  dimnames(mat) <- list(sorted_ids, sorted_ids)
  df_mat <- as.data.frame(mat)
  
  if (heatmap) {
    library(tidyverse)
    
    df_long <- df_mat %>%
      rownames_to_column("factor1") %>%
      pivot_longer(-factor1, names_to = "factor2", values_to = "jaccard")
    
    # Zachowaj kolejność numericzną na osiach
    df_long <- df_long %>%
      mutate(
        factor1 = factor(factor1, levels = sorted_ids),
        factor2 = factor(factor2, levels = sorted_ids)
      )
    
    p <- ggplot(df_long, aes(x = factor1, y = factor2, fill = jaccard)) +
      geom_tile(color = "white") +
      scale_fill_viridis_c(name = "Jaccard index", na.value = "white") +
      coord_fixed() +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()
      )
    
    print(p)
  }
  
  return(df_mat)
}


# ##############################################################################
# ---- summary factors ----
# ##############################################################################
factors_rsidGenes_P1e4tss100kb$factor_id %>% unique() %>% length()
factors_rsidGenes_P1e4tss100kb$rsID %>% unique() %>% length()
factors_rsidGenes_P1e4tss100kb$gene_symbol %>% unique() %>% length()

factors_rsidGenes_P1e4tss100kb$gene_symbol %>% unique() %>% length()

factors_rsidGenes_P1e4tss100kb %>%
  group_by(factor_id) %>% 
  nest %>% 
  mutate(n_rsid = map(data, ~ .x$rsID %>% unique %>% length)) %>% 
  mutate(n_genes = map(data, ~ .x$gene_symbol %>% unique %>% length)) %>% 
  unnest(c(n_rsid, n_genes)) %>% 
  mutate(mean_rsIDperGene = n_rsid/n_genes) %>% 
  select(-data) %>% 
  as.data.frame() 

factors_rsidGenes_P1e4tss100kb %>%
  group_by(factor_id) %>%
  summarise(rsid_set = list(unique(rsID))) %>%
  deframe() %>% 
  unname() %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("gene_symbol", "freq")) %>% 
  arrange(desc(freq)) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("shared_by_n_factors", "n_rsID")) %>% 
  mutate(perc_rsID = n_rsID/559630*100)


factors_rsidGenes_P1e4tss100kb %>%
  group_by(factor_id) %>%
  summarise(gene_set = list(unique(gene_symbol))) %>%
  deframe() %>% 
  unname() %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("gene_symbol", "freq")) %>% 
  # filter(freq == 1) %>% 
  arrange(desc(freq)) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  as.data.frame() %>% 
  set_colnames(c("shared_by_n_factors", "n_genes")) %>% 
  mutate(percent_genes = n_genes / 17508)


# Wywołanie funkcji
jaccard_df <- calculate_jaccard_similarity_matrix(
  factors_rsidGenes_P1e4tss100kb %>%
    group_by(factor_id) %>%
    summarise(gene_set = list(unique(gene_symbol))) %>%
    deframe(), heatmap = T)

# Wynik
jaccard_df %>% 
  rownames_to_column("factor1") %>%
  pivot_longer(-factor1, names_to = "factor2", values_to = "jaccard") %>%
  drop_na() %>%
  filter(factor1 != factor2) %>%
  mutate(
    f1 = pmin(factor1, factor2),
    f2 = pmax(factor1, factor2)
  ) %>%
  distinct(f1, f2, .keep_all = TRUE) %>%
  select(factor1 = f1, factor2 = f2, jaccard) %>% 
  mutate(
    similarity_class = case_when(
      jaccard <= 0.05 ~ "none",
      jaccard <= 0.15 ~ "low",
      jaccard <= 0.30 ~ "moderate",
      jaccard <= 0.50 ~ "high",
      jaccard >  0.50 ~ "very high"
    )
  ) %>%
  filter(jaccard > 0.15) %>% 
  as.data.frame()

