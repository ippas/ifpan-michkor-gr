regulation_master_without_norm_df %>%
  filter(size_list == 50,
         rank_criterion == "log2ratio",
         regulation == "up") %>% 
  .$hgnc_symbol %>% 
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "up",
                           size_list = 50,
                           rank_criterion = "log2ratio") -> ja_tissue_up_without_norm

regulation_norm_master_df %>% 
  filter(normalize == "n_list + 0",
         size_list == 50,
         rank_criterion == "log2ratio",
         regulation == "up") %>% 
  .$hgnc_symbol %>% 
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "up",
                           size_list = 50,
                           rank_criterion = "log2ratio") -> ja_tissue_up_norm



regulation_master_double_rank_score_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  filter(regulation == "up") %>% 
  .$hgnc_symbol %>%
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "up",
                           size_list = 50,
                           rank_criterion = "log2ratio") -> ja_tissue_up_double_rs


regulation_master_without_norm_df %>%
  filter(size_list == 50,
         rank_criterion == "log2ratio",
         regulation == "down") %>% 
  .$hgnc_symbol %>% 
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "down",
                           size_list = 50,
                           rank_criterion = "log2ratio") -> ja_tissue_down_without_norm

regulation_norm_master_df %>% 
  filter(normalize == "n_list + 0",
         size_list == 50,
         rank_criterion == "log2ratio",
         regulation == "down") %>% 
  .$hgnc_symbol %>% 
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "down",
                           size_list = 50,
                           rank_criterion = "log2ratio") -> ja_tissue_down_norm



regulation_master_double_rank_score_lists %>% 
  bind_rows(.id = "label") %>% 
  separate(label, into = c("regulation", "metric", "size_list")) %>% 
  filter(metric == "log2ratio") %>% 
  filter(regulation == "down") %>% 
  .$hgnc_symbol %>%
  jaccard_index_per_tissue(gene_list =  .,
                           tissue_data = tissue_master_df,
                           regulation = "down",
                           size_list = 50,
                           rank_criterion = "log2ratio") -> ja_tissue_down_double_rs

rbind(
  `rownames<-`(ja_tissue_up_without_norm, "Up, without normalize"),
  `rownames<-`(ja_tissue_up_norm, "Up, normalize by n lists"),
  `rownames<-`(ja_tissue_up_double_rs, "Up, double rank score")
) %>% pheatmap(cluster_rows = F, cluster_cols = F, display_numbers = TRUE, main = "Jaccard Index: master upregulated genes vs tissues lists") -> ja_heatmap_up

rbind(
  `rownames<-`(ja_tissue_down_without_norm, "Down, without normalize"),
  `rownames<-`(ja_tissue_down_norm, "Down, normalize by n lists"),
  `rownames<-`(ja_tissue_down_double_rs, "Down, double rank score")
)  %>% pheatmap(cluster_rows = F, cluster_cols = F, display_numbers = TRUE, main = "Jaccard Index: master downregulated genes vs tissues lists") -> ja_heatmap_down

ja_heatmap_up  
ja_heatmap_down


rbind(
  `rownames<-`(ja_tissue_up_without_norm, "Up, without normalize"),
  `rownames<-`(ja_tissue_up_norm, "Up, normalize by n lists"),
  `rownames<-`(ja_tissue_up_double_rs, "Up, double rank score")
) %>%  apply(., 1, function(row) {
  c(
    mean = mean(row, na.rm = TRUE),
    median = median(row, na.rm = TRUE),
    sd = sd(row, na.rm = TRUE)
  )
}) %>% t()


rbind(
  `rownames<-`(ja_tissue_down_without_norm, "Down, without normalize"),
  `rownames<-`(ja_tissue_down_norm, "Down, normalize by n lists"),
  `rownames<-`(ja_tissue_down_double_rs, "Down, double rank score")
) %>%  apply(., 1, function(row) {
  c(
    mean = mean(row, na.rm = TRUE),
    median = median(row, na.rm = TRUE),
    sd = sd(row, na.rm = TRUE)
  )
}) %>% t()

# Function to set the upper triangle to NA
lower_triangular <- function(mat) {
  mat[upper.tri(mat)] <- NA
  return(mat)
}

tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  jaccard_matrix(gene_list = .) %>% pheatmap(display_numbers = T, main = "upregulated") 



tissue_master_df %>% 
  filter(regulation == "down", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol))  %>% as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  jaccard_matrix(gene_list = .) %>% pheatmap(display_numbers = T, main = "downregulated") 


################################################################################
tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "up") %>% 
      .$hgnc_symbol})) %>% 
  jaccard_matrix(gene_list = .) %>% pheatmap(display_numbers = T, main = "Jaccard similarity index, upregulated genes") 


tissue_master_df %>% 
  filter(regulation == "down", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "down") %>% 
      .$hgnc_symbol})) %>% 
  jaccard_matrix(gene_list = .) %>% pheatmap(display_numbers = T, main = "Jaccard similarity index, downregulated genes") 



tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "up") %>% 
      .$hgnc_symbol})) %>% lapply(., unique) -> master_lists_down
   
  
  perform_chi2_tests(datasets = ., total_genes = {tissue_master_df %>% 
          filter(regulation == "up", 
                 rank_criterion == "log2ratio",
                 size_list == 50)}) -> master_lists_chi2


perform_chi2_tests(c(master_lists_down), hgnc_symbols_vector_v110) -> master_lists_chi2


tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  filter(hgnc_symbol %in% master_lists_down$master_up_double_rs) %>% as.data.frame() %>% 
  filter(hgnc_symbol %in% master_genes_from_lung_down) %>% .$hgnc_symbol %>% table

master_lists_chi2$number_overlap_matrix %>% pheatmap()



tissue_master_df %>% 
  filter(regulation == "up", 
         rank_criterion == "log2ratio",
         size_list == 50) %>% 
  ungroup(hgnc_symbol) %>%
  select(c(tissue, hgnc_symbol)) %>% 
  as.data.frame() %>% 
  split(.$tissue, .$hgnc_symbol ) %>% lapply(., function(x) {x$hgnc_symbol}) %>%  
  append(list(master_up_double_rs = {regulation_master_double_rank_score_lists %>% 
      bind_rows(.id = "label") %>% 
      separate(label, into = c("regulation", "metric", "size_list")) %>% 
      filter(metric == "log2ratio") %>% 
      filter(regulation == "down") %>% 
      .$hgnc_symbol})) %>% lapply(., unique) -> master_lists_up


master_lists_up$

