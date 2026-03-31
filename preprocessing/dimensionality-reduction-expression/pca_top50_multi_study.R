gene_list_log2ratio_top_50 %>% head


gene_list_log2ratio_top_50 %>% 
  select(-data) %>% 
  unnest() %>% 
  ungroup %>% 
  select(c(label_regulation, simple_tissue, hgnc_symbol, treatment, log2ratio)) %>% head

pca_df <- gene_list_log2ratio_top_50 %>%
  select(-data) %>%
  unnest() %>%
  ungroup() %>%
  select(label_regulation, hgnc_symbol, log2ratio) %>%
  pivot_wider(names_from = hgnc_symbol, values_from = log2ratio) %>%
  column_to_rownames("label_regulation") %>%
  as.data.frame()


pca_df %>% dim

pca_mat <- pca_df %>%
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(across(everything(), ~replace(.x, is.infinite(.x), 0))) %>%
  select(where(~ sd(.) > 0)) 


pca_mat %>% dim

pca_mat %>% .[1:130, 1:1000] %>% 
  as.matrix

tmp <- summary(is.na(pca_mat))   # ile NA?
summary(is.infinite(as.matrix(pca_mat)))  # ile Inf lub -Inf?


pca_res <- prcomp(pca_mat, scale. = TRUE)

pca_coords <- as.data.frame(pca_res$x)
pca_coords$label <- rownames(pca_coords)

fig <- plotly::plot_ly(
  data = pca_coords,
  x = ~PC1,
  y = ~PC2,
  type = 'scatter',
  mode = 'markers',
  text = ~label,
  hoverinfo = 'text',
  marker = list(
    size = 6,
    color = 'rgba(33, 150, 243, 0.6)',
    line = list(width = 1, color = 'rgba(33, 150, 243, 1)')
  )
)

fig <- fig %>%
  plotly::layout(
    title = "Interactive PCA of GR-dependent gene signatures",
    xaxis = list(
      title = paste0("PC1 (", round(summary(pca_res)$importance[2, 1] * 100, 1), "%)"),
      range = c(-15, 15)
    ),
    yaxis = list(
      title = paste0("PC2 (", round(summary(pca_res)$importance[2, 2] * 100, 1), "%)"),
      range = c(-10, 10)
    ),
    hovermode = 'closest'
  )

fig



# factor analysis
# Macierz: wiersze = sygnatury, kolumny = geny
# Wszystkie NA już wcześniej zastąpione (np. zerami)
# Kolumny ze stałą wartością usunięte
fa_mat <- pca_mat # to, co poszło do prcomp()

# install.packages("psych")  # jeśli nie masz jeszcze
fa_res <- psych::fa(fa_mat, nfactors = 3, rotate = "varimax", fm = "ml")


loadings <- as.data.frame(unclass(fa_res$loadings))
loadings$gene <- rownames(loadings)

# Najmocniejsze geny w czynniku 1
loadings %>% arrange(desc(abs(ML1))) %>% head(15)

# Albo sortować dla ML2, ML3 itd.


scores <- as.data.frame(fa_res$scores)
scores$label <- rownames(fa_mat)

# Wykres czynnikowy: np. czynnik 1 vs 2
ggplot2::ggplot(scores, ggplot2::aes(x = ML1, y = ML2)) +
  ggplot2::geom_point() +
  ggplot2::geom_text(ggplot2::aes(label = label), size = 2.5, vjust = -0.4) +
  ggplot2::theme_minimal() +
  ggplot2::labs(title = "FA: GR response factors (ML1 vs ML2)")

