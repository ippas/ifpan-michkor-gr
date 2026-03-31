icd_df <- data.frame(
  icdcategory = c("F0", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9", "neurodegenerative"),
  n_associations = c(4, 7, 4, 17, 1, 0, 0, 0, 3, 0, 9),
  n_genes = c(40, 44, 39, 58, 10, 0, 0, 0, 15, 0, 65),
  n_sources = c(3, 2, 1, 2, 1, 0, 0, 0, 1, 0, 2),
  stringsAsFactors = FALSE
)


# ---- prepare data ----
# ---- prepare data ----
df_plot <- icd_df %>% 
  mutate(
    n_phenotypesAll = c(24, 44, 27, 70, 27, 12, 4, 11, 17, 8, 29),
    n_genesAll      = c(4474, 6157, 5555, 7558, 4973, 3295, 1381, 1539, 2834, 1099, 4868)
  ) %>% 
  mutate(
    n_norm_phenotypes = n_associations / n_phenotypesAll,
    n_norm_genes      = n_genes / n_genesAll,
    n_norm_source     = n_sources / 3
  ) %>% 
  mutate(
    # 0–1 scaling każdej kolumny osobno
    n_norm_phenotypes_01 =
      (n_norm_phenotypes - min(n_norm_phenotypes, na.rm = TRUE)) /
      (max(n_norm_phenotypes, na.rm = TRUE) - min(n_norm_phenotypes, na.rm = TRUE)),
    
    n_norm_genes_01 =
      (n_norm_genes - min(n_norm_genes, na.rm = TRUE)) /
      (max(n_norm_genes, na.rm = TRUE) - min(n_norm_genes, na.rm = TRUE)),
    
    n_norm_source_01 =
      (n_norm_source - min(n_norm_source, na.rm = TRUE)) /
      (max(n_norm_source, na.rm = TRUE) - min(n_norm_source, na.rm = TRUE))
  ) %>% 
  mutate(
    score = (n_norm_phenotypes_01 + n_norm_genes_01 + n_norm_source_01) / 3
  ) %>%
  mutate(
    score_01 =
      (score - min(score, na.rm = TRUE)) /
      (max(score, na.rm = TRUE) - min(score, na.rm = TRUE))
  ) %>%
  # ---- jedna kolejność osi Y dla wszystkich wykresów: wg score ----
arrange(score) %>%  # low -> high (high będzie na górze)
  mutate(icdcategory = factor(icdcategory, levels = icdcategory)) %>% 
  mutate(sum_score = n_norm_genes_01 + n_norm_phenotypes_01 + n_norm_source_01)


# kolory / styl
bar_fill  <- "#552c17ff"   # earth-tone brown
bar_edge  <- "#4A3A28"
bar_width <- 0.7
bar_lwd   <- 1

title_pheno <- "phenotypes"
title_genes <- "genes"
title_src   <- "sources"
title_score <- "mean of component"

x_scale <- scale_x_continuous(
  limits = c(0, 1),
  breaks = c(0, 0.5, 1),
  position = "top"
)

# wspólny motyw: wszystko czarne
theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# ---- plots ----
p1 <- ggplot(df_plot, aes(x = n_norm_phenotypes, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  scale_x_continuous(
    position = "top"
  ) +
  # labs(x = "0–1", y = NULL, title = title_pheno) +
  theme_classic() +
  theme_black_text +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p2 <- ggplot(df_plot, aes(x = n_associations, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  scale_x_continuous(
    position = "top"
  ) +
  # labs(x = "0–1", y = NULL, title = title_genes) +
  theme_classic() +
  theme_black_text

p3 <- ggplot(df_plot, aes(x = n_genes, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  scale_x_continuous(
    position = "top"
  ) +
  # labs(x = "0–1", y = NULL, title = title_pheno) +
  theme_classic() +
  theme_black_text +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p4 <- ggplot(df_plot, aes(x = n_norm_genes, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  scale_x_continuous(
    position = "top"
  ) +
  # labs(x = "0–1", y = NULL, title = title_pheno) +
  theme_classic() +
  theme_black_text +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )


# ---- combine ----
fig <- (p1 | p2 | p3 | p4)

# ---- save ----
out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures"
out_file <- "summaryAssociations_grSystemic_test1_05.02.2026.svg"

svg(
  filename  = file.path(out_dir, out_file),
  width     = 12,
  height    = 6,
  pointsize = 10
)
print(fig)
dev.off()




# ##############################################################################
# ---- summary ----
gene_list <- c(gene_list_disgenet, gene_list_genebass, gene_list_GWASCatalog)

keep <- grepl("(^|_)(G3x|F[0-9]x)(_|$)", names(gene_list)) &
  !grepl("F[0-9]x_.*F[0-9]x", names(gene_list))


keep <-
  (
    grepl("(^|_)(G3x|F[0-9]x)(_|$)", names(gene_list)) &
      !grepl("F[0-9]x_.*F[0-9]x", names(gene_list))
  ) |
  grepl("neurodegenerative", names(gene_list), ignore.case = TRUE)

gene_list[grepl("G3x|F0x|F1x|F2x|F3x|F4x|F5x|F6x|F7x|F8x|F9x", names(gene_list))]

gene_list_filtered <- gene_list[keep]

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df

gene_list_filtered[grepl("neurodegenerative",  names(gene_list_filtered))] %>% names

gene_list_filtered[grepl("neurodegenerative",  names(gene_list_filtered))] %>% unname %>% unlist %>% unique %>% length()

gene_list_filtered[grepl("F3x",  names(gene_list_filtered))] %>% length()

GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>%
  filter(
    grepl("(^|_)(G3x|F[0-9]x)(_|$)", Var1),
    !grepl("F[0-9]x_.*F[0-9]x", Var1)
  )


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  filter(p_value < 0.05) %>% 
  filter(log2_odds_ratio > 0) %>% 
  filter(gene_overlap_count >= 3) %>% 
  filter(grepl("F3x", Var1)) %>% .$p_value %>% summary


GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>%
  filter(p_value < 0.05) %>% 
  filter(log2_odds_ratio > 0) %>% 
  filter(gene_overlap_count >= 3) %>% 
  filter(grepl("neurodegenerative", Var1)) %>% 
  .$overlap_genes %>% strsplit(",") %>% unlist %>% unique() %>% length()
