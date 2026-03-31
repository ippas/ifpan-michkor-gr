
# ---- input ----
icd_df <- data.frame(
  icdcategory = c("F0", "F1", "F2", "F3", "F4", "F5", "F6", "F7", "F8", "F9"),
  n_associations = c(11, 7, 4, 17, 1, 0, 0, 0, 3, 0),
  n_genes = c(84, 44, 39, 58, 10, 0, 0, 0, 15, 0),
  n_sources = c(3, 2, 1, 2, 1, 0, 0, 0, 1, 0),
  stringsAsFactors = FALSE
)

# ---- prepare data (tylko do rankingu / score) ----
df_plot <- icd_df %>%
  mutate(
    n_phenotypesAll = c(36, 44, 27, 70, 27, 12, 4, 11, 17, 8),
    n_genesAll      = c(5743, 6157, 5555, 7558, 4973, 3295, 1381, 1539, 2834, 1099)
  ) %>%
  mutate(
    n_norm_phenotypes = n_associations / n_phenotypesAll,
    n_norm_genes      = n_genes / n_genesAll,
    n_norm_source     = n_sources / 3
  ) %>%
  mutate(
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
  arrange(score) %>%  # wspólna kolejność osi Y
  mutate(icdcategory = factor(icdcategory, levels = icdcategory))

# ---- style ----
bar_fill  <- "#552c17ff"
bar_edge  <- "#4A3A28"
bar_width <- 0.7
bar_lwd   <- 1

title_pheno <- "phenotypes (raw)"
title_genes <- "genes (raw)"
title_src   <- "sources (raw)"
title_score <- "mean component"

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# ---- RAW x-scales ----
x_raw_pheno  <- scale_x_continuous(position = "top", expand = expansion(mult = c(0, 0.05)))
x_raw_genes  <- scale_x_continuous(position = "top", expand = expansion(mult = c(0, 0.05)))
x_raw_source <- scale_x_continuous(position = "top", breaks = 0:3, limits = c(0, 3))

# ---- plots ----
p1 <- ggplot(df_plot, aes(x = n_genes, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  x_raw_genes +
  labs(x = "count", y = NULL, title = title_genes) +
  theme_classic() +
  theme_black_text

p2 <- ggplot(df_plot, aes(x = n_associations, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  x_raw_pheno +
  labs(x = "count", y = NULL, title = title_pheno) +
  theme_classic() +
  theme_black_text +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p3 <- ggplot(df_plot, aes(x = n_sources, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 0.7) +
  x_raw_source +
  labs(x = "count", y = NULL, title = title_src) +
  theme_classic() +
  theme_black_text +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank()
  )

p4 <- ggplot(df_plot, aes(x = score_01, y = icdcategory)) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge, size = bar_lwd, alpha = 1) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1), position = "top") +
  labs(x = "0–1", y = NULL, title = title_score) +
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
out_file <- "summaryAssociations_grSystemic_raw_v1_04.02.2026.svg"

svg(
  filename  = file.path(out_dir, out_file),
  width     = 12,
  height    = 6,
  pointsize = 10
)
print(fig)
dev.off()




# ##############################################################################
# ===============================
# ścieżka i nazwa pliku
# ===============================
out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures"
out_file <- "summaryAssociations_nNormGenes_v1_04.02.2026.svg"

dev.off()
# ===============================
# otwórz urządzenie SVG
# ===============================
svg(
  filename  = file.path(out_dir, out_file),
  width     = 4,
  height    = 6,
  pointsize = 10
)

# ===============================
# RYSUJ WYKRES
# ===============================
p_genes <- ggplot(df_plot, aes(x = n_norm_genes, y = icdcategory)) + 
  geom_col(width = bar_width, fill = bar_fill, color = "black", size = bar_lwd, alpha = 0.7) + 
  scale_x_continuous(position = "top" ) +
    # labs(x = "0–1", y = NULL, title = title_genes) + 
  theme_classic() + 
  theme_black_text

print(p_genes)

# ===============================
# zamknij urządzenie
# ===============================
dev.off()

