library(dplyr)
library(tidyr)
library(pheatmap)

# =========================================
# 1. wybór traitów powiązanych z ieu-a-1187
# =========================================
selected_traits <- ldsc_results_annotated %>%
  filter(p1_id == "ieu-a-1187", summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()

selected_traits <- unique(c("ieu-a-1187", selected_traits))

# =========================================
# 2. pobranie wszystkich wyników dla wybranych traitów
# =========================================
pairwise_rg <- ldsc_results_annotated %>%
  filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
  filter(!is.na(summary_rg)) %>%
  select(p1_id, p2_id, summary_rg, summary_rg_p)

# =========================================
# 3. sprowadzenie A-B i B-A do jednej pary
# =========================================
pairwise_unique <- pairwise_rg %>%
  rowwise() %>%
  mutate(
    trait_min = sort(c(p1_id, p2_id))[1],
    trait_max = sort(c(p1_id, p2_id))[2]
  ) %>%
  ungroup() %>%
  group_by(trait_min, trait_max) %>%
  summarise(
    summary_rg = first(na.omit(summary_rg)),
    summary_rg_p = first(na.omit(summary_rg_p)),
    .groups = "drop"
  )

# =========================================
# 4. rozpisanie na obie strony
# =========================================
pairwise_sym <- bind_rows(
  pairwise_unique %>%
    transmute(
      p1_id = trait_min,
      p2_id = trait_max,
      summary_rg = summary_rg
    ),
  pairwise_unique %>%
    transmute(
      p1_id = trait_max,
      p2_id = trait_min,
      summary_rg = summary_rg
    )
)

# =========================================
# 5. pełna siatka + przekątna = 1
# =========================================
full_grid <- expand.grid(
  p1_id = selected_traits,
  p2_id = selected_traits,
  stringsAsFactors = FALSE
) %>%
  as_tibble()

pairwise_rg_full <- full_grid %>%
  left_join(pairwise_sym, by = c("p1_id", "p2_id")) %>%
  mutate(
    summary_rg = ifelse(p1_id == p2_id, 1, summary_rg)
  )

# =========================================
# 6. etykiety: id | trait
# =========================================
trait_labels_df <- metadata_ieu_EURsampleSize10000 %>%
  select(id, trait) %>%
  distinct() %>%
  mutate(label = paste0(id, " | ", trait))

label_vector <- trait_labels_df$label
names(label_vector) <- trait_labels_df$id

missing_ids <- setdiff(selected_traits, names(label_vector))
label_vector[missing_ids] <- missing_ids

# =========================================
# 7. macierz rg
# =========================================
rg_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg) %>%
  as.data.frame()

rownames(rg_matrix) <- rg_matrix$p1_id
rg_matrix$p1_id <- NULL
rg_matrix <- as.matrix(rg_matrix)

# nazwy osi
rownames(rg_matrix) <- label_vector[rownames(rg_matrix)]
colnames(rg_matrix) <- label_vector[colnames(rg_matrix)]

# =========================================
# 8. macierz do klastrowania
# brakujące wartości tylko technicznie do dendrogramu
# =========================================
rg_matrix_for_clustering <- rg_matrix
rg_matrix_for_clustering[is.na(rg_matrix_for_clustering)] <- 0

hc <- hclust(dist(rg_matrix_for_clustering))

# =========================================
# 9. heatmapa z dendrogramami
# =========================================
pheatmap(
  mat = rg_matrix,
  cluster_rows = hc,
  cluster_cols = hc,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-1, 1, length.out = 101),
  na_col = "grey85",
  border_color = "white",
  main = "Genetic correlation heatmap",
  angle_col = "90"
)

