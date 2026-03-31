library(dplyr)
library(tidyr)
library(ggplot2)

# =========================================
# 1. wybór traitów powiązanych z ieu-a-1187
# =========================================
ldsc_results_annotated %>% 
  filter(p1_subcategory == "Psychiatric / neurological") %>% 
  select(p1_id, p1_trait) %>% unique



selected_traits <- ldsc_results_annotated %>%
  filter(p1_id == "ieu-a-806", summary_rg_p < 0.05) %>%
  pull(p2_id) %>%
  unique()


selected_traits <- unique(c("ieu-a-806", selected_traits))

# =========================================
# 2. pobranie wszystkich wyników dla wybranych traitów
# =========================================
pairwise_rg <- ldsc_results_annotated %>%
  filter(p1_id %in% selected_traits, p2_id %in% selected_traits) %>%
  filter(!is.na(summary_rg)) %>%
  select(p1_id, p2_id, summary_rg, summary_rg_p)

# =========================================
# 3. sprowadzenie A-B i B-A do jednej wspólnej pary
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
# 4. rozpisanie z powrotem na obie strony
# =========================================
pairwise_sym <- bind_rows(
  pairwise_unique %>%
    transmute(
      p1_id = trait_min,
      p2_id = trait_max,
      summary_rg,
      summary_rg_p
    ),
  pairwise_unique %>%
    transmute(
      p1_id = trait_max,
      p2_id = trait_min,
      summary_rg,
      summary_rg_p
    )
)

# =========================================
# 5. pełna siatka wszystkich kombinacji
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
    summary_rg = ifelse(p1_id == p2_id, 1, summary_rg),
    summary_rg_p = ifelse(p1_id == p2_id, NA, summary_rg_p)
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
# 7. macierz do klastrowania
#    brakujące wartości tylko technicznie na 0
# =========================================
rg_matrix <- pairwise_rg_full %>%
  select(p1_id, p2_id, summary_rg) %>%
  pivot_wider(names_from = p2_id, values_from = summary_rg) %>%
  as.data.frame()

rownames(rg_matrix) <- rg_matrix$p1_id
rg_matrix$p1_id <- NULL
rg_matrix <- as.matrix(rg_matrix)

rg_matrix_for_clustering <- rg_matrix
rg_matrix_for_clustering[is.na(rg_matrix_for_clustering)] <- 0

# jedno wspólne klastrowanie
hc <- hclust(dist(rg_matrix_for_clustering))
trait_order <- hc$labels[hc$order]

# =========================================
# 8. dane do wykresu
# =========================================
pairwise_rg_plot <- pairwise_rg_full %>%
  mutate(
    p1_id = factor(p1_id, levels = trait_order),
    p2_id = factor(p2_id, levels = trait_order)
  )

# =========================================
# 9. heatmapa
# =========================================
ggplot(pairwise_rg_plot, aes(x = p1_id, y = p2_id, fill = summary_rg)) +
  geom_tile(color = "white") +
  scale_x_discrete(labels = label_vector) +
  scale_y_discrete(labels = label_vector) +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    na.value = "grey85"
  ) +
  coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank()
  ) +
  labs(
    fill = "rg",
    title = "Genetic correlation heatmap"
  )

