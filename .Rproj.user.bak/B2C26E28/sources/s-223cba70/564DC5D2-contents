plot_df <- chea_df_raw %>%
  mutate(
    TF = stringr::str_split(Term, " ", simplify = TRUE)[, 1],
    logP = -log10(Adjusted.P.value)
  ) %>%
  group_by(cluster) %>%
  arrange(Adjusted.P.value) %>%
  slice_head(n = 10) %>%
  mutate(TF_rank = paste0("top", row_number())) %>%
  ungroup() %>% 
  mutate(
    nuclear_receptors = TF %in% nuclear_receptors$hgnc_symbol,
    logP = ifelse(logP > 10, 10, logP)
  ) %>% 
  filter(Adjusted.P.value < 0.1)

plot_df %>% head

plot_df %>% .$P.value %>% max

TF_mean_sal_df <- read.delim("data/review-dextis/TF-mean-sal.tsv",
           sep = "\t")

TF_mean_sal_df %>% head
  
left_join(
  plot_df,
  TF_mean_sal_df,
  by = c("TF" = "original", "cluster" = "cluster")
) %>% 
  # filter(gene_name == "NR3C1") %>% as.data.frame()
  mutate(threshold = ifelse(mean > 8, T, F)) %>% 
  group_by(cluster, gene_name) %>% 
  nest %>% 
  mutate(sum_threshold = map_dbl(data, ~ sum(.x$threshold)))  %>% 
  unnest(sum_threshold) %>% 
  # filter(sum_threshold == 0) %>% 
  mutate(tissue_agree = ifelse(sum_threshold == 0, FALSE, TRUE)) %>% 
  unnest(data) -> plot_df2



# Uzupełnienie siatki, by zachować puste klastry
all_clusters <- unique(chea_df_raw$cluster)
complete_df <- expand.grid(
  cluster = all_clusters,
  TF_rank = paste0("top", 1:10)
)

plot_df_complete <- complete_df %>%
  left_join(plot_df2, by = c("cluster", "TF_rank"))

# Finalny wykres
tf_plot  <- ggplot(plot_df_complete, aes(x = cluster, y = TF_rank, fill = logP)) +
  geom_tile(width = 0.95, height = 0.95) +
  # geom_point(
  #   data = subset(plot_df_complete, tissue_agree == TRUE),
  #   aes(x = cluster, y = TF_rank),
  #   shape = 15, size = 6, color = "#69995D",
  #   position = position_nudge(y = -0.2, x = 0.12),
  #   na.rm = TRUE
  # ) +
  geom_point(
    data = subset(plot_df_complete, tissue_agree == TRUE),
    aes(x = cluster, y = TF_rank),
    shape = 16, size = 6, color = "#1e4620",
    position = position_nudge(x = 0.3, y = 0.3),
    na.rm = TRUE
  )+
  # geom_point(
  #   data = subset(plot_df_complete, tissue_agree5 == TRUE),
  #   aes(x = cluster, y = TF_rank),
  #   shape = 16, size = 6, color = "#449e48",
  #   position = position_nudge(x = 0.1, y = 0.3),
  #   na.rm = TRUE
  # )+
  geom_tile(
    data = subset(plot_df_complete, nuclear_receptors == TRUE),
    aes(x = cluster, y = TF_rank),
    fill = NA,                      # brak wypełnienia
    color = "#ab87c9",              # fioletowa ramka
    size = 1,                # grubość ramki
    width = 0.95, height = 0.95
  ) +
  geom_tile(
    data = subset(plot_df_complete, TF == "NR3C1"),
    aes(x = cluster, y = TF_rank),
    fill = NA,                      # brak wypełnienia
    color = "#451F55",              # fioletowa ramka
    size = 1,                # grubość ramki
    width = 0.95, height = 0.95
  ) +
  
  geom_text(aes(label = TF), size = 5, color = "black",
            position = position_nudge(y = 0, x = 0),
            na.rm = TRUE
  ) +
  scale_y_discrete(limits = paste0("top", 10:1)) +
  scale_x_discrete(position = "top", limits = sort(all_clusters)) +
  
  # 🔥 KLUCZOWA ZMIANA TUTAJ 🔥 #
  scale_fill_gradient(
    low = "#ffffff", high = "#b2182b",
    limits = c(0, 10), # Skala od 0 do 10
    name = expression(-log[10](italic(p))),
    na.value = "white"
  ) +
  
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 1.5),
    axis.line = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0, hjust = 0, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank(),
    legend.position = "bottom"
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))


dev.off()

# svg("results/figures/review-dextis/enrichr_GR_signatures_tissue_agree.svg", width = 20, height = 12.5)
tf_plot
dev.off()
F###

get_legend <- function(myplot) {
  library(ggplot2)
  tmp <- ggplotGrob(myplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  if (length(leg) > 0) {
    return(tmp$grobs[[leg]])
  } else {
    return(NULL)
  }
}

# Wyciągnij legendę z wykresu
legend_only <- get_legend(tf_plot)


svg("results/figures/review-dextis/legendR.svg", width = 5, height = 1.5)
legend_only
dev.off()
