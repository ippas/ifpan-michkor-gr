
# Dane + wykres + zapis do pliku SVG
svg("results_v2/summary_geneLists_brain/figures/piechart_summaryGeneListFilterN10ShortTime_perList.svg", 
    width = 6, height = 6)

tibble(
  freq = 1:15,
  n_genes = c(3503, 1605, 771, 343, 161, 73, 47, 21, 8, 15, 9, 9, 6, 3, 1)
) %>%
  mutate(freq_group = ifelse(freq >= 10, "≥10", as.character(freq))) %>%
  group_by(freq_group) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  arrange(desc(as.numeric(gsub("≥", "", freq_group)))) %>%
  { 
    pie(
      rev(.$n_genes),
      labels = paste0(rev(.$freq_group), " (", rev(.$n_genes), ")"),
      col = rep("grey90", nrow(.)),  # jasno szary kolor
      border = "black",              # czarna obwódka
      init.angle = 90,
      clockwise = FALSE,
      cex = 0.8,
      radius = 1
    )
  }

dev.off()  # ✅ zamyka plik SVG i zapisuje



# Tworzenie i zapis piechartu
svg("results_v2/summary_geneLists_brain/figures/piechart_summaryGeneListAll_perList.svg",
    width = 6, height = 6)

gene_summary_list$all$summary_all$freq_genes_per_list %>%
  mutate(freq_group = ifelse(freq >= 10, "≥10", as.character(freq))) %>%
  group_by(freq_group) %>%
  summarise(n_genes = sum(n_genes), .groups = "drop") %>%
  arrange(desc(as.numeric(gsub("≥", "", freq_group)))) %>%
  {
    pie(
      rev(.$n_genes),
      labels = paste0(rev(.$freq_group), " (", rev(.$n_genes), ")"),
      col = rep("grey90", nrow(.)),  # jasno szary kolor
      border = "black",              # czarna obwódka
      init.angle = 90,
      clockwise = FALSE,
      cex = 0.8,
      radius = 1
    )
  }

dev.off()  # ✅ zapisuje i zamyka plik SVG
