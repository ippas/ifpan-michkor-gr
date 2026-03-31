# boxplot
filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(hgnc_occurence > 35) %>%
  # filter(hgnc_symbol %in% c("PLAU", "ENC1", "ARL4C", "RGS3", "HBEGF", "GEM", 
  #                           "PTGS2", "PTPRE", "IL11", "IER3", "IL6", "PLK2", 
  #                           "IL1B", "CCL2", "JUN", "LIF", "G0S2", "ADAMTS1", 
  #                           "SOX4", "DUSP5", "PER2")) %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  mutate(log2ratio = as.numeric(log2ratio),
         hgnc_symbol = fct_reorder(hgnc_symbol, hgnc_occurence, .desc = TRUE)) %>%  # Sortowanie wg hgnc_occurence malejąco
  ggplot(aes(x = hgnc_symbol, y = log2ratio)) +
  geom_boxplot(outlier.shape = NA, size = 1) +  # Pogrubiony boxplot
  geom_jitter(width = 0.2, alpha = 1, color = "black") +  # Punkty
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "firebrick", size = 1) +  # Linia na y = 0.5 (czerwona)
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "blue4", size = 1) +  # Linia na y = -0.5 (niebieska)
  theme_minimal() +
  labs(title = "Ekspresja genów (log2ratio)",
       x = "Gen",
       y = "log2ratio") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 





################################################################################


plot_log2ratio_gr_database <- function(data, hgnc_symbol) {
  data %>%
    filter(regulation %in% c("down", "up")) %>% 
    filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
    filter(hgnc_symbol == !!hgnc_symbol) %>%
    mutate(log2ratio = as.numeric(log2ratio)) %>% 
    mutate(label = paste(label, treatment, dose, time, treatment_type, environment, comparison, sep = "_")) %>%
    select(hgnc_symbol, regulation, label, log2ratio) %>%
    ggplot(aes(x = log2ratio, y = reorder(label, log2ratio), color = regulation)) +
    geom_point(size = 4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = paste("Wykres log2ratio dla", hgnc_symbol, "w różnych eksperymentach"),
      x = "log2ratio",
      y = "Label",
      color = "Regulation"
    ) +
    scale_color_manual(values = c("down" = "blue4", "up" = "firebrick")) +
    theme_minimal()
}


plot_log2ratio_gr_database <- function(data, hgnc_symbols, nrow = 1, point_size = 4, order_by_gene = NULL) {
  
  # Filtrowanie i przygotowanie danych
  filtered_data <- data %>%
    filter(regulation %in% c("down", "up")) %>% 
    filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
    filter(!is.na(log2ratio)) %>% 
    filter(hgnc_symbol %in% hgnc_symbols) %>%
    mutate(log2ratio = as.numeric(log2ratio)) %>% 
    mutate(full_label = paste(label, treatment, dose, time, treatment_type, environment, comparison, sep = "_")) 
  
  # Ustalanie kolejności etykiet dla osi y, jeśli podano gen do sortowania
  if (!is.null(order_by_gene) && order_by_gene %in% filtered_data$hgnc_symbol) {
    ref_order <- filtered_data %>%
      filter(hgnc_symbol == order_by_gene) %>%
      arrange(desc(log2ratio), is.na(log2ratio)) %>%
      pull(full_label) %>%
      unique()
    
    all_labels   <- unique(filtered_data$full_label)
    other_labels <- sort(setdiff(all_labels, ref_order))
    new_levels   <- c(ref_order, other_labels)
    
    filtered_data <- filtered_data %>%
      mutate(full_label = factor(full_label, levels = new_levels))
  } else {
    new_levels <- unique(filtered_data$full_label)
  }
  
  # Tworzenie wykresu z legendą poniżej i odwróconą osią y
  ggplot(filtered_data, aes(x = log2ratio, y = full_label, color = regulation)) +
    geom_point(size = point_size) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
    labs(
      title = "Wykres log2ratio dla wybranych genów w różnych eksperymentach",
      x = "log2ratio",
      y = "Etykieta eksperymentu",
      color = "Regulacja"
    ) +
    scale_color_manual(values = c("down" = "blue4", "up" = "firebrick")) +
    facet_wrap(~ hgnc_symbol, nrow = nrow) +
    scale_y_discrete(limits = rev(new_levels)) +
    theme_minimal() +
    theme(legend.position = "bottom")
}


plot_log2ratio_gr_database(
  data = filt_1_12_gr_database, 
  hgnc_symbol = c("FKBP5", "TSC22D3", "KLF9", "PER1", "PDK4"), 
  point_size = 2,
  order_by_gene = "FKBP5"
)

svg("results/figures/gene-database-summary/filtered-1-12-gene-database-summary/expression-top5-genes.svg", 
    height = 14,
    width = 20)
plot_log2ratio_gr_database(
  data = filt_1_12_gr_database, 
  hgnc_symbol = c("FKBP5", "TSC22D3", "KLF9", "PER1", "PDK4"), 
  point_size = 2,
  order_by_gene = "FKBP5"
)
dev.off()


plot_log2ratio_gr_database(
  data = filt_1_12_gr_database, 
  hgnc_symbol = c("PTGS2", "PER2", "SOX4", "DUSP5", "IL6"), 
  point_size = 2,
  order_by_gene = "DUSP5" 
)


################################################################################
# text plot

# Parametry
order_by_gene <- "FKBP5"
hgnc_symbols <- c("FKBP5", "TSC22D3", "KLF9", "PER1", "PDK4")

# Krok 1: Filtrowanie i przygotowanie danych, ustalanie kolejności etykiet
filtered_data <- filt_1_12_gr_database %>%
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(!is.na(log2ratio)) %>% 
  filter(hgnc_symbol %in% hgnc_symbols) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  mutate(full_label = paste(label, treatment, dose, time, treatment_type, environment, comparison, sep = "_"))

if (!is.null(order_by_gene) && order_by_gene %in% filtered_data$hgnc_symbol) {
  ref_order <- filtered_data %>%
    filter(hgnc_symbol == order_by_gene) %>%
    arrange(desc(log2ratio), is.na(log2ratio)) %>%
    pull(full_label) %>%
    unique()
  
  all_labels   <- unique(filtered_data$full_label)
  other_labels <- sort(setdiff(all_labels, ref_order))
  new_levels   <- c(ref_order, other_labels)
  
  filtered_data <- filtered_data %>%
    mutate(full_label = factor(full_label, levels = new_levels))
} else {
  new_levels <- unique(filtered_data$full_label)
}

# Wyciągamy kolejność etykiet dla 'label' – przyjmujemy, że dla każdego label istnieje odpowiadający full_label
labels_order <- filtered_data %>%
  distinct(label, full_label) %>%
  arrange(match(full_label, new_levels)) %>%
  pull(label)

# Wyznaczamy etykiety do filtrowania – tylko te, które dotyczą naszych genów
labels_to_filt <- filtered_data %>%
  filter(hgnc_symbol %in% hgnc_symbols) %>% 
  pull(label) %>% 
  unique()

# Krok 2: Przetwarzanie danych z filt_1_12_gr_database przy użyciu nest() i map()
final_data <- filt_1_12_gr_database %>%
  filter(regulation %in% c("down", "up"),
         !log2ratio %in% c("", "NA", "Inf", "-Inf"),
         !is.na(log2ratio)) %>%
  group_by(label) %>%
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  unique %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  nest() %>%
  mutate(
    top3_up = map(data, ~ .x %>%
                    slice_max(order_by = log2ratio, n = 3, with_ties = FALSE) %>%
                    select(hgnc_symbol, log2ratio)),
    top3_down = map(data, ~ .x %>%
                      slice_min(order_by = log2ratio, n = 3, with_ties = FALSE) %>%
                      select(hgnc_symbol, log2ratio)),
    top3_up_genes = map_chr(data, ~ .x %>%
                              slice_max(order_by = log2ratio, n = 3, with_ties = FALSE) %>%
                              pull(hgnc_symbol) %>%
                              paste(collapse = ", ")),
    top3_down_genes = map_chr(data, ~ .x %>%
                                slice_min(order_by = log2ratio, n = 3, with_ties = FALSE) %>%
                                pull(hgnc_symbol) %>%
                                paste(collapse = ", "))
  ) %>% unique
  select(label, top3_up_genes, top3_down_genes) %>%
  filter(label %in% labels_to_filt) %>%
  mutate(label = factor(label, levels = labels_order)) %>%
  arrange(label)

final_data


# Krok 1. Rozdzielenie top3_up_genes i top3_down_genes na osobne kolumny
final_data_sep <- final_data %>%
  separate(top3_up_genes, into = c("up1", "up2", "up3"), sep = ", ") %>%
  separate(top3_down_genes, into = c("down1", "down2", "down3"), sep = ", ")

# Krok 2. Przekształcenie do formatu long
final_data_long <- final_data_sep %>%
  pivot_longer(
    cols = c(up1, up2, up3, down1, down2, down3),
    names_to = "dir_rank",
    values_to = "gene"
  ) %>%
  # Wydzielamy informację o kierunku (Up/Down) oraz pozycji (1, 2, 3)
  mutate(
    direction = if_else(grepl("^up", dir_rank), "Up", "Down"),
    rank = sub("^(up|down)", "", dir_rank),
    xaxis = paste(direction, rank)  # tworzymy etykietę np. "Up 1", "Down 2", ...
  )

# Ustalamy kolejność dla osi x (6 kolumn)
x_levels <- c("Up 1", "Up 2", "Up 3", "Down 1", "Down 2", "Down 3")
final_data_long <- final_data_long %>%
  mutate(xaxis = factor(xaxis, levels = x_levels))

# Krok 3. Dodanie kolumny z etykietą dla genów z wyrażeniem plotmath:
bold_genes <- c("FKBP5", "TSC22D3", "PER", "PDK4", "KLF9")
final_data_long <- final_data_long %>%
  mutate(gene_label = if_else(gene %in% bold_genes,
                              paste0("bold('", gene, "')"),
                              paste0("'", gene, "'"))
  )


svg("results/figures/gene-database-summary/filtered-1-12-gene-database-summary/expression-top3-text.svg", 
    height = 14,
    width = 20)
# Krok 4. Tworzenie wykresu tekstowego z białym tłem, wyrównaniem do lewej, parse=TRUE aby interpretować plotmath
ggplot(final_data_long, aes(x = xaxis, y = label, label = gene_label)) +
  geom_text(hjust = 0, size = 4, parse = TRUE) +
  labs(x = "", y = "Label", title = "Top 3 Genes Up i Down (rozdzielone na 6 kolumn)") +
  theme_classic() +
  scale_y_discrete(limits = rev(levels(final_data$label)))

dev.off()
