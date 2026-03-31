# summary database, summary performed on HGNC symbols

papers_data_preprocessing %>% 
  # .$hgnc_symbol %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% # clean treatment 
  mutate(treatment = ifelse(treatment == "corticoterone", "corticosterone", treatment)) %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  .$hgnc_symbol %>% unique %>% length()

# 1. count genes
papers_data_preprocessing %>% 
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  arrange(desc(Freq)) %>%
  .$Freq %>% 
  table() %>% 
  as.data.frame() %>% 
  set_colnames(c("occurence_gene", "number_genes")) %>% 
  arrange(desc(number_genes)) %>% 
  mutate(cum_sum = cumsum(number_genes)) %>% 
  mutate(cum_perc = cum_sum/13251) %>%
  mutate(perc_point = number_genes/13251) 


# Przetwarzanie danych z uproszczeniem górnych 5%
processed_data <- papers_data_preprocessing %>%
  .$hgnc_symbol %>%
  table() %>%
  as.data.frame() %>%
  .$Freq %>%
  table() %>%
  as.data.frame() %>%
  setNames(c("occurence_gene", "number_genes")) %>%
  arrange(desc(as.numeric(occurence_gene))) %>%
  mutate(cum_sum = cumsum(number_genes)) %>%
  mutate(cum_perc = cum_sum / sum(number_genes)) %>%
  mutate(category = ifelse(cum_perc < 0.055, "Top 5%", as.character(occurence_gene))) %>%
  group_by(category) %>%
  summarise(number_genes = sum(number_genes)) %>%
  arrange(desc(as.numeric(category))) 

# Generowanie histogramu
ggplot(processed_data, aes(x = factor(category, levels = c(1:14, "Top 5%")), y = number_genes)) +
  geom_bar(stat = "identity", fill = "blue", alpha = 0.7) +
  labs(
    title = "Histogram of Gene Occurrences (Top 5% Simplified)",
    x = "Occurrence of Gene",
    y = "Number of Genes"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12)
  )


# 2. regulation: up, down
papers_data_preprocessing %>%  .$regulation %>% table 
  
papers_data_preprocessing %>%
  filter(regulation %in% c("down", "up")) %>%
  select(hgnc_symbol, regulation) %>%
  distinct() %>%
  group_by(hgnc_symbol) %>%
  dplyr::summarize(regulations = paste(sort(unique(regulation)), collapse = ", ")) %>%
  mutate(category = case_when(
    regulations == "down" ~ "only_down",
    regulations == "up" ~ "only_up",
    regulations == "down, up" ~ "both"
  )) %>%
  count(category)

venn_data <- papers_data_preprocessing %>%
  filter(regulation %in% c("down", "up")) %>%
  select(hgnc_symbol, regulation) %>%
  distinct() %>%
  group_by(regulation) %>%
  summarise(genes = list(hgnc_symbol)) %>%
  deframe()

# Sprawdzenie danych
print(venn_data)

venn_data <- lapply(venn_data, function(x) x[!is.na(x)])

# Diagram Venn
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("Up", "Down"),
  filename = NULL, # Wyświetlenie w RStudio
  fill = c("white", "white"), # Odcienie szarości
  alpha = 0.5,
  cex = 1.5, # Rozmiar tekstu
  cat.cex = 1.5, # Rozmiar nazw kategorii
  cat.col = "black", # Kolor nazw kategorii
  margin = 0.1
)

grid::grid.draw(venn.plot)

# 3. treatment (GC)
papers_data_preprocessing %>% 
  filter(treatment %in%  c("dexamethasone", "corticosterone", "cortisol", "prednisone", "prednisolone", "hydrocortisone", "budesonide")) %>% 
  select(c(hgnc_symbol, treatment)) %>% unique %>% 
  mutate(treatment = ifelse(treatment == "prednisone", "prednisone/prednisolone", treatment)) %>% 
  mutate(treatment = ifelse(treatment == "prednisolone", "prednisone/prednisolone", treatment)) %>% 
  mutate(treatment = ifelse(treatment == "cortisol", "cortisol/corticosterone", treatment)) %>% 
  mutate(treatment = ifelse(treatment == "corticosterone", "cortisol/corticosterone", treatment)) %>%
  unique 

  

# Podział danych na grupy
grouped_data <- papers_data_preprocessing %>%
  filter(treatment %in% c("dexamethasone", "corticosterone", "cortisol", "prednisone", "prednisolone", "hydrocortisone", "budesonide")) %>%
  select(hgnc_symbol, treatment) %>%
  unique() %>%
  mutate(treatment = ifelse(treatment == "prednisone", "prednisone/prednisolone", treatment)) %>%
  mutate(treatment = ifelse(treatment == "prednisolone", "prednisone/prednisolone", treatment)) %>%
  mutate(treatment = ifelse(treatment == "cortisol", "cortisol/corticosterone", treatment)) %>%
  mutate(treatment = ifelse(treatment == "hydrocortisone", "cortisol/corticosterone", treatment)) %>%
  mutate(treatment = ifelse(treatment == "corticosterone", "cortisol/corticosterone", treatment)) %>%
  unique()

grouped_list <- split(grouped_data$hgnc_symbol, grouped_data$treatment)

grouped_list <- lapply(grouped_list, function(x) x[!is.na(x)])

# Rysowanie wykresu Venn
venn.plot <- venn.diagram(
  x = grouped_list,
  # category.names = c("A", "B", "C", "D"),
  filename = NULL, # NULL, aby wyświetlić w RStudio bez zapisu do pliku
  fill = c("red", "green", "blue", "yellow"),
  alpha = 0.5,
  cex = 1,
  cat.cex = 1,
  cat.col = c("red", "green", "blue", "yellow"),
  margin = 0.1
)

lapply(grouped_list, length)


grid::grid.draw(venn.plot)


grouped_data %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_gc = map_int(data, ~ length(.x$treatment))) %>% 
  arrange(desc(n_gc)) %>% 
  filter(n_gc == 5) %>% 
  select(-data) %>% 
  as.data.frame()

papers_data_preprocessing %>% 
  filter(treatment_type %in% c("acute", "chronic", "stress_induction")) %>% 
  select(hgnc_symbol, treatment_type) %>% 
  unique 


# Przekształcenie danych na listy
venn_data <- papers_data_preprocessing %>%
  filter(treatment_type %in% c("acute", "chronic", "stress_induction")) %>%
  select(hgnc_symbol, treatment_type) %>%
  unique() %>%
  group_by(treatment_type) %>%
  summarise(genes = list(hgnc_symbol)) %>%
  deframe()

# Wygląd listy dla VennDiagram
print(venn_data)

venn_data <- lapply(venn_data, function(x) x[!is.na(x)])

# Diagram Venn w czerni i bieli
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("Acute", "Chronic", "Stress Induction"),
  filename = NULL, # Bez zapisu do pliku
  fill = c("gray80", "gray80", "gray80"), # Odcienie szarości
  alpha = 0.5, # Przezroczystość
  cex = 1, # Rozmiar tekstu
  cat.cex = 1, # Rozmiar nazw kategorii
  cat.col = "black", # Kolor nazw kategorii
  margin = 0.1
)

grid::grid.draw(venn.plot)

papers_data_preprocessing %>% 
  select(hgnc_symbol, simple_tissue) %>% 
  unique %>% 
  .$simple_tissue %>% table


