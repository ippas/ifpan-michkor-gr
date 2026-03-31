papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% 
  filter(regulation %in% c("down", "up")) %>% 
  select(c(hgnc_symbol, regulation, log2ratio)) %>% 
  unique %>% 
  filter(log2ratio != "NA") %>%
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  # filter(!log2ratio %in% c("", "NA")) %>% 
  select(c(hgnc_symbol, regulation)) %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~nrow(.x))) %>% 
  unnest(n_regulation) -> regulation_genes_df 

regulation_genes_df %>% 
  filter(n_regulation == 2) %>% 
  .$hgnc_symbol -> bidirectionally_regulated_genes

regulation_genes_df %>% 
  filter(n_regulation == 1) %>% 
  unnest(data) %>% 
  filter(regulation == "down") %>% 
  .$hgnc_symbol -> down_genes

regulation_genes_df %>% 
  filter(n_regulation == 1) %>% 
  unnest(data) %>% 
  filter(regulation == "up") %>% 
  .$hgnc_symbol -> up_genes

bidirectionally_regulated_genes %>% length()
up_genes %>% length()
down_genes %>% length()


papers_data_preprocessing %>% 
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>% 
  select(regulation, log2ratio) %>%  
  filter(regulation %in% c("up", "down")) %>% 
  .$regulation %>%
  table

papers_data_preprocessing %>% 
  filter(log2ratio != "NA") %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%   
  filter(regulation %in% c("down", "up")) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  mutate(group = ifelse(hgnc_symbol %in% up_genes, "up_only", regulation)) %>% 
  mutate(group = ifelse(hgnc_symbol %in% down_genes, "down_only", group)) %>% 
  mutate(group = factor(group, levels = c("down_only", "down", "up", "up_only"))) -> processed_gene_regulation

processed_gene_regulation %>%  
  select(group, log2ratio) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  ggplot(aes(x = group, y = log2ratio)) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Regulacja genów", y = "Log2(FC)")

processed_gene_regulation %>% 
  select(group, log2ratio) %>% 
  .$group %>% table
  
processed_gene_regulation %>%
  group_by(group) %>%
  dplyr::summarize(median_log2ratio = median(log2ratio, na.rm = TRUE))

processed_gene_regulation %>% 
  select(group, hgnc_symbol) %>%
  unique() %>% 
  .$group %>% table


processed_gene_regulation %>% 
  select(group, log2ratio) %>% 
  filter(group %in% c("up_only", "up")) %>%
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  # mutate(log2ratio = log10(as.numeric(log2ratio))) %>% 
  wilcox.test(log2ratio ~ group, data = .)

processed_gene_regulation %>% 
  select(group, log2ratio) %>% 
  filter(group %in% c("down_only", "down")) %>%
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  # mutate(log2ratio = log10(as.numeric(log2ratio))) %>% 
  wilcox.test(log2ratio ~ group, data = .)

processed_gene_regulation %>% 
  select(group, log2ratio) %>% 
  filter(group %in% c("up_only", "up")) %>%
  mutate(log2ratio = as.numeric(log2ratio)) %>%
  group_by(group) %>%
  summarise(values = list(log2ratio)) %>%
  summarise(ks_test_result = list(ks.test(values[[1]], values[[2]]))) %>%
  pull(ks_test_result) %>%
  print()


papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  dim


papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%   
  filter(regulation %in% c("down", "up")) %>% 
  select(hgnc_symbol, regulation) %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(n_regulation = map(data, ~nrow(.x))) %>% 
  unnest(n_regulation) %>% 
  filter(n_regulation == 2) %>% dim

papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%   
  filter(regulation %in% c("down")) %>% 
  select(hgnc_symbol, regulation) %>% unique() %>% dim
  

papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%   
  filter(regulation %in% c("up")) %>% 
  select(hgnc_symbol, regulation) %>% unique() %>% dim

# both 6230
# all down 9972, only down 3742
# all up 9455, only up 3225

# both 6230
# down 3742
# up 3225


# Załadowanie pakietu
library(VennDiagram)

# Definicja liczebności zbiorów
both <- 6230/2
down <- 3742
up <- 3225

# Ustalenie tej samej wielkości dla obu zbiorów
equal_size <- max(down + both, up + both)  # Można też użyć stałej wartości

dev.off()
svg("results/figures/reports-gr-database/venn-diagram-raw-up-donw.svg", width = 6, height = 6)

# Tworzenie diagramu Venna z równymi kołami
venn.plot <- draw.pairwise.venn(
  area1 = equal_size,   # Równa wielkość pierwszego zbioru
  area2 = equal_size,   # Równa wielkość drugiego zbioru
  cross.area = both,    # Wspólny obszar
  category = c("spadek poziomu\nekspresji genów", "wzorst poziomu\nekspresji genów"),
  fill = c("white", "white"),
  alpha = 0.1,
  cex = 2,
  cat.cex = 2,
  lwd = 2
)
  
dev.off()



papers_data_preprocessing %>% 
  # filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%   
  .$hgnc_symbol %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>%
  arrange(desc(freq)) %>% 
  .$freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("freq_genes", "n_genes")) %>% 
  mutate(propability = n_genes/13251)  %>% 
  mutate(cumsum = cumsum(propability)) %>% 
  mutate(freq_genes = ifelse(cumsum > 0.95, "more_than_12", as.character(freq_genes))) %>%
  group_by(freq_genes) %>% 
  summarise(n_genes = sum(n_genes), propability = sum(propability)) %>%
  ungroup() %>%
  mutate(freq_genes = factor(freq_genes, levels = c(as.character(1:12), "more_than_12"))) %>%
  arrange(freq_genes) %>% 
  mutate(cumsum = cumsum(propability)) -> data

data

# Otwarcie urządzenia graficznego do zapisu w formacie SVG
svg("results/figures/reports-gr-database/piechart-occurence-genes-database.svg", width = 8, height = 8) 

# Tworzenie wykresu
ggplot(data, aes(x = "", y = n_genes, fill = freq_genes)) +
  geom_bar(stat = "identity", width = 1, color = "black") +
  coord_polar("y", start = 0) +
  theme_void() +
  scale_fill_brewer(palette = "Paired") +
  labs(fill = "Gene Frequency") +
  ggtitle("Distribution of Gene Frequency")

# Zamknięcie urządzenia graficznego i zapis pliku
dev.off()


papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  # select(hgnc_symbol, treatment_type, source, dose, treatment, time) %>% unique() %>% 
  filter(treatment != "GRKO") %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  select(hgnc_symbol, treatment_type) %>% 
  na.omit() %>% 
  unique() %>% 
  .$treatment_type %>% 
  table
  

  

