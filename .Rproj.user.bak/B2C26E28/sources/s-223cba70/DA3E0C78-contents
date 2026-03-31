# dane do wykresu venna, dla wszystkich genów które są up lub down
filt_1_gr_database %>% 
  # filter(hgnc_occurence > 8) %>% 
  filter(regulation %in% c("down", "up")) %>%
  # filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  select(hgnc_symbol, regulation) %>% 
  unique %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~nrow(.x))) %>% 
  unnest() %>% 
  unique() %>% 
  mutate(group = ifelse(n_regulation == 2, "up_down", regulation)) %>% 
  select(hgnc_symbol, group) %>% 
  unique %>% 
  .$group %>% 
  table

venn_df <- filt_1_gr_database %>%
  filter(regulation %in% c("down", "up")) %>%
  select(hgnc_symbol, regulation)

# Tworzymy listy genów dla poszczególnych regulacji:
up_genes <- venn_df %>% 
  filter(regulation == "up") %>% 
  pull(hgnc_symbol) %>% 
  unique()

down_genes <- venn_df %>% 
  filter(regulation == "down") %>% 
  pull(hgnc_symbol) %>% 
  unique()

# Przygotowanie listy dla funkcji venn.diagram
venn_list <- list(Up = up_genes, Down = down_genes)

# Tworzenie wykresu Venn’a
venn_plot <- venn.diagram(
  x = venn_list,               # lista zawierająca zbiory genów
  filename = NULL,             # nie zapisujemy do pliku, tylko tworzymy obiekt graficzny
  fill = c("red", "blue"),     # kolory wypełnienia poszczególnych kół
  alpha = c(0.5, 0.5),         # przezroczystość kolorów
  cex = 2,                     # rozmiar tekstu w obszarach wykresu
  cat.cex = 2,                 # rozmiar etykiet kategorii
  cat.col = c("red", "blue")   # kolory etykiet kategorii
)

# Wyświetlenie wykresu Venn’a
grid.newpage()
grid.draw(venn_plot)  


################################################################################  
# dane do wykresu venna, dla wszystkich genów które mają wartości log2ratio
filt_1_gr_database %>% 
  # filter(hgnc_occurence > 8) %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  select(hgnc_symbol, regulation) %>% 
  unique %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~nrow(.x))) %>% 
  unnest() %>% 
  unique() %>% 
  mutate(group = ifelse(n_regulation == 2, "up_down", regulation)) %>% 
  select(hgnc_symbol, group) %>% 
  unique %>% 
  .$group %>% 
  table

# obserwacja/wniosek, im więcej genów się uwzględni tym więcej genów ulega aktywacji i hamowaniu


venn_df <- filt_1_gr_database %>%
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  select(hgnc_symbol, regulation)

# Tworzymy listy genów dla poszczególnych regulacji:
up_genes <- venn_df %>% 
  filter(regulation == "up") %>% 
  pull(hgnc_symbol) %>% 
  unique()

down_genes <- venn_df %>% 
  filter(regulation == "down") %>% 
  pull(hgnc_symbol) %>% 
  unique()

# Przygotowanie listy dla funkcji venn.diagram
venn_list <- list(Up = up_genes, Down = down_genes)

# Tworzenie wykresu Venn’a
venn_plot <- venn.diagram(
  x = venn_list,               # lista zawierająca zbiory genów
  filename = NULL,             # nie zapisujemy do pliku, tylko tworzymy obiekt graficzny
  fill = c("red", "blue"),     # kolory wypełnienia poszczególnych kół
  alpha = c(0.5, 0.5),         # przezroczystość kolorów
  cex = 2,                     # rozmiar tekstu w obszarach wykresu
  cat.cex = 2,                 # rozmiar etykiet kategorii
  cat.col = c("red", "blue")   # kolory etykiet kategorii
)

# Wyświetlenie wykresu Venn’a
grid.newpage()
grid.draw(venn_plot)  


