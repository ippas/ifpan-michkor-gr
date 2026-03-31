papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  filter(treatment != "GRKO") %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  select(hgnc_symbol, treatment_type) %>% 
  na.omit() %>% 
  unique()

# treatment_type venn diagram
# Twój przetworzony zbiór danych
venn_data <- papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  filter(treatment != "GRKO") %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  select(hgnc_symbol, treatment_type) %>% 
  na.omit() %>% 
  unique()

# Podział genów na trzy grupy
chronic_genes <- venn_data %>% filter(treatment_type == "chronic") %>% pull(hgnc_symbol)
acute_genes <- venn_data %>% filter(treatment_type == "acute") %>% pull(hgnc_symbol)
stress_genes <- venn_data %>% filter(treatment_type == "stress_induction") %>% pull(hgnc_symbol)

# Tworzenie listy dla wykresu Venn'a
venn_list <- list(
  Chronic = chronic_genes,
  Acute = acute_genes,
  Stress_Induction = stress_genes
)

# Rysowanie wykresu Venn'a
venn.plot <- venn.diagram(
  x = venn_list,
  filename = NULL,
  fill = c("white", "white", "white"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.5,
  cat.col = c("blue", "red", "green"),
  margin = 0.1
)

grid.newpage()
grid.draw(venn.plot)
svg("results/figures/reports-gr-database/venn-treatment-type.svg", width = 8, height = 8) 
grid.draw(venn.plot)
dev.off()


venn_list %>% unname() %>% unlist %>% table %>% as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  filter(freq == 3) %>% .$hgnc_symbol %>% as.character() -> genes_common_treatment_type

papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  filter(treatment != "GRKO") %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  filter(hgnc_symbol %in% genes_common_treatment_type) %>% 
  select(hgnc_symbol, regulation, treatment_type, tissue, dose, source, treatment, time) %>% 
  unique() %>% 
  .$hgnc_symbol %>% table %>% as.data.frame() %>% 
  set_colnames(c("hgnc_symbol", "freq")) %>% 
  arrange(desc(freq))

papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  filter(treatment != "GRKO") %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  filter(hgnc_symbol %in%  c("FKBP5", "KLF9", "PDK4", "TSC22D3", "DDIT4", "BCL6", "ZBTB16", "SGK1", "MERTK", "NFKBIA")) %>% 
  select(hgnc_symbol, source, treatment_type) %>% 
  arrange(hgnc_symbol)

papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%  
  filter(treatment != "GRKO") %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  filter(hgnc_symbol %in% genes_common_treatment_type) %>% 
  select(hgnc_symbol, regulation, treatment_type, info) -> tmp_data_treatment 
  
tmp_data_treatment %>% 
  filter(treatment_type == "chronic") %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(freq = map(data, ~nrow(.x))) %>% 
  unnest(freq) %>% 
  arrange(desc(freq)) %>% 
  filter(freq > 2) %>% 
  select(-data) %>% 
  mutate(type = "chronic") ->tmp_chronic


tmp_data_treatment %>% 
  filter(treatment_type == "acute") %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(freq = map(data, ~nrow(.x))) %>% 
  unnest(freq) %>% 
  arrange(desc(freq)) %>% 
  filter(freq > 20) %>% 
  select(-data) %>% 
  mutate(type = "acute") -> tmp_acute


tmp_data_treatment %>% 
  filter(treatment_type == "stress_induction") %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(freq = map(data, ~nrow(.x))) %>% 
  unnest(freq) %>% 
  arrange(desc(freq)) %>% 
  filter(freq > 2) %>%
  select(-data) %>% 
  mutate(type = "stress") %>% 
  as.data.frame() -> tmp_stress
  # filter(hgnc_symbol == "FKBP5")

rbind(tmp_acute, tmp_chronic, tmp_stress) %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_type = map(data, ~nrow(.x))) %>% unnest(n_type) %>% 
  arrange(desc(n_type)) %>% 
  select(-data) %>% 
  as.data.frame() %>% 
  filter(n_type > 1)

  
  # 
  # select(hgnc_symbol, regulation, treatment_type, tissue, dose, source, treatment, time) %>% 
  # unique() %>% 
  # .$hgnc_symbol %>% table %>% as.data.frame() %>% 
  # set_colnames(c("hgnc_symbol", "freq")) %>% 
  # arrange(desc(freq))

