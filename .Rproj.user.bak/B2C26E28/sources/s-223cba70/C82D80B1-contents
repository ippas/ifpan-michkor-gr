papers_data_preprocessing2 %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% # clean treatment 
  mutate(treatment = ifelse(treatment == "corticoterone", "corticosterone", treatment)) %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(hgnc_occurence = map(data, ~nrow(.x))) %>% 
  unnest(hgnc_occurence) %>% 
  filter(hgnc_occurence > 12) %>% 
  filter(!is.na(hgnc_symbol)) %>% 
  unnest(data) %>% 
  # head %>% 
  # filter(source != "pmid:34272384") %>% 
  # filter(source != "pmid:24926665") %>% 
  # filter(source != "pmid:36699046") %>% 
  ungroup() %>% 
  mutate(label = paste(label, treatment, dose, time, treatment_type, environment, comparison, sep = "_")) %>%
  as.data.frame() -> filt_1_12_gr_database

filt_1_12_gr_database %>% 
  filter(is.na(log2ratio))



filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>%
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(!is.na(log2ratio)) %>% 
  group_by(hgnc_symbol, hgnc_occurence) %>% 
  nest() %>% 
  mutate(log2ratio_hgnc_occurence = map(data, ~nrow(.x))) %>% 
  unnest(log2ratio_hgnc_occurence) %>% 
  mutate(variance_expression = map(data, ~var(.x$log2ratio))) %>%
  mutate(mean_expression = map_dbl(data, ~ mean(as.numeric(as.character(.x$log2ratio)), na.rm = TRUE))) %>%
  unnest(variance_expression) %>%
  unnest(mean_expression) %>%
  mutate(cv = variance_expression/mean_expression) %>% 
  mutate(n_paper = map(data, ~length(unique(.x$source)))) %>% 
  mutate(n_tissues = map(data, ~length(unique(.x$simple_tissue))))  %>% 
  mutate(collapse_tissues = map(data, ~paste(unique(.x$simple_tissue), collapse = "|")))  %>%
  unnest(n_paper, n_tissues, collapse_tissues) %>% 
  mutate(hgnc_label = paste(hgnc_symbol,
                            log2ratio_hgnc_occurence,
                            n_paper,
                            n_tissues,
                            round(mean_expression, 2),
                            round(variance_expression, 2),
                            round(cv, 2),
                            sep = "\n")) %>%
  # mutate(hgnc_label = paste0(
  #   "bold('", hgnc_symbol, "')", "\n",
  #   log2ratio_hgnc_occurence, "\n",
  #   n_paper, "\n",
  #   n_tissues, "\n",
  #   round(mean_expression, 2), "\n",
  #   round(variance_expression, 2), "\n",
  #   round(cv, 2)
  # )) %>% 
  unnest(data) %>% 
  as.data.frame() -> filt_1_12_gr_database
  
