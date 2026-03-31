papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% # clean treatment 
  mutate(treatment = ifelse(treatment == "corticoterone", "corticosterone", treatment)) %>% 
  mutate(treatment_type = ifelse(time == "3weeks", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time == "3m", "chronic", treatment_type)) %>% 
  mutate(treatment_type = ifelse(time %in% c("1h", "2h", "3h", "4h", "5h", "6h", "12h", "18h", "24h") & is.na(treatment_type), "acute", treatment_type)) %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(hgnc_occurence = map(data, ~nrow(.x))) %>% 
  unnest(hgnc_occurence) %>% 
  filter(hgnc_occurence != 1) %>% 
  unnest(data) %>% 
  ungroup() %>% 
  as.data.frame() -> filt_1_gr_database

  