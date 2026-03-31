papers_data_preprocessing %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%   
  filter(regulation %in% c("down", "up")) %>% 
  select(c(hgnc_symbol, regulation, log2ratio)) %>% 
  unique %>% 
  # filter(log2ratio != "NA") %>%
  select(c(hgnc_symbol, regulation)) %>% 
  unique() %>% 
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~nrow(.x))) %>% 
  unnest(n_regulation) %>% 
  filter(n_regulation == 2) %>% 
  .$hgnc_symbol -> bidirectionally_regulated_genes

bidirectionally_regulated_genes %>% length()

papers_data_preprocessing %>% 
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>% 
  select(regulation, log2ratio) %>%  
  filter(regulation %in% c("up", "down")) %>% 
  .$regulation %>%
  table

papers_data_preprocessing %>% 
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>% 
  select(regulation, log2ratio) %>% 
  na.omit() %>%
  filter(log2ratio != "NA") %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>%
  # mutate(log2ratio = ifelse(log2ratio == Inf, max_log2ratio_value, log2ratio)) %>% 
  # mutate(log2ratio = ifelse(log2ratio == -Inf, min_log2ratio_value, log2ratio)) %>% 
  ggplot(aes(x = regulation, y = log2ratio)) +
  geom_boxplot() +
  # geom_jitter(size = 0.5)
  theme_minimal() +
  labs(x = "Regulacja genów", y = "Log2(FC)")



max_log2ratio_value <- 22
min_log2ratio_value <- -16


papers_data_preprocessing %>% 
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>% 
  select(hgnc_symbol, regulation, log2ratio) %>% 
  na.omit() %>%
  filter(log2ratio != "NA") %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  group_by(hgnc_symbol) %>% 
  nest %>% 
  mutate(max_log2ratio = map(data, ~max(.x$log2ratio))) %>% 
  mutate(min_log2ratio = map(data, ~min(.x$log2ratio))) %>% 
  unnest(c(max_log2ratio, min_log2ratio)) %>% 
  mutate(max_log2ratio = ifelse(max_log2ratio == Inf, max_log2ratio_value, max_log2ratio)) %>% 
  mutate(min_log2ratio = ifelse(min_log2ratio == -Inf, min_log2ratio_value, min_log2ratio)) %>% 
  mutate(diff_regulation = max_log2ratio - min_log2ratio) %>% 
  arrange(desc(diff_regulation)) %>% head(20) %>% tail


papers_data_preprocessing %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% 
  select(c(hgnc_symbol, regulation)) %>% 
  unique %>%
  group_by(hgnc_symbol) %>% 
  nest() %>% 
  mutate(n_regulation = map(data, ~nrow(.x))) %>% 
  unnest(n_regulation) %>% 
  filter(n_regulation == 2) %>% 
  .$hgnc_symbol -> bidirectionally_regulated_genes

papers_data_preprocessing %>% 
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>% 
  # select(hgnc_symbol, regulation, log2ratio, time) %>% 
  # na.omit() %>%
  # filter(log2ratio != "NA") %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  mutate(log2ratio = ifelse(log2ratio == Inf, max_log2ratio_value, log2ratio)) %>% 
  mutate(log2ratio = ifelse(log2ratio == -Inf, min_log2ratio_value, log2ratio)) -> data_to_check_regulation
  


data_to_check_regulation %>% 
  filter(regulation == "up") %>% .$time %>% table

data_to_check_regulation %>% 
  filter(regulation == "down") %>% .$time %>% table

data_to_check_regulation %>% 
  filter(regulation == "up") %>% .$simple_tissue %>% table

data_to_check_regulation %>% 
  filter(regulation == "down") %>% .$simple_tissue %>% table

data_to_check_regulation %>%
  count(simple_tissue, regulation) %>%  # Zlicza liczbę przypadków dla każdej kombinacji
  ggplot(aes(x = simple_tissue, y = n, fill = regulation)) +
  geom_col(position = "dodge") +  # Słupki obok siebie dla up/down
  theme_minimal() +
  labs(
       x = "Tkanina",
       y = "Liczba wystąpień",
       fill = "Regulacja") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




gene_tissue_regulation <- papers_data_preprocessing %>% 
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>%
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% 
  group_by(hgnc_symbol) %>%
  summarise(
    n_tissues_up = n_distinct(simple_tissue[regulation == "up"]),
    n_tissues_down = n_distinct(simple_tissue[regulation == "down"]),
    same_tissue = any(simple_tissue[regulation == "up"] %in% simple_tissue[regulation == "down"]),
    different_tissue = any(!simple_tissue[regulation == "up"] %in% simple_tissue[regulation == "down"]) &
      any(!simple_tissue[regulation == "down"] %in% simple_tissue[regulation == "up"])
  ) %>%
  mutate(
    regulation_type = case_when(
      same_tissue & !different_tissue ~ "same_tissue_only",
      different_tissue & !same_tissue ~ "different_tissue_only",
      same_tissue & different_tissue ~ "both_same_and_different",
      TRUE ~ "unknown" # Na wszelki wypadek, gdyby jakiś gen nie pasował
    )
  )


gene_tissue_regulation %>% .$regulation_type %>% table

gene_tissue_regulation %>%
  summarise(
    same_tissue_count = sum(same_tissue, na.rm = TRUE),
    different_tissue_count = sum(different_tissue, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = everything(), names_to = "Typ regulacji", values_to = "Liczba genów")


papers_data_preprocessing %>% 
  filter(regulation %in% c("up", "down")) %>%
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>%
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% .$hgnc_symbol %>% unique()
  
papers_data_preprocessing %>% 
  filter(regulation %in% c("up", "down")) %>%
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>%
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% 
  group_by(hgnc_symbol) %>% 
  summarise(
    n_treatments_up = n_distinct(treatment[regulation == "up"]),
    n_treatments_down = n_distinct(treatment[regulation == "down"]),
    same_treatment = any(treatment[regulation == "up"] %in% treatment[regulation == "down"]),
    different_treatment = any(!treatment[regulation == "up"] %in% treatment[regulation == "down"]) &
      any(!treatment[regulation == "down"] %in% treatment[regulation == "up"])
  ) %>%
  mutate(
    regulation_type_treatment = case_when(
      same_treatment & !different_treatment ~ "same_treatment_only",
      different_treatment & !same_treatment ~ "different_treatment_only",
      same_treatment & different_treatment ~ "both_same_and_different",
      TRUE ~ "unknown" # Na wszelki wypadek, gdyby jakiś gen nie pasował
    )
  )


papers_data_preprocessing %>% 
  filter(regulation %in% c("up", "down")) %>%
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>%
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>% 
  mutate(time = recode(time, 
                       "1" = "1h", 
                       "10" = "10h", 
                       "2" = "2h", 
                       "4" = "4h", 
                       "24" = "24h", 
                       "3m" = "3months",
                       "8weeks" = "8weeks",
                       "18h" = "18h")) %>%  # Dodane "18h"
  group_by(hgnc_symbol) %>% 
  summarise(
    # Liczba unikalnych treatment, tissue, time i dose dla up/down
    n_treatments_up = n_distinct(treatment[regulation == "up"]),
    n_treatments_down = n_distinct(treatment[regulation == "down"]),
    n_tissues_up = n_distinct(simple_tissue[regulation == "up"]),
    n_tissues_down = n_distinct(simple_tissue[regulation == "down"]),
    n_times_up = n_distinct(time[regulation == "up"]),
    n_times_down = n_distinct(time[regulation == "down"]),
    n_doses_up = n_distinct(dose[regulation == "up"]),
    n_doses_down = n_distinct(dose[regulation == "down"]),
    
    # Sprawdzenie, czy wszystkie treatment/tissue/time/dose są takie same w up i down
    same_treatment = all(treatment[regulation == "up"] %in% treatment[regulation == "down"]) &
      all(treatment[regulation == "down"] %in% treatment[regulation == "up"]),
    
    same_tissue = all(simple_tissue[regulation == "up"] %in% simple_tissue[regulation == "down"]) &
      all(simple_tissue[regulation == "down"] %in% simple_tissue[regulation == "up"]),
    
    same_time = all(time[regulation == "up"] %in% time[regulation == "down"]) &
      all(time[regulation == "down"] %in% time[regulation == "up"]),
    
    same_dose = all(dose[regulation == "up"] %in% dose[regulation == "down"]) &
      all(dose[regulation == "down"] %in% dose[regulation == "up"]),
    
    # Sprawdzenie częściowej zgodności
    partial_treatment = any(treatment[regulation == "up"] %in% treatment[regulation == "down"]) |
      any(treatment[regulation == "down"] %in% treatment[regulation == "up"]),
    
    partial_tissue = any(simple_tissue[regulation == "up"] %in% simple_tissue[regulation == "down"]) |
      any(simple_tissue[regulation == "down"] %in% simple_tissue[regulation == "up"]),
    
    partial_time = any(time[regulation == "up"] %in% time[regulation == "down"]) |
      any(time[regulation == "down"] %in% time[regulation == "up"]),
    
    partial_dose = any(dose[regulation == "up"] %in% dose[regulation == "down"]) |
      any(dose[regulation == "down"] %in% dose[regulation == "up"])
  ) %>%
  mutate(
    regulation_type_treatment = case_when(
      same_treatment ~ "same_treatment_only",
      partial_treatment ~ "both_same_and_different_treatment",
      TRUE ~ "different_treatment"
    ),
    regulation_type_tissue = case_when(
      same_tissue ~ "same_tissue_only",
      partial_tissue ~ "both_same_and_different_tissue",
      TRUE ~ "different_tissue"
    ),
    regulation_type_time = case_when(
      same_time ~ "same_time_only",
      partial_time ~ "both_same_and_different_time",
      TRUE ~ "different_time"
    ),
    regulation_type_dose = case_when(
      same_dose ~ "same_dose_only",
      partial_dose ~ "both_same_and_different_dose",
      TRUE ~ "different_dose"
    )
  ) -> comparision_up_down

comparision_up_down %>% dim
  filter(regulation_type_tissue %in% c("same_tissue_only")) %>% 
  filter(regulation_type_treatment %in% c("same_treatment_only")) %>% 
  filter(regulation_type_time %in% c("same_time_only")) %>% .$hgnc_symbol -> genes_to_check

comparision_up_down %>% 
  select(hgnc_symbol, regulation_type_treatment, regulation_type_tissue, regulation_type_time) %>% 
  mutate(label = paste0(regulation_type_treatment, "_", regulation_type_tissue, "_", regulation_type_time)) %>% 
  .$regulation_type_treatment %>% table

comparision_up_down %>% 
  select(hgnc_symbol, regulation_type_treatment, regulation_type_tissue, regulation_type_time, regulation_type_dose) %>% 
  mutate(label = paste0(regulation_type_treatment, "_", regulation_type_tissue, "_", regulation_type_time)) %>% 
  .$regulation_type_time %>% table

comparision_up_down %>% 
  select(hgnc_symbol, regulation_type_treatment, regulation_type_tissue, regulation_type_time, regulation_type_dose) %>% 
  .$regulation_type_dose %>% table 

comparision_up_down %>% 
  select(hgnc_symbol, regulation_type_treatment, regulation_type_tissue, regulation_type_time, regulation_type_dose) %>% 
  .$regulation_type_tissue %>% table


comparision_up_down %>% 
  filter(regulation_type_tissue == "different_tissue" | regulation_type_treatment == "different_type_treatment" | regulation_type_time == "different_time" | regulation_type_dose == "different_dose") %>% dim

  
comparision_up_down %>% 
  filter(regulation_type_tissue == "same_tissue_only" & regulation_type_treatment == "same_treatment_only" & regulation_type_time == "same_time_only" & regulation_type_dose == "same_dose_only") %>% dim
  .$hgnc_symbol -> genes_to_check
  
  
    
comparision_up_down %>% dim
  filter(regulation_type_tissue != "different_tissue" & regulation_type_treatment != "different_type_treatment" & regulation_type_time != "different_time") %>% dim

comparision_up_down %>% 
  filter(regulation_type_tissue != "different_tissue" & regulation_type_treatment != "different_type_treatment" & regulation_type_time != "different_time") 
  
  



papers_data_preprocessing %>% 
  filter(regulation %in% c("up", "down")) %>%
  filter(hgnc_symbol %in% bidirectionally_regulated_genes) %>%
  filter(!(treatment %in% c("NA", "PHA", "formoterol", "aldosterone", "vehicle-DMSO", "vitamin-d3", "mifepristone", "eplerenone"))) %>%
  mutate(time = ifelse(time == "1", "1h", time)) %>% 
  mutate(time = ifelse(time == "10", "10h", time)) %>% 
  mutate(time = ifelse(time == "2", "2h", time)) %>% 
  mutate(time = ifelse(time == "4", "4h", time)) %>% 
  mutate(time = ifelse(time == "24", "24h", time)) %>% 
  mutate(time = ifelse(time == "3m", "3months", time)) %>% 
  filter(hgnc_symbol %in% genes_to_check) %>% .$regulation %>% table
  select(c(hgnc_symbol, gene_name, source, tissue, cell, log2ratio, fdr,  species, environment, treatment, dose, time, treatment_type, regulation, simple_tissue, dose)) %>%
  arrange(hgnc_symbol) %>% view

