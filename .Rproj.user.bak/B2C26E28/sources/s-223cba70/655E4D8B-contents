# dane do wykresu venna, dla wszystkich genów które są up lub down
filt_1_12_gr_database %>% 
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


################################################################################  
# dane do wykresu venna, dla wszystkich genów które mają wartości log2ratio
filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(!is.na(log2ratio)) %>% 
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


# Przykładowe geny FKBP5 (barplot?)
filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(hgnc_symbol == "FKBP5") %>%
  mutate(log2ratio = as.numeric(log2ratio)) %>%  # Konwersja log2ratio do numerycznego
  ggplot(aes(x = regulation, y = log2ratio, color = regulation)) +
  geom_boxplot(outlier.shape = NA, size = 1) +  # Pogrubienie boxplotu
  geom_jitter(width = 0.2, alpha = 1, color = "black") +  # Punkty na czarno
  theme_minimal() +
  labs(title = "Ekspresja FKBP5 (log2ratio) w zależności od regulation",
       x = "Regulation",
       y = "log2ratio") +
  scale_color_manual(values = c("down" = "blue4", "up" = "firebrick")) +
  theme(legend.position = "none")
  
# Zdefiniuj geny, które chcesz uwzględnić

filt_1_12_gr_database %>% 
  select(hgnc_symbol, hgnc_occurence) %>% 
  unique %>% 
  filter(hgnc_occurence > 35) %>% 
  na.omit %>% 
  .$hgnc_symbol -> genes_of_interest 

genes_of_interest <- c("FKBP5", "GENE2", "GENE3")  # Dodaj tutaj interesujące Cię geny

filt_1_12_gr_database %>% 
  filter(regulation %in% c("down", "up")) %>% 
  filter(!log2ratio %in% c("", "NA", "Inf", "-Inf")) %>%
  filter(hgnc_occurence > 35) %>%  
  filter(!is.na(hgnc_symbol)) %>% 
  mutate(log2ratio = as.numeric(log2ratio)) %>% 
  ggplot(aes(x = hgnc_symbol, y = log2ratio)) +
  geom_boxplot(outlier.shape = NA, size = 1) +  # Pogrubiony boxplot
  geom_jitter(width = 0.2, alpha = 1, color = "black") +  # Punkty
  theme_minimal() +
  labs(title = "Ekspresja genów (log2ratio)",
       x = "Gen",
       y = "log2ratio")






