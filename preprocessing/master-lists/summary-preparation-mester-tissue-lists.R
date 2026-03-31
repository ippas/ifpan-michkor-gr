################################################################################
# prepare data frame with parameters to prepare master lists for tissues
data.frame(tissue = papers_data_preprocessing %>% 
             filter(!(simple_tissue %in% c("Other"))) %>% 
             .$simple_tissue %>% unique() %>% rep(each = 8),
           metric = rep(c("log2ratio", "fdr"), 56),
           sort_order = rep(c("decrease", "increase"), 56),
           top = rep(rep(c(10, 25, 50, 100), each = 2), 14)
) -> parameters_df_tissues

rbind({parameters_df_tissues %>% mutate(regulation = "up")}, 
      {parameters_df_tissues %>% mutate(regulation = "down")}) %>% 
  mutate(sort_order = ifelse(regulation == "down", "increase", sort_order)) -> parameters_df_tissues


parameters_df_tissues %>% 
  filter(top == 50) %>% 
  filter(metric == "log2ratio") -> parameters_df_tissues_log2ratio_top50


lapply(1:nrow(parameters_df_tissues_log2ratio_top50), function(i){
  print(parameters_df_tissues_log2ratio_top50[i, ])
  
  tryCatch({
    create_ranked_master_list(data =  papers_data_preprocessing,
                              arrange_by = parameters_df_tissues_log2ratio_top50[i, "metric"],
                              simple_tissue == parameters_df_tissues_log2ratio_top50[i, "tissue"],
                              regulation == parameters_df_tissues_log2ratio_top50[i, "regulation"],
                              !(treatment %in% treatment_to_remove),
                              columns = c(
                                "source",
                                "tissue",
                                "cell",
                                "dose",
                                "treatment",
                                "treatment_type",
                                "regulation",
                                "comparison",
                                "environment"
                              ),
                              keep_column =  c("simple_tissue", "method"),
                              size_list = parameters_df_tissues_log2ratio_top50[i, "top"],
                              max_rank = parameters_df_tissues_log2ratio_top50[i, "top"],
                              top_n = parameters_df_tissues_log2ratio_top50[i, "top"],
                              sort_order = parameters_df_tissues_log2ratio_top50[i, "sort_order"]
    )
  }, error = function(e) NULL)
  
}) -> tissue_master_lists_log2ratio_top50

# set names for results
names(tissue_master_lists_log2ratio_top50) <- parameters_df_tissues_log2ratio_top50 %>% mutate(names = paste(tissue, regulation, top, metric, sep = "_")) %>% .$names 

# assessment after elimination
lapply(tissue_master_lists_log2ratio_top50, function(x){
  x$gene_lists 
}) %>% unlist(recursive = F) %>% 
  unname() %>% 
  unlist %>% 
  unique() %>% length()

lapply(tissue_master_lists_log2ratio_top50, function(x){
  x$gene_lists 
}) %>% unlist(recursive = F) %>% 
  unname() %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  .$Freq %>% 
  table %>%
  as.data.frame() %>% 
  set_colnames(c("times_occurence", "freq")) %>% 
  mutate(percent = freq/10901) %>% 
  mutate(times_occurence = as.character(times_occurence)) %>%
  mutate(times_occurence = as.numeric(times_occurence)) %>% 
  filter(times_occurence > 9) %>% .$freq %>% sum


# summary tissue gene list
lapply(tissue_master_lists_log2ratio_top50, function(x){
  x$master_df$hgnc_symbol
}) %>% 
  unname() %>% unlist %>% unique() %>% length()


lapply(tissue_master_lists_log2ratio_top50, function(x){
  x$master_df$hgnc_symbol
}) %>% 
  unname() %>% unlist %>%
  table %>% 
  as.data.frame() %>% 
  .$Freq %>% 
  table %>% 
  as.data.frame() %>% .$Freq %>% cat(sep = "\n")


# Ładujemy pakiet ggplot2
library(ggplot2)

# Tworzymy ramkę danych
data <- data.frame(
  frequency = c("once", "twice", "three", "more"),
  number_of_genes = c(872, 133, 26, 15)
)

# Tworzymy ramkę danych z określoną kolejnością poziomów dla osi x
data <- data.frame(
  frequency = factor(c("once", "twice", "three", "more"), levels = c("once", "twice", "three", "more")),
  number_of_genes = c(872, 133, 26, 15)
)

ggplot(data, aes(x = frequency, y = number_of_genes)) +
  geom_bar(stat = "identity", fill = "gray70", width = 0.65) + # kontrolujemy szerokość słupków
  labs(x = "Frequency of Genes", y = "Number of Genes") +
  theme_classic(base_size = 25) +
  theme(
    axis.title.x = element_text(size = 100, margin = margin(t = 50)), # odstęp nad tytułem osi X
    axis.title.y = element_text(size = 100, margin = margin(r = 50)), # odstęp po prawej od tytułu osi Y
    axis.text = element_text(size = 80)
  )


lapply(tissue_master_lists_log2ratio_top50, function(x){
  x$master_df$hgnc_symbol
}) %>% 
  unname() %>% unlist %>%
  table %>% 
  as.data.frame() %>% 
  filter(Freq > 3)


################################################################################
# charakterystyka bazy, bez filtrowania list
lapply(1:nrow(parameters_df_tissues_log2ratio_top50), function(i){
  print(parameters_df_tissues_log2ratio_top50[i, ])
  
  tryCatch({
    create_ranked_master_list(data =  papers_data_preprocessing,
                              arrange_by = parameters_df_tissues_log2ratio_top50[i, "metric"],
                              simple_tissue == parameters_df_tissues_log2ratio_top50[i, "tissue"],
                              regulation == parameters_df_tissues_log2ratio_top50[i, "regulation"],
                              # !(treatment %in% treatment_to_remove),
                              columns = c(
                                "source",
                                "tissue",
                                "cell",
                                "dose",
                                "treatment",
                                "treatment_type",
                                "regulation",
                                "comparison",
                                "environment"
                              ),
                              keep_column =  c("simple_tissue", "method"),
                              size_list = 0,
                              max_rank = parameters_df_tissues_log2ratio_top50[i, "top"],
                              top_n = parameters_df_tissues_log2ratio_top50[i, "top"],
                              sort_order = parameters_df_tissues_log2ratio_top50[i, "sort_order"]
    )
  }, error = function(e) NULL)
  
}) -> tissue_master_lists_log2ratio_all

# set names for results
names(tissue_master_lists_log2ratio_all) <-  parameters_df_tissues_log2ratio_top50  %>% mutate(names = paste(tissue, regulation, sep = "_")) %>% .$names 


tissue_master_lists_log2ratio_all$brain_up$gene_lists


lapply(tissue_master_lists_log2ratio_all, function(x){
  x$gene_lists 
}) %>% unlist(recursive = F) %>% 
  unname() %>% unlist %>% unique() %>% length()


# raw list
papers_data_preprocessing %>% 
  # drop_na(log2ratio) %>%
  # filter(log2ratio != "NA") %>%
  # filter(!(treatment %in% treatment_to_remove)) %>%
  # filter(regulation == "up") %>%
  refine_gene_lists(data = .,
                    columns = c(
                      "source",
                      "tissue",
                      "cell",
                      "dose",
                      "treatment",
                      "treatment_type",
                      "regulation",
                      "comparison",
                      "environment"
                    ),
                    cumsum_thresholds = cumsum_thresholds,
                    freq_thresholds = freq_thresholds,
                    keep_column =  c("simple_tissue", "method"))  -> raw_list


raw_list$refine_gene_lists$gene_lists %>% 
  unname() %>% 
  unlist %>% 
  table %>% 
  as.data.frame() %>% 
  .$Freq %>% 
  table %>% 
  as.data.frame() %>% 
  set_colnames(c("gene_frequency", "number_of_genes")) %>% 
  mutate(proportion = number_of_genes / 13610) %>% 
  mutate(gene_frequency = as.character(gene_frequency)) %>% 
  mutate(gene_frequency = as.numeric(gene_frequency)) %>% 
  mutate(gene_frequency = ifelse(gene_frequency > 9, "more", as.character(gene_frequency))) %>% 
  group_by(gene_frequency) %>% 
  summarise(
    number_of_genes = sum(number_of_genes),
    proportion = sum(proportion)
  ) %>% 
  as.data.frame() -> data

multiple <- 5
data %>%
  mutate(
    gene_frequency2 = factor(
      c("once", "twice", "three", "four", "five", "six", "seven", "eight", "nine", "more"),
      levels = c("once", "twice", "three", "four", "five", "six", "seven", "eight", "nine", "more")
    )
  ) %>% 
  ggplot(aes(x = gene_frequency2, y = number_of_genes)) +
  geom_bar(stat = "identity", fill = "gray70", width = 0.65) +
  labs(x = "Frequency of Genes", y = "Number of Genes") +
  theme_classic(base_size = 25) +
  theme(
    axis.title.x = element_text(size = 10 * multiple, margin = margin(t = 0 * multiple)),
    axis.title.y = element_text(size = 10 * multiple, margin = margin(r = 5 * multiple)),
    axis.text = element_text(size = 8 * multiple),
    axis.text.x = element_text(angle = 45, hjust = 1) # kąt 45 stopni dla tekstu osi X
  )

