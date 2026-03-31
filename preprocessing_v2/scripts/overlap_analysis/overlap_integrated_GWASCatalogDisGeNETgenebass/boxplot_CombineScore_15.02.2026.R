# ##############################################################################
# ---- p < 0.05 max & n overlap genes, p > 0.05 & max n overlap genes ----
# ##############################################################################


library(tidyverse)
boxplot_df <-
  GWASCatalogDisGeNETgenebass_GrSystemic_overlapChi2$processed$original_data$df %>% 
  mutate(log10_pvalue = -log10(p_value)) %>% 
  mutate(log10_geneOverlap = log10(gene_overlap_count + 1)) %>% 
  mutate(combine_score = log10_pvalue * log10_geneOverlap) %>% 
  arrange(desc(combine_score)) %>% 
  filter(!(Var1 %in% c("DisGeNET_NA", "genebass_NA", "GWASCatalog_NA",
                       "genebass_F9x_Date_F99_first_reported_(mental_disorder,_not_otherwise_specified)"))) %>% 
  
  rowwise() %>%
  mutate(
    matched = list(str_extract_all(Var1, paste(patterns, collapse="|"))[[1]]),
    n_unique_patterns = n_distinct(matched)
  ) %>%
  ungroup() %>%
  filter(n_unique_patterns < 2, n_unique_patterns != 0) %>% 
  
  mutate(
    icd10_category = map_chr(matched, ~ if(length(.x)==1) .x else NA_character_)
  ) %>% 
  select(-c(fdr, fdr_value, matched, n_unique_patterns)) %>% 
  
  group_by(icd10_category) %>% 
  slice_max(order_by = combine_score, n = 20, with_ties = FALSE) %>%  # ↑ możesz zmienić ile topów
  ungroup() %>% 
  
  filter(combine_score > 0)


order_levels <-
  boxplot_df %>%
  group_by(icd10_category) %>%
  summarise(mean_cs = mean(combine_score, na.rm = TRUE)) %>%
  arrange(desc(mean_cs)) %>%
  pull(icd10_category)

boxplot_df$icd10_category <- factor(boxplot_df$icd10_category,
                                    levels = order_levels)

order_levels <-
  boxplot_df %>%
  group_by(icd10_category) %>%
  summarise(mean_cs = mean(combine_score, na.rm = TRUE)) %>%
  arrange(desc(mean_cs)) %>%
  pull(icd10_category)

boxplot_df$icd10_category <- factor(boxplot_df$icd10_category,
                                    levels = rev(order_levels))  # rev → najwyżej u góry

library(ggplot2)

library(ggplot2)

p <- ggplot(boxplot_df,
            aes(y = icd10_category,
                x = combine_score)) +
  
  geom_boxplot(width = 0.6,
               fill = "#552c17ff",
               alpha = 0.6,
               outlier.shape = NA) +
  
  geom_jitter(height = 0.15,
              size = 1.6,
              alpha = 0.5,
              color = "black") +
  
  scale_x_continuous(position = "top") +   # <- oś X na górze
  
  labs(y = "ICD-10 category",
       x = "Combined score") +
  
  theme_classic(base_size = 14)

p

dev.off()
out_dir  <- "/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_GWASCatalogDisGeNETgenebass/figures"
out_file <- paste0("boxplot_combinedScore_grSystemic_", "15.02.2026", ".svg")

svg(filename = file.path(out_dir, out_file), width = 5, height = 4)
print(p)
dev.off()


p2 <- ggplot(boxplot_df,
             aes(y = icd10_category,
                 x = combine_score,
                 fill = Var2)) +
  
  geom_boxplot(width = 0.65,
               alpha = 0.7,
               outlier.shape = NA,
               position = position_dodge(width = 0.75)) +
  
  geom_point(aes(color = Var2),
             size = 1.4,
             alpha = 0.45,
             position = position_jitterdodge(
               jitter.height = 0.15,
               dodge.width = 0.75
             )) +
  
  scale_x_continuous(position = "top") +
  
  labs(y = "ICD-10 category",
       x = "Combined score",
       fill = "Source",
       color = "Source") +
  
  theme_classic(base_size = 14) +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal"
  )

p2

install.packages("rstatix")   # tylko raz
library(rstatix)

wilcox_results <-
  boxplot_df %>%
  group_by(icd10_category) %>%
  t_test(combine_score ~ Var2)

wilcox_results

boxplot_df %>% 
  filter(icd10_category == "F2x") %>% 
  filter(Var2 == "global_GR_genes_globalUp5TissuesDerivedCells") %>% 
  .$combine_score

t.test(c(7.21908514, 4.18649616, 2.03946314, 1.11319019, 0.40649980, 0.34811396, 0.32195109, 0.15818486, 0.13253126, 0.13253126, 0.09669748),
       c(3.87413349, 0.71970180, 0.48617988, 0.34054680, 0.33944648, 0.33828592, 0.26519477, 0.08989490, 0.08792266))

