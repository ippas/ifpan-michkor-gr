source("preprocessing/functions/R/install-load-packages.R")

read_genes_from_phenotype_models(directory_path = "data/prs-models-pan-biobank-uk/") -> genes_phenotypes_PanUkBiobank

################################################################################
# Read and preprocess the data
categories_biobankuk <- read.csv("data/phenotypes/panukbiobank-phenotype-category.csv") %>%
  filter(category_type == "320_Origin_Categories") %>%
  group_by(category_name) %>%
  nest() %>%
  mutate(n = map_int(data, nrow)) %>%
  filter(n >= 10, n <= 300) %>%
  pull(category_name) 


overlap_resutls_random_dependent_fdr0.01_new <- analyze_uk_biobank_categories_gene_overlap(
  categories = categories_biobankuk,
  genes_phenotypes = genes_phenotypes_PanUkBiobank,
  # papers_gene_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  papers_gene_list = papers_gene_list,
  hgnc_symbols_vector = hgnc_symbols_vector_v110,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  random_genes_type = "dependent",
  fdr_threshold = 0.01, # Default value, adjust as needed
  overlap_threshold = 3, # Default value, adjust as needed
  permutations = 1000, # Default value, adjust as needed
  seed = 123, # Default value, adjust as needed
  num_cores = 30 # Adjust based on your system's capabilities
)

overlap_resutls_random_independent_fdr0.05 <- analyze_uk_biobank_categories_gene_overlap(
  categories = categories_biobankuk,
  genes_phenotypes = genes_phenotypes_PanUkBiobank,
  papers_gene_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],
  hgnc_symbols_vector = hgnc_symbols_vector_v110,
  path_metafile = "data/phenotypes/panukbiobank-phenotype-category.csv",
  random_genes_type = "dependent",
  fdr_threshold = 0.05, # Default value, adjust as needed
  overlap_threshold = 3, # Default value, adjust as needed
  permutations = 1000, # Default value, adjust as needed
  seed = 123, # Default value, adjust as needed
  num_cores = 30 # Adjust based on your system's capabilities
)

pie_chart_palette <-c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A", "4" = "#984EA3", 
                      "5" = "#FF7F00", "more" = "#999999")
pie_chart_palette <-  generate_palette(6, "Pastel1")



generate_random_gene_list_dependent(genes_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))], 
                                    genes_vector = hgnc_symbols_vector_v110) -> random_genes_dependent

generate_random_gene_list_independent(genes_vector = hgnc_symbols_vector_v110, 
                                      size_vector = sapply(papers_gene_list[!grepl("marpiech", names(papers_gene_list))], length)) -> random_genes_independent


create_gene_pie_chart(gene_list = papers_gene_list[!grepl("marpiech", names(papers_gene_list))],threshold_percent = 95, 
                      palette = pie_chart_palette) -> p1

create_gene_pie_chart(gene_list = random_genes_dependent, threshold_percent = 95, 
                      palette = pie_chart_palette) -> p2


create_gene_pie_chart(gene_list = random_genes_independent, threshold_percent = 95,
                      palette = pie_chart_palette) -> p3


p1 + p3 + p2 + plot_layout(guides = "collect") & 
  plot_annotation(theme = theme(legend.position = "bottom"))


# example code testing
overlap_results$`Mental health`$significant_uniq_data$overlap_genes -> overlap_genes_mental_health

overlap_results$Depression$significant_uniq_data$overlap_genes -> overlap_genes_depression

overlap_results$`Psychosocial factors`$significant_uniq_data$overlap_genes -> overlap_genes_psychosocial_factors

overlap_genes_psychosocial_factors %>% length()

intersect(overlap_genes_depression, overlap_genes_psychosocial_factors) %>% length()

intersect(overlap_genes_mental_health, overlap_genes_psychosocial_factors) %>% length()

lapply(overlap_results, function(x){intersect(x$significant_uniq_data$overlap_gene, overlap_genes_psychosocial_factors) %>% length()/length(overlap_genes_psychosocial_factors)})

lapply(overlap_results, function(x){intersect(x$significant_uniq_data$overlap_gene, overlap_genes_depression) %>% length()/length(overlap_genes_depression)})

lapply(overlap_results, function(x){intersect(x$significant_uniq_data$overlap_gene, overlap_genes_mental_health) %>% length()/length(overlap_genes_mental_health)})


data.frame(
  psychosocial = sapply(overlap_results, function(x) {
    length(intersect(x$significant_uniq_data$overlap_gene, overlap_genes_psychosocial_factors)) / length(overlap_genes_psychosocial_factors)
  }),
  depression = sapply(overlap_results, function(x) {
    length(intersect(x$significant_uniq_data$overlap_gene, overlap_genes_depression)) / length(overlap_genes_depression)
  }),
  mental_health = sapply(overlap_results, function(x) {
    length(intersect(x$significant_uniq_data$overlap_gene, overlap_genes_mental_health)) / length(overlap_genes_mental_health)
  })
) -> results_df

results_df %>%
  arrange(desc(psychosocial), depression, mental_health)

