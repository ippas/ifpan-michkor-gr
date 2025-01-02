# Install yaml package if not already installed
if (!require("yaml")) install.packages("yaml")

# Load the yaml library
library(yaml)


# Directory path where the YAML files are stored
dir_path <- "data/prs-models-pan-biobank-uk/"

# Vector with file names
significant_phenotypes_prs <- c(
  "biobankuk-anti_gout_agent_microtuble_polymerization_inhibitor-both_sexes--na-EUR-1e-08.yml",
  "biobankuk-allopurinol-both_sexes--na-EUR-1e-08.yml",
  "biobankuk-30280-both_sexes--immature_reticulocyte_fraction-EUR-1e-08.yml",
  "biobankuk-30090-both_sexes--platelet_crit-EUR-1e-08.yml",
  "biobankuk-egfrcreacys-both_sexes--estimated_glomerular_filtration_rate_cystain_c-EUR-1e-08.yml",
  "biobankuk-20151-both_sexes--forced_vital_capacity_fvc_best_measure-EUR-1e-08.yml",
  "biobankuk-3063-both_sexes--forced_expiratory_volume_in_1_second_fev1_-EUR-1e-08.yml"
)

# Initialize a list to store the contents of the YAML files
yaml_contents <- list()

# Loop through the vector and read each YAML file
for (file_name in significant_phenotypes_prs) {
  file_path <- paste0(dir_path, file_name) # Create the full file path
  yaml_contents[[file_name]] <- yaml.load_file(file_path) # Read and store the file content
}


significant_models_variants_df <- map_df(significant_phenotypes_prs, function(model_name) {
  file_path <- paste0(dir_path, model_name)
  model_data <- yaml.load_file(file_path)
  variants_df <- model_data$score_model$variants %>% 
    lapply(as.data.frame) %>% 
    bind_rows(.id = "rsid")
  model_clean_name <- gsub("-1e-08.yml", "", model_name) # Remove the file extension for the model name
  variants_df <- variants_df %>% mutate(model = model_clean_name) # Add the model name as a column
  return(variants_df)
}, .id = "model_name") %>% 
  select(model, everything()) %>% 
  select(-model_name) %>% 
  separate(gnomadid, into = c("chromosome", "position_with_extra"), sep = ":") %>%
  mutate(position = str_split(position_with_extra, "_", simplify = TRUE)[,1]) %>%
  select(-position_with_extra) %>% 
  select(c(model, rsid, effect_allele, ref, effect_size, symbol, chromosome, position))


write.table(significant_models_variants_df, 
            "data/1kg-vcf/signif-sportsmen-phenotypes-rsid.tsv",
            row.names = F,
            col.names = T,
            quote = F,
            sep = "\t")


