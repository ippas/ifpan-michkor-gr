# Define the path to the folder containing the R scripts to be sourced
script_folder <- "preprocessing/preprocessing-genes-papers"

# List of script names to be sourced
script_names <- c("adrenal-gland.R", "arc.R", "asm.R", "embryos-scRNAseq.R", 
                  "hipocamp-fs.R", "install-load-packages.R", "lung.R", 
                  "nervous-cells.R", "npscs-ko-cav1.R", "npscs.R", 
                  "placenta-pnss.R", "prefrontal-cortex.R")

# Loop through each script name to construct the full path and source it
for (script in script_names) {
  script_path <- file.path(script_folder, script)
  source(script_path)
}

# Create an empty list called gene_paper_preprocessing_database
gene_paper_preprocessing_database <- list()

# Create a named vector that maps list names to existing DataFrame variables
name_mapping <- c("placenta-PNSS" = "placanta_to_database", 
                  "adrenal-gland" = "adrenal_gland_database", 
                  "NPSCs-KO-Cav1" = "npscs_ko_cav1_database", 
                  "embryos-scRNA-seq" = "embryos_scRNAseq_database", 
                  "hipocamp-fs" = "hippocampus_fs_database",
                  "prefrontal-cortex" = "prefrontal_cortex_database",
                  "ASM" = "asm_database",
                  "nervous-cells" = "nervous_cells_database",
                  "NPSCs" = "npscs_database",
                  "lung" = "lung_database",
                  "ARC" = "arc_database")  # Assuming ARC is already in your workspace


# Loop to populate gene_paper_preprocessing_database
for (new_name in names(name_mapping)) {
  old_name <- name_mapping[new_name]
  
  # Check if the DataFrame exists
  if (exists(old_name)) {
    gene_paper_preprocessing_database[[new_name]] <- get(old_name)
    
  } else {
    warning(paste("Variable", old_name, "not found. Skipping."))
  }
}


save_data_to_tsv(data_list = gene_paper_preprocessing_database, path = "data/supplement-genes-papers/preprocessing-genes-papers/")
