# Load necessary packages and functions for preprocessing
source("preprocessing/preprocessing-genes-papers/install-load-packages.R")

# Read the ARC gene data from a TSV file into a data frame named arc_df
read.table("data/supplement-genes-papers/ARC/arc-genes.tsv", header = TRUE) -> arc_df

# Create a new data frame (arc_database) that matches the desired template
# The data frame will include constants and values sourced from arc_df
data.frame(
  # Information about the source article
  article_source = "https://www.sciencedirect.com/science/article/pii/S2212877819302522",  # URL of the article
  pmid = 31176677,  # PubMed ID of the article
  
  # Gene-specific information
  gene_name = arc_df$gene_name,  # Extract gene names from arc_df
  ensembl_id = NA,  # Ensembl IDs (not available)
  log2ratio = NA,  # Log2 ratio (not available)
  regulation = arc_df$regulation,  # Extract regulation direction from arc_df
  
  # Experimental details
  species = "mouse",  # Species used in the experiment
  tissue = "ARC",  # Tissue where the experiment was performed
  cell = "glia|neurons",  # Cell types involved in the experiment
  environment = "in-vivo",  # Environmental context of the experiment
  treatment = "corticosterone",  # Treatment applied
  dose = "75μg/mL",  # Dose of the treatment
  time = arc_df$time,  # Extract time from arc_df
  
  # Analysis details
  fdr_threshold = 0.1,  # False Discovery Rate threshold
  method = "RNA-seq",  # Method used for gene expression analysis
  comparison = "CORT_vs_Ethanol",  # Conditions being compared
  strain = "C57BL/6",  # Strain of the animals used
  statistical_method = "DESeq2",  # Statistical method used for analysis
  treatment_type = "acute"  # Type of treatment (acute/chronic)
) -> arc_database  # Save the resulting data frame as arc_database
