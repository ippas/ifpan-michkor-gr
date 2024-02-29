affymetrix_to_genomic_annotations <- function(affy_ids, ensembl_version = "default", species) {

  # Translates Affymetrix probe IDs to gene symbols and Ensembl IDs using BioMart.
  #
  # This function takes a vector of Affymetrix probe IDs, the desired Ensembl database version,
  # and the species (human, mouse, or rat) as inputs. It queries the appropriate BioMart
  # database to retrieve the corresponding gene symbols and Ensembl gene IDs.
  #
  # Args:
  #   affy_ids: A vector of Affymetrix probe IDs.
  #   ensembl_version: A character string indicating the Ensembl version to use. If "default",
  #                    the current default version of the Ensembl database is used. Otherwise,
  #                    this should specify the host name for a specific Ensembl archive version.
  #   species: A character string indicating the species. Valid options are "human", "mouse", or "rat".
  #
  # Returns:
  #   A dataframe with the Affymetrix probe IDs, corresponding gene symbols, and Ensembl gene IDs.
    
  # Map common species names to the appropriate BioMart dataset names
  species_dataset_name <- switch(species,
                                 "human" = "hsapiens_gene_ensembl",
                                 "mouse" = "mmusculus_gene_ensembl",
                                 "rat" = "rnorvegicus_gene_ensembl",
                                 stop("Unsupported species"))  # Stop execution if an unsupported species is provided
  
  # Debugging print statement (can be removed in production code)
  print(species_dataset_name)
  
  # Initialize the BioMart dataset based on the provided Ensembl version and species
  if (ensembl_version == "default") {
    mart <- useMart("ensembl", dataset = species_dataset_name)
  } else {
    # Use a specific Ensembl archive if a version other than "default" is specified
    mart <- useMart(biomart = "ensembl", dataset = species_dataset_name, host = ensembl_version)
  }
  
  # Conditional execution based on species to set up the correct attributes and filters for the BioMart query
  if (species == "human") {
    # Placeholder for human-specific logic; needs to be implemented based on the actual Affymetrix array used
  } else if (species == "mouse") {
    # Perform the BioMart query for mouse, using attributes and filters specific to the mouse Affymetrix array
    results <- getBM(attributes = c('affy_mouse430_2', 'external_gene_name', 'ensembl_gene_id'),
                     filters = 'affy_mouse430_2',
                     values = affy_ids,
                     mart = mart)
  } else if (species == "rat") {
    # Perform the BioMart query for rat, using attributes and filters specific to the rat Affymetrix array
    results <- getBM(attributes = c('affy_rat230_2', 'external_gene_name', 'ensembl_gene_id'),
                     filters = 'affy_rat230_2',
                     values = affy_ids,
                     mart = mart)
  }
  
  # Return the results of the BioMart query
  return(results)
}
