# ---------------------------------------------------------
# disgenet_gene2disease_chunks()
# ---------------------------------------------------------
# Input: character vector of gene symbols
# Output: list of results from disgenet2r::gene2disease()
# ---------------------------------------------------------

library(disgenet2r)
disgenet_gene2disease_chunks <- function(
    genes,
    chunk_size = 100,
    database = "CURATED",
    vocabulary = "HGNC",
    verbose = TRUE,
    columns_to_exclude = c("gene_symbol", 
                           "geneid", 
                           "ensemblid", 
                           "geneNcbiType", 
                           "geneDSI", 
                           "geneDPI",
                           "genepLI",
                           "uniprotids",
                           "protein_classid",
                           "protein_class_name",
                           "evidence_index", 
                           "yearInitial", 
                           "yearFinal", 
                           "score",
                           "numberPmidsWithChemsFiltered",
                           "numNCTSWithChemsIncludedInEvidences",
                           "evidence_level",
                           "chemsIncludedInEvidenceBySource", 
                           "numChemsIncludedInEvidences", 
                           "numPMIDSWithChemsIncludedInEvidences",
                           "numPMIDs",
                           "disease_prevalence_class",
                           "disease_prevalence_geo_area",
                           "disease_prevalence_type",
                           "disease_inheritance",
                           "numCTsupportingAssociation",
                           "numPMIDs",
                           "diseaseType"
    ),
    return_only_combine = TRUE,
    save_to_rds = NULL
) {
  # Pakiety wymagane
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("disgenet2r", quietly = TRUE)) stop("Package 'disgenet2r' is required.")
  if (!requireNamespace("progress", quietly = TRUE)) stop("Package 'progress' is required.")
  
  # Split genes into chunks
  gene_chunks <- split(genes, ceiling(seq_along(genes) / chunk_size))
  n_chunks <- length(gene_chunks)
  if (verbose) message("🧬 Total genes: ", length(genes), " → ", n_chunks, " chunks")
  
  results_list <- vector("list", n_chunks)
  
  # Progress bar setup
  pb <- progress::progress_bar$new(
    format = "⏳ Processing chunk :current/:total [:bar] :percent ETA: :eta",
    total = n_chunks,
    width = 60,
    clear = FALSE
  )
  
  # Iterate over chunks
  for (i in seq_along(gene_chunks)) {
    pb$tick()
    chunk_genes <- gene_chunks[[i]]
    if (verbose) message("\n[", i, "/", n_chunks, "] Processing ", length(chunk_genes), " genes")
    
    res <- tryCatch(
      disgenet2r::gene2disease(
        gene = chunk_genes,
        database = database,
        vocabulary = vocabulary,
        verbose = FALSE
      ),
      error = function(e) {
        if (verbose) message("⚠️ Error in chunk ", i, ": ", e$message)
        return(NULL)
      }
    )
    
    if (!is.null(res) && !is.null(res@qresult)) {
      df <- res@qresult %>%
        dplyr::select(!dplyr::any_of(columns_to_exclude)) %>%
        unique()
      results_list[[i]] <- df
    } else {
      results_list[[i]] <- NULL
      if (verbose) message("⚠️ Chunk ", i, " returned no valid results.")
    }
  }
  
  # Combine all results
  combine <- results_list %>%
    purrr::compact() %>%
    dplyr::bind_rows() %>%
    unique()
  
  # Optional save
  if (!is.null(save_to_rds)) {
    if (verbose) message("💾 Saving combined results to RDS: ", save_to_rds)
    tryCatch({
      saveRDS(combine, file = save_to_rds)
      if (verbose) message("✅ File saved successfully.")
    }, error = function(e) {
      message("❌ Failed to save file: ", e$message)
    })
  }
  
  # Return results
  if (return_only_combine) {
    if (verbose) message("\n✅ Returning combined results only.")
    return(combine)
  } else {
    if (verbose) message("\n✅ Returning both chunks and combined results.")
    return(list(
      chunks = results_list,
      combine = combine
    ))
  }
}


# =========================================================
# 💻 Example usage
# =========================================================
# disgenet_metadataDiseases <- disgenet_gene2disease_chunks(
#   genes = hgnc_symbols_vector_v110,
#   chunk_size = 100,
#   verbose = F,
#   return_only_combine = T,
#   save_to_rds = "data/databases/disgenet/disgenet_metadataDiseases.rds"
# )




