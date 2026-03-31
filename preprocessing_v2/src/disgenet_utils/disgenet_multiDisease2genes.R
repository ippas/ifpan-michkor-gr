disgenet_multiDisease2genes <- function(
    disease_ids,
    database = "CURATED",
    score = c(0, 1),
    filter_genes = NULL,
    verbose = TRUE
) {
  
  # -------------------------------------------
  # NORMALIZATION OF DISEASE IDs
  # -------------------------------------------
  disease_ids <- ifelse(
    grepl("^UMLS_", disease_ids),
    disease_ids,
    paste0("UMLS_", disease_ids)
  )
  
  if (verbose)
    message("🧬 Running multiDisease2genes on ", length(disease_ids), " diseases...")
  
  # start time
  start_time <- Sys.time()
  
  # -------------------------------------------
  # PREPARE OUTPUT LIST
  # -------------------------------------------
  results_list <- list()
  
  # Columns we keep
  columns_keep <- c(
    "gene_symbol", "geneid", "ensemblid", "uniprotids",
    "protein_classid", "protein_class_name",
    "diseaseVocabularies", "disease_name", "diseaseUMLSCUI",
    "diseaseClasses_MSH", "diseaseClasses_UMLS_ST",
    "diseaseClasses_DO", "diseaseClasses_HPO",
    "score", "evidence_index"
  )
  
  # time tracking
  iter_times <- c()
  last_iter_time <- Sys.time()
  
  # -------------------------------------------
  # MAIN LOOP
  # -------------------------------------------
  for (i in seq_along(disease_ids)) {
    
    d_id <- disease_ids[i]
    if (verbose) message("\n[", i, "/", length(disease_ids), "] ", d_id)
    
    # --- ETA + TOTAL ELAPSED ---
    if (i > 1) {
      iter_times <- c(iter_times, as.numeric(difftime(Sys.time(), last_iter_time, units = "secs")))
      mean_t <- mean(iter_times)
      sd_t   <- ifelse(length(iter_times) > 1, sd(iter_times), 0)
      rem    <- (length(disease_ids) - i + 1) * mean_t
      
      total_elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      
      if (verbose) {
        message(sprintf(
          "⏳ mean: %.2fs ± %.2f | ETA: %ds | total elapsed: %ds",
          mean_t, sd_t, round(rem), round(total_elapsed)
        ))
      }
    }
    last_iter_time <- Sys.time()
    
    # -------------------------------------------
    # QUERY DISGENET
    # -------------------------------------------
    res <- tryCatch(
      disgenet2r::disease2gene(
        disease  = d_id,
        database = database,
        score    = score
      ),
      error = function(e) {
        if (verbose) message("⚠️ Error for ", d_id, ": ", e$message)
        return(NULL)
      }
    )
    
    # Case 1: "no results for the query"
    if (is.character(res)) {
      if (verbose) message("⚠️ No results for ", d_id)
      next
    }
    
    # Case 2: invalid object
    if (!methods::is(res, "DataGeNET.DGN")) {
      if (verbose) message("⚠️ Invalid object for ", d_id)
      next
    }
    
    # Case 3: empty qresult
    if (nrow(res@qresult) == 0) {
      if (verbose) message("⚠️ Empty qresult for ", d_id)
      next
    }
    
    # -------------------------------------------
    # DATA + FILTERS
    # -------------------------------------------
    df <- as.data.frame(res@qresult)
    
    if (!is.null(filter_genes)) {
      df <- dplyr::filter(df, gene_symbol %in% filter_genes)
      if (nrow(df) == 0) {
        if (verbose) message("⚠️ All genes filtered out for ", d_id)
        next
      }
    }
    
    curated_name <- gsub(" ", "_", df$disease_name[1])
    
    df <- df %>%
      dplyr::select(any_of(columns_keep)) %>%
      dplyr::mutate(
        curated_disease_name = curated_name,
        input_disease_id = d_id
      )
    
    results_list[[curated_name]] <- df
  }
  
  return(results_list)
}



# ##############################################################################
# ---- example ----
# ##############################################################################
# ids <- c("C1832916", "C0005586", "C0036341", "C0525045")
# 
# res <- disgenet_multiDisease2genes(
#   disease_ids  = ids,
#   filter_genes = hgnc_symbols_vector_v110,
#   verbose = T
# )
# 
# 
# res$Bipolar_Disorder
