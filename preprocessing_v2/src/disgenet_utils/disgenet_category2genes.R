# ============================================================
# disgenet_category2genes()
# ------------------------------------------------------------
# High-level wrapper for extracting DisGeNET gene associations
# for diseases belonging to one or more category terms.
#
# Features:
#   • Filters diseases based on user-defined category column
#   • Retrieves disease-associated genes from DisGeNET
#   • Returns:
#         - Raw gene lists (score 0–1)
#         - Gene lists filtered by multiple minimum score thresholds
#   • Supports optional gene filtering (e.g., only HGNC-validated genes)
#   • Optional saving of the complete output object to an RDS file
#   • Compatible with disgenet_multiDisease2genes()
#
# Arguments:
#   msh_categories   Character vector of category terms to match.
#   metadata_df      Data frame containing disease metadata
#                    (e.g. disgenet_metadataDiseases).
#   category_column  Column inside metadata_df used to match categories.
#                    Must be a list-column (e.g. diseaseClasses_MSH).
#   database         DisGeNET database source ("CURATED", "ALL", etc.).
#   score            Score range passed to DisGeNET API.
#   filter_genes     Optional vector of gene symbols to keep.
#   score_thresholds Numeric vector of minimum scores to generate
#                    filtered gene lists (e.g. c(0, 0.5, 0.7, 0.9)).
#   save_to_rds      Optional file path. If not NULL, the full output
#                    object is saved as an .rds file.
#   verbose          Logical. If TRUE, print progress messages.
#
# Returns:
#   A named list with the structure:
#       $raw_disgenet
#       $geneLists_scoreMin0
#       $geneLists_scoreMin0.5
#       $geneLists_scoreMin0.6
#       $geneLists_scoreMin0.7
#       $geneLists_scoreMin0.8
#       $geneLists_scoreMin0.9
#
# Each element is itself a list of data frames, named by curated_disease_name.
# ============================================================

disgenet_category2genes <- function(
    msh_categories,
    metadata_df,
    category_column = diseaseClasses_MSH,
    database = "CURATED",
    score = c(0, 1),
    filter_genes = NULL,
    score_thresholds = c(0, 0.5, 0.6, 0.7, 0.8, 0.9),
    subset_n = NULL,
    save_to_rds = NULL,
    verbose = TRUE
) {
  
  # -------------------------------------------
  # 1) Select diseases based on category column
  # -------------------------------------------
  if (verbose) {
    message("🔍 Selecting disease IDs for category/ies: ",
            paste(msh_categories, collapse = ", "))
  }
  
  col <- rlang::enquo(category_column)
  
  subset_df <- metadata_df %>%
    dplyr::filter(
      purrr::map_lgl(!!col, ~ any(.x %in% msh_categories))
    )
  
  if (nrow(subset_df) == 0) {
    stop("❌ No diseases matched provided categories in column: ",
         rlang::quo_text(col))
  }
  
  disease_ids <- unique(subset_df$diseaseid)
  
  if (!is.null(subset_n)) {
    if (verbose) message("🧪 TEST MODE: using only first ", subset_n, " diseases")
    disease_ids <- head(disease_ids, subset_n)
  }
  
  if (verbose) {
    message("📌 Number of disease IDs to process: ", length(disease_ids))
  }
  
  
  # -------------------------------------------
  # 2) Fetch raw DisGeNET results (score 0–1)
  # -------------------------------------------
  if (verbose) message("\n🎯 Fetching RAW DisGeNET results...")
  
  raw_list <- disgenet_multiDisease2genes(
    disease_ids  = disease_ids,
    database     = database,
    score        = score,
    filter_genes = filter_genes,
    verbose      = verbose
  )
  
  
  # -------------------------------------------
  # 3) Construct output structure
  # -------------------------------------------
  out <- list()
  out$raw_disgenet <- raw_list
  
  
  # -------------------------------------------
  # 4) Generate gene-symbol vectors for score thresholds
  # -------------------------------------------
  for (th in score_thresholds) {
    
    if (verbose)
      message(sprintf("\n🎯 Generating unique gene lists for score ≥ %.1f ...", th))
    
    filtered_list <- purrr::map(raw_list, function(df) {
      unique(df$gene_symbol[df$score >= th])
    })
    
    out[[paste0("geneLists_scoreMin", th)]] <- filtered_list
  }
  
  
  # -------------------------------------------
  # 5) Save output to RDS
  # -------------------------------------------
  if (!is.null(save_to_rds)) {
    if (verbose) message("💾 Saving output to RDS: ", save_to_rds)
    saveRDS(out, file = save_to_rds)
  }
  
  return(out)
}



# ##############################################################################
# ---- example ----
# ##############################################################################
# disgenet_mentalDisordersF03 <- disgenet_category2genes(
#   msh_categories = "Mental Disorders (F03)",
#   metadata_df    = disgenet_metadataDiseases,
#   filter_genes   = hgnc_symbols_vector_v110,
#   subset_n = 5,
#   save_to_rds    = "data/databases/disgenet/disgenet_mentalDisordersF03.rds"
# )
