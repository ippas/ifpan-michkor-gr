multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression <- function(
    gene_symbol,
    verbose = TRUE,
    keep_unmapped = FALSE,
    show_progress = TRUE,
    quiet_single = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
  })
  
  # ============================================================
  # Validate input
  # ============================================================
  if (is.null(gene_symbol) || length(gene_symbol) == 0) {
    stop("gene_symbol must be a non-empty character vector.")
  }
  
  gene_symbol <- as.character(gene_symbol)
  gene_symbol <- unique(gene_symbol)
  gene_symbol <- gene_symbol[!is.na(gene_symbol) & gene_symbol != ""]
  
  if (show_progress && verbose) {
    message(sprintf(
      "GTEx: processing %d gene(s)...",
      length(gene_symbol)
    ))
  }
  
  # ============================================================
  # Iterate over genes
  # ============================================================
  res_list <- purrr::map(gene_symbol, function(gs) {
    
    if (show_progress && verbose) {
      message(" - Processing gene: ", gs)
    }
    
    # Safe execution: failure of one gene does not stop the whole run
    out <- tryCatch(
      {
        downloadMedianSubtissuesGTEx_tissueSummaryExpression(
          gene_symbol   = gs,
          verbose       = if (quiet_single) FALSE else verbose,
          keep_unmapped = keep_unmapped
        )
      },
      error = function(e) {
        # Return a single-row placeholder with error information
        tibble(
          ontologyId = NA_character_,
          datasetId  = NA_character_,
          gencodeId  = NA_character_,
          geneSymbol = NA_character_,
          unit       = NA_character_,
          tissue     = NA_character_,
          n_subtissues = NA_integer_,
          median_mean  = NA_real_,
          median_sd    = NA_real_,
          median_min   = NA_real_,
          median_max   = NA_real_,
          query_gene_symbol = gs,
          error_msg = conditionMessage(e)
        )
      }
    )
    
    # Add information about the queried gene symbol
    if (!("query_gene_symbol" %in% colnames(out))) {
      out <- out %>%
        mutate(query_gene_symbol = gs, .before = 1)
    }
    
    out
  })
  
  # ============================================================
  # Combine results
  # ============================================================
  result <- bind_rows(res_list)
  
  if (verbose) {
    n_ok <- sum(is.na(result$error_msg))
    n_err <- sum(!is.na(result$error_msg))
    message(sprintf(
      "GTEx: completed. Successful genes: %d | Failed genes: %d",
      n_ok, n_err
    ))
  }
  
  return(result)
}


multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression <- function(
    gene_symbol,
    verbose = TRUE,
    keep_unmapped = FALSE,
    show_progress = TRUE,
    quiet_single = TRUE
) {
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
  })

  # ============================================================
  # Validate input
  # ============================================================
  if (is.null(gene_symbol) || length(gene_symbol) == 0) {
    stop("gene_symbol must be a non-empty character vector.")
  }

  gene_symbol <- as.character(gene_symbol)
  gene_symbol <- unique(gene_symbol)
  gene_symbol <- gene_symbol[!is.na(gene_symbol) & gene_symbol != ""]

  n_total <- length(gene_symbol)
  i <- 0

  if (show_progress && verbose) {
    message(sprintf("GTEx: processing %d gene(s)...", n_total))
  }

  # ============================================================
  # Iterate over genes
  # ============================================================
  res_list <- purrr::map(gene_symbol, function(gs) {

    # Progress counter
    i <<- i + 1
    if (show_progress && verbose) {
      message(sprintf(
        "GTEx: %d/%d (%.1f%%) — %s",
        i, n_total, 100 * i / n_total, gs
      ))
    }

    # Safe execution: failure of one gene does not stop the whole run
    out <- tryCatch(
      {
        downloadMedianSubtissuesGTEx_tissueSummaryExpression(
          gene_symbol   = gs,
          verbose       = if (quiet_single) FALSE else verbose,
          keep_unmapped = keep_unmapped
        )
      },
      error = function(e) {
        tibble(
          ontologyId = NA_character_,
          datasetId  = NA_character_,
          gencodeId  = NA_character_,
          geneSymbol = NA_character_,
          unit       = NA_character_,
          tissue     = NA_character_,
          n_subtissues = NA_integer_,
          median_mean  = NA_real_,
          median_sd    = NA_real_,
          median_min   = NA_real_,
          median_max   = NA_real_,
          query_gene_symbol = gs,
          error_msg = conditionMessage(e)
        )
      }
    )

    # Add information about the queried gene symbol
    if (!("query_gene_symbol" %in% colnames(out))) {
      out <- out %>%
        mutate(query_gene_symbol = gs, .before = 1)
    }

    # Ensure error_msg exists for successful calls
    if (!("error_msg" %in% colnames(out))) {
      out <- out %>%
        mutate(error_msg = NA_character_)
    }

    out
  })

  # ============================================================
  # Combine results
  # ============================================================
  result <- bind_rows(res_list)

  if (verbose) {
    n_ok  <- sum(is.na(result$error_msg))
    n_err <- sum(!is.na(result$error_msg))
    message(sprintf(
      "GTEx: completed. Successful genes: %d | Failed genes: %d",
      n_ok, n_err
    ))
  }

  return(result)
}


