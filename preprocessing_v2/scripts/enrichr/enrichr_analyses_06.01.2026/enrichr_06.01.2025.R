sig_names <- c(
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp",
  "minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_LungCellsDown",
  "minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown",
  "global_GR_genes_globalUp5TissuesDerivedCells",
  "global_GR_genes_globalDown5TissuesDerivedCells"
)

flat_allGrSignatures_31.10.2025[sig_names] %>% 
  lapply(., length)



names_mapper <- c(
  minusGlobalUpDown5TissuesDerivedCells_BloodCellsUp    = "blood_up",
  minusGlobalUpDown5TissuesDerivedCells_BloodCellsDown  = "blood_down",
  minusGlobalUpDown5TissuesDerivedCells_LungCellsUp     = "lung_up",
  minusGlobalUpDown5TissuesDerivedCells_LungCellsDown   = "lung_down",
  minusGlobalUpDown5TissuesDerivedCells_NeuralCellsUp   = "brain_up",
  minusGlobalUpDown5TissuesDerivedCells_NeuralCellsDown = "brain_down",
  global_GR_genes_globalUp5TissuesDerivedCells          = "global_up",
  global_GR_genes_globalDown5TissuesDerivedCells        = "global_down"
)

run_enrichr_multi <- function(
    gene_lists,
    database = "ChEA_2022",
    min_overlap_genes = 3,
    fdr_threshold = 0.05,
    names_mapper = NULL,
    order_signatures = c(
      "global_up", "global_down",
      "brain_up", "brain_down",
      "blood_up", "blood_down",
      "lung_up", "lung_down"
    ),
    drop_old_pvals = TRUE,
    verbose = TRUE,
    xlsx_file = NULL   # default: no saving; if provided, write XLSX
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(tibble)
    library(stringr)
  })
  
  if (!is.list(gene_lists) || length(gene_lists) == 0) {
    stop("gene_lists must be a non-empty list.")
  }
  
  # ============================================================
  # 0. APPLY NAMES MAPPER (GLOBAL, ONCE)
  # ============================================================
  if (!is.null(names_mapper)) {
    
    if (is.null(names(gene_lists))) {
      stop("gene_lists must be a named list when using names_mapper.")
    }
    
    missing_map <- setdiff(names(gene_lists), names(names_mapper))
    if (length(missing_map) > 0) {
      stop(
        "names_mapper missing entries for: ",
        paste(missing_map, collapse = ", ")
      )
    }
    
    names(gene_lists) <- names_mapper[names(gene_lists)]
  }
  
  # ============================================================
  # ORDER INPUT LISTS (AFFECTS OUTPUT ORDER + SHEET ORDER)
  # ============================================================
  if (!is.null(order_signatures)) {
    present <- intersect(order_signatures, names(gene_lists))
    rest    <- setdiff(names(gene_lists), present)
    gene_lists <- gene_lists[c(present, rest)]
  }
  
  # ============================================================
  # HELPERS
  # ============================================================
  count_unique_associated_genes <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(0L)
    if (!("Genes" %in% colnames(df))) return(0L)
    
    df$Genes %>%
      as.character() %>%
      strsplit(";", fixed = TRUE) %>%
      unlist() %>%
      str_trim() %>%
      (\(x) x[nzchar(x)])() %>%
      unique() %>%
      length() %>%
      as.integer()
  }
  
  standardize_enrichr <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    
    # rename to pvalue / FDR
    if (all(c("P.value", "Adjusted.P.value") %in% colnames(df))) {
      df <- df %>%
        dplyr::rename(
          pvalue = P.value,
          FDR    = Adjusted.P.value
        )
    }
    
    # drop old pvals
    if (drop_old_pvals) {
      drop_cols <- intersect(
        c("Old.P.value", "Old.Adjusted.P.value"),
        colnames(df)
      )
      if (length(drop_cols) > 0) {
        df <- df %>% dplyr::select(-all_of(drop_cols))
      }
    }
    
    df
  }
  
  # ============================================================
  # RUN ENRICHR
  # ============================================================
  results <- imap(gene_lists, function(gv, nm) {
    
    gv <- unique(as.character(gv))
    gv <- gv[!is.na(gv) & nzchar(gv)]
    
    if (verbose) message("Enrichr: ", nm, " (n_genes=", length(gv), ")")
    
    raw <- if (length(gv) == 0) tibble() else run_enrichr(gv, database)
    raw <- standardize_enrichr(raw)
    
    # export results: overlap filter ONLY (no FDR filtering)
    overlap_only <- if ("n_genes" %in% colnames(raw)) {
      raw %>% dplyr::filter(.data$n_genes >= min_overlap_genes)
    } else {
      tibble()
    }
    
    # compute filtered for summary + coloring
    filtered <- if ("FDR" %in% colnames(overlap_only)) {
      overlap_only %>% dplyr::filter(.data$FDR <= fdr_threshold)
    } else {
      tibble()
    }
    
    list(
      input_genes_n = length(gv),
      raw = raw,
      overlap_only = overlap_only,  # <-- written to XLSX
      filtered = filtered
    )
  })
  
  # ============================================================
  # SUMMARY TABLE (ORDERED)
  # ============================================================
  summary_df <- imap_dfr(results, function(res, nm) {
    
    n_genes_list  <- res$input_genes_n
    n_assoc       <- if (nrow(res$filtered) == 0) 0L else nrow(res$filtered)
    n_assoc_genes <- count_unique_associated_genes(res$filtered)
    prop_assoc    <- if (n_genes_list > 0) n_assoc_genes / n_genes_list else NA_real_
    
    tibble(
      signature = nm,
      n_genes_in_list = n_genes_list,
      n_associations = n_assoc,
      n_associated_genes = n_assoc_genes,
      prop_associated_genes = prop_assoc
    )
  })
  
  if (!is.null(order_signatures)) {
    summary_df <- summary_df %>%
      mutate(signature = factor(signature, levels = c(order_signatures, setdiff(unique(signature), order_signatures)))) %>%
      arrange(signature) %>%
      mutate(signature = as.character(signature))
  }
  
  # ============================================================
  # GLOBAL SUMMARY (ALL SIGNATURES)
  # - counts after overlap + FDR
  # - top rows:
  #   * min pvalue across ALL signatures
  #   * max n_genes across ALL signatures
  #   preference: filtered_all; fallback: overlap_only_all
  # ============================================================
  overlap_all <- bind_rows(map(results, "overlap_only"), .id = "signature")
  filtered_all <- bind_rows(map(results, "filtered"), .id = "signature")
  
  n_associations_all <- if (nrow(filtered_all) == 0) 0L else as.integer(nrow(filtered_all))
  
  n_unique_genes_all <- {
    if (nrow(filtered_all) == 0 || !("Genes" %in% colnames(filtered_all))) {
      0L
    } else {
      filtered_all$Genes %>%
        as.character() %>%
        strsplit(";", fixed = TRUE) %>%
        unlist() %>%
        str_trim() %>%
        (\(x) x[nzchar(x)])() %>%
        unique() %>%
        length() %>%
        as.integer()
    }
  }
  
  top_source <- if (nrow(filtered_all) > 0) filtered_all else overlap_all
  
  global_min_pvalue_row <- {
    if (nrow(top_source) == 0 || !("pvalue" %in% colnames(top_source))) {
      tibble()
    } else {
      top_source %>%
        mutate(pvalue = suppressWarnings(as.numeric(.data$pvalue))) %>%
        filter(!is.na(.data$pvalue)) %>%
        arrange(.data$pvalue) %>%
        slice(1) %>%
        as_tibble()
    }
  }
  
  global_max_ngenes_row <- {
    if (nrow(top_source) == 0 || !("n_genes" %in% colnames(top_source))) {
      tibble()
    } else {
      top_source %>%
        mutate(n_genes = suppressWarnings(as.numeric(.data$n_genes))) %>%
        filter(!is.na(.data$n_genes)) %>%
        arrange(desc(.data$n_genes)) %>%
        slice(1) %>%
        as_tibble()
    }
  }
  
  out <- list(
    summary = summary_df,
    global_summary = tibble(
      n_associations_all_signatures = n_associations_all,
      n_unique_associated_genes_all_signatures = n_unique_genes_all
    ),
    global_top = list(
      min_pvalue_row = global_min_pvalue_row,   # full Enrichr row + signature
      max_n_genes_row = global_max_ngenes_row   # full Enrichr row + signature
    ),
    enrichr = list(
      raw = map(results, "raw"),
      overlap_only = map(results, "overlap_only"),
      filtered = map(results, "filtered")
    )
  )
  
  # ============================================================
  # OPTIONAL: WRITE XLSX
  # - write overlap_only results (overlap>=min_overlap_genes), NO FDR filtering
  # - one sheet per signature; sheet order = order_signatures / input order
  # - freeze first row + first column
  # - first column auto-width (Term)
  # - format:
  #   Odds.Ratio, Combined.Score -> 3 decimals
  #   pvalue, FDR -> scientific, 3 decimals
  # - highlight rows passing FDR threshold in light green
  # ============================================================
  if (!is.null(xlsx_file)) {
    
    suppressPackageStartupMessages({
      library(openxlsx)
    })
    
    wb <- openxlsx::createWorkbook()
    
    sheet_names <- names(gene_lists)
    
    for (sn in sheet_names) {
      
      df <- out$enrichr$overlap_only[[sn]]
      if (is.null(df)) df <- tibble()
      
      # Excel sheet name constraints
      sheet_safe <- sn
      sheet_safe <- gsub("[:\\\\/\\?\\*\\[\\]]", "_", sheet_safe)
      if (nchar(sheet_safe) > 31) sheet_safe <- substr(sheet_safe, 1, 31)
      
      openxlsx::addWorksheet(wb, sheet_safe)
      
      # write
      openxlsx::writeData(
        wb,
        sheet = sheet_safe,
        x = df,
        withFilter = TRUE
      )
      
      # freeze first row AND first column
      openxlsx::freezePane(
        wb,
        sheet = sheet_safe,
        firstRow = TRUE,
        firstCol = TRUE
      )
      
      # auto-width for first column (Term)
      openxlsx::setColWidths(
        wb,
        sheet = sheet_safe,
        cols = 1,
        widths = "auto"
      )
      
      # nothing else to style if empty
      if (nrow(df) == 0 || ncol(df) == 0) next
      
      coln <- colnames(df)
      idx_or  <- match("Odds.Ratio", coln)
      idx_cs  <- match("Combined.Score", coln)
      idx_p   <- match("pvalue", coln)
      idx_fdr <- match("FDR", coln)
      
      nR <- nrow(df)
      nC <- ncol(df)
      data_rows <- 2:(nR + 1)  # header row is 1
      
      style_dec3 <- openxlsx::createStyle(numFmt = "0.000")
      style_sci3 <- openxlsx::createStyle(numFmt = "0.000E+00")
      
      if (!is.na(idx_or)) {
        openxlsx::addStyle(
          wb, sheet_safe, style = style_dec3,
          rows = data_rows, cols = idx_or,
          gridExpand = TRUE, stack = TRUE
        )
      }
      
      if (!is.na(idx_cs)) {
        openxlsx::addStyle(
          wb, sheet_safe, style = style_dec3,
          rows = data_rows, cols = idx_cs,
          gridExpand = TRUE, stack = TRUE
        )
      }
      
      if (!is.na(idx_p)) {
        openxlsx::addStyle(
          wb, sheet = sheet_safe, style = style_sci3,
          rows = data_rows, cols = idx_p,
          gridExpand = TRUE, stack = TRUE
        )
      }
      
      if (!is.na(idx_fdr)) {
        openxlsx::addStyle(
          wb, sheet = sheet_safe, style = style_sci3,
          rows = data_rows, cols = idx_fdr,
          gridExpand = TRUE, stack = TRUE
        )
      }
      
      # conditional formatting: highlight rows with FDR <= threshold (light green)
      if (!is.na(idx_fdr) && nR > 0) {
        fdr_col_letter <- openxlsx::int2col(idx_fdr)
        rule <- paste0("$", fdr_col_letter, "2<=", fdr_threshold)
        
        style_green <- openxlsx::createStyle(
          fgFill = "#C6EFCE"
        )
        
        openxlsx::conditionalFormatting(
          wb, sheet = sheet_safe,
          cols = 1:nC, rows = data_rows,
          type = "expression",
          rule = rule,
          style = style_green
        )
      }
    }
    
    openxlsx::saveWorkbook(wb, file = xlsx_file, overwrite = TRUE)
    out$xlsx_file <- xlsx_file
  }
  
  return(out)
}

# ============================================================
# EXAMPLE RUN (no save)
# ============================================================
# res <- run_enrichr_multi(
#   gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
#   database = "ChEA_2022",
#   min_overlap_genes = 3,
#   fdr_threshold = 0.05,
#   names_mapper = names_mapper
# )

# ============================================================
# EXAMPLE RUN (save)
# ============================================================
res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "ChEA_2022",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_CHEA2022_GRsignatures_06.01.2026.xlsx"
)

res$global_top


res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "GO_Biological_Process_2025",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_GOBiologicalProcess_GRsignatures_06.01.2026.xlsx"
)

res$global_top

res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "COMPARTMENTS_Experimental_2025",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_COMPARTMENTSExperimental2025_GRsignatures_06.01.2026.xlsx"
)

res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "KEGG_2021_Human",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_KEGG2021Human_GRsignatures_06.01.2026.xlsx"
)


res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "DrugMatrix",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_DrugMatrix_GRsignatures_07.01.2026.xlsx"
)


res$enrichr$overlap_only %>% 
  lapply(., head, 10)


res <- run_enrichr_multi(
  gene_lists = flat_allGrSignatures_31.10.2025[sig_names],
  database = "Sciplex_Drug_Perturbation_Signatures_2025",
  min_overlap_genes = 3,
  fdr_threshold = 0.05,
  names_mapper = names_mapper,
  xlsx_file = "results_v2/overlap/enrichr/enrichr_SciplexDrugPerturbationSignatures2025_GRsignatures_07.01.2026.xlsx"
)
