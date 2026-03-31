downloadMedianSubtissuesGTEx_tissueSummaryExpression <- function(
    gene_symbol,
    verbose = TRUE,
    keep_unmapped = FALSE
) {
  suppressPackageStartupMessages({
    library(dplyr)
  })
  
  # ============================================================
  # INTERNAL MAPPER (subtissue -> tissue)
  # ============================================================
  tissue_detailed_mapper <- list(
    ## BRAIN
    Brain = c(
      "Brain_Amygdala",
      "Brain_Anterior_cingulate_cortex_BA24",
      "Brain_Caudate_basal_ganglia",
      "Brain_Cerebellar_Hemisphere",
      "Brain_Cerebellum",
      "Brain_Cortex",
      "Brain_Frontal_Cortex_BA9",
      "Brain_Hippocampus",
      "Brain_Hypothalamus",
      "Brain_Nucleus_accumbens_basal_ganglia",
      "Brain_Putamen_basal_ganglia",
      "Brain_Spinal_cord_cervical_c-1",
      "Brain_Substantia_nigra"
    ),
    
    ## Adipose
    Adipose_Tissue = c(
      "Adipose_Subcutaneous",
      "Adipose_Visceral_Omentum"
    ),
    
    Blood_Vessel = c(
      "Artery_Aorta",
      "Artery_Coronary",
      "Artery_Tibial",
      "Cells_EBV-transformed_lymphocytes"
    ),
    
    Cervix_Uteri = c(
      "Cervix_Ectocervix",
      "Cervix_Endocervix"
    ),
    
    Colon = c(
      "Colon_Sigmoid",
      "Colon_Transverse"
    ),
    
    Esophagus = c(
      "Esophagus_Gastroesophageal_Junction",
      "Esophagus_Mucosa",
      "Esophagus_Muscularis"
    ),
    
    Heart = c(
      "Heart_Atrial_Appendage",
      "Heart_Left_Ventricle"
    ),
    
    Kidney = c(
      "Kidney_Cortex",
      "Kidney_Medulla"
    ),
    
    Skin = c(
      "Skin_Not_Sun_Exposed_Suprapubic",
      "Skin_Sun_Exposed_Lower_leg"
    ),
    
    Adrenal_Gland        = c("Adrenal_Gland"),
    Bladder              = c("Bladder"),
    Fallopian_Tube       = c("Fallopian_Tube"),
    Liver                = c("Liver"),
    Lung                 = c("Lung"),
    Minor_Salivary_Gland = c("Minor_Salivary_Gland"),
    Muscle_Skeletal      = c("Muscle_Skeletal"),
    Nerve_Tibial         = c("Nerve_Tibial"),
    Ovary                = c("Ovary"),
    Pancreas             = c("Pancreas"),
    Pituitary            = c("Pituitary"),
    Prostate             = c("Prostate"),
    Spleen               = c("Spleen"),
    Stomach              = c("Stomach"),
    Testis               = c("Testis"),
    Thyroid              = c("Thyroid"),
    Uterus               = c("Uterus"),
    Vagina               = c("Vagina"),
    Whole_Blood          = c("Whole_Blood")
  )
  
  tissue_mapper_vec <- unlist(
    lapply(names(tissue_detailed_mapper), function(tissue) {
      setNames(
        rep(tissue, length(tissue_detailed_mapper[[tissue]])),
        tissue_detailed_mapper[[tissue]]
      )
    }),
    use.names = TRUE
  )
  
  # ============================================================
  # 1) gene_symbol -> gencodeId
  # ============================================================
  g <- gtexr::get_gene_search(gene_symbol)
  if (is.null(g) || nrow(g) == 0) {
    stop(sprintf("GTEx: nie znaleziono genu '%s'.", gene_symbol))
  }
  if (!("gencodeId" %in% colnames(g))) {
    stop("GTEx: get_gene_search() nie zwróciło kolumny 'gencodeId'.")
  }
  
  idx <- 1L
  if ("geneSymbol" %in% colnames(g)) {
    hit <- which(toupper(g$geneSymbol) == toupper(gene_symbol))
    if (length(hit) > 0) idx <- hit[1]
  }
  
  gencode_id <- g$gencodeId[idx]
  if (verbose) message(sprintf("GTEx: %s -> %s", gene_symbol, gencode_id))
  
  # ============================================================
  # 2) medians per subtissue
  # ============================================================
  med <- gtexr::get_median_gene_expression(gencodeIds = gencode_id)
  if (is.null(med) || nrow(med) == 0) {
    stop("GTEx: get_median_gene_expression() zwróciło pusty wynik.")
  }
  if (!("tissueSiteDetailId" %in% colnames(med))) {
    stop("GTEx: brak kolumny 'tissueSiteDetailId' w medianach.")
  }
  
  median_col <- intersect(c("median", "medianExpression", "median_expression"), colnames(med))
  if (length(median_col) == 0) {
    stop(paste0(
      "GTEx: nie znalazłem kolumny z medianą. Kolumny: ",
      paste(colnames(med), collapse = ", ")
    ))
  }
  median_col <- median_col[1]
  
  med2 <- med %>%
    mutate(
      tissue = unname(tissue_mapper_vec[tissueSiteDetailId])
    )
  
  if (!keep_unmapped) {
    med2 <- med2 %>% filter(!is.na(tissue))
  }
  
  # ============================================================
  # 3) Tissue summary (required meta + 5 summary columns)
  # ============================================================
  required_meta <- c("ontologyId", "datasetId", "gencodeId", "geneSymbol", "unit")
  missing_meta <- setdiff(required_meta, colnames(med2))
  if (length(missing_meta) > 0) {
    stop(sprintf(
      "GTEx: brakuje wymaganych kolumn meta w medianach: %s",
      paste(missing_meta, collapse = ", ")
    ))
  }
  
  meta_row <- med2 %>%
    select(all_of(required_meta)) %>%
    distinct() %>%
    slice(1)
  
  tissue_summary <- med2 %>%
    group_by(tissue) %>%
    summarise(
      n_subtissues = dplyr::n(),
      median_mean  = mean(.data[[median_col]], na.rm = TRUE),
      median_sd    = sd(.data[[median_col]], na.rm = TRUE),
      median_min   = min(.data[[median_col]], na.rm = TRUE),
      median_max   = max(.data[[median_col]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(median_mean)) %>%
    mutate(
      ontologyId = meta_row$ontologyId,
      datasetId  = meta_row$datasetId,
      gencodeId  = meta_row$gencodeId,
      geneSymbol = meta_row$geneSymbol,
      unit       = meta_row$unit
    ) %>%
    select(
      ontologyId, datasetId, gencodeId, geneSymbol, unit,
      tissue,
      n_subtissues, median_mean, median_sd, median_min, median_max
    )
  
  return(tissue_summary)
}


# tmp_gtex <- downloadMedianSubtissuesGTEx_tissueSummaryExpression(
#     gene_symbol = c("NR3C1", "AR"),
#     verbose = TRUE,
#     keep_unmapped = FALSE
# ) 
