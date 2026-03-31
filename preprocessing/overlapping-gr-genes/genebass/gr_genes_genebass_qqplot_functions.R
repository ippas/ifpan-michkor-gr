# ##############################################################################
# ---- prepare data ----
# ##############################################################################
lite_grSignatures 

genebass_mentalHealth_skat <-read.delim("data/genebass/mentalHealth_Pvalue_SKAT_0.05.tsv.bgz")
genebass_mentalHealth_skat_all <-read.delim("data/genebass/genebass_mental_health_all_SKAT.tsv.bgz")

genebass_mentalHealth_skat %>% head

# ##############################################################################
# ---- prepare function ----
# ##############################################################################
qqplot_genebass_by_phenocode <- function(
    genebass_data,
    phenocode,
    pvalue_column = "pvalue",
    annotation_filter = NULL,
    title_prefix = "QQ-plot:",
    highlight_genes = NULL,
    highlight_color = "blue"
) {
  # Filtrowanie danych
  df <- genebass_data %>%
    filter(phenocode == !!phenocode) %>%
    filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  # Filtrowanie po typach mutacji
  if (!is.null(annotation_filter)) {
    df <- df %>% filter(annotation %in% annotation_filter)
  }
  
  # Sprawdzenie dostępności danych
  if (nrow(df) == 0) {
    stop(paste("Brak danych dla phenocode:", phenocode))
  }
  
  # Metadane
  description <- df$description[1] %||% paste("Phenocode", phenocode)
  n_variants <- nrow(df)
  n_genes <- n_distinct(df$gene_symbol)
  annotation_types <- df$annotation %>% unique() %>% sort() %>% paste(collapse = ", ")
  
  # Przygotowanie danych do wykresu
  pvalues <- df[[pvalue_column]]
  order_idx <- order(pvalues)
  
  qq_df <- data.frame(
    expected = -log10(ppoints(length(pvalues))),
    observed = -log10(sort(pvalues)),
    mutation_type = df$annotation[order_idx],
    gene_symbol = df$gene_symbol[order_idx],
    highlight = FALSE
  )
  
  if (!is.null(highlight_genes)) {
    qq_df$highlight <- qq_df$gene_symbol %in% highlight_genes
  }
  
  # Wykres
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(data = subset(qq_df, !highlight), size = 1) +
    geom_point(data = subset(qq_df, highlight), color = highlight_color, size = 1.5) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = paste0(title_prefix, " ", description, " (ID: ", phenocode, ")"),
      subtitle = paste0(
        "n genes = ", n_genes,
        ", n variants = ", n_variants,
        ", ", annotation_types
      )
    ) +
    theme_minimal()
}
qqplot_genebass_by_phenocode <- function(
    genebass_data,
    phenocode,
    pvalue_column = "pvalue",
    annotation_filter = NULL,
    title_prefix = "QQ-plot:",
    highlight_genes = NULL,
    highlight_color = "blue",
    label_box_padding = 0.7,
    label_point_padding = 0.5,
    title_text_size = 14,
    subtitle_text_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    label_text_size = 3
) {
  # Wymagany pakiet
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Pakiet 'ggrepel' jest wymagany. Zainstaluj go przez install.packages('ggrepel')")
  }
  
  # Filtrowanie danych
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  if (!is.null(annotation_filter)) {
    df <- df %>% dplyr::filter(annotation %in% annotation_filter)
  }
  
  if (nrow(df) == 0) {
    stop(paste("Brak danych dla phenocode:", phenocode))
  }
  
  # Metadane
  description <- df$description[1] %||% paste("Phenocode", phenocode)
  n_variants <- nrow(df)
  n_genes <- dplyr::n_distinct(df$gene_symbol)
  annotation_types <- df$annotation %>% unique() %>% sort() %>% paste(collapse = ", ")
  
  # QQ-plot data
  pvalues <- df[[pvalue_column]]
  order_idx <- order(pvalues)
  
  qq_df <- data.frame(
    expected = -log10(ppoints(length(pvalues))),
    observed = -log10(sort(pvalues)),
    mutation_type = df$annotation[order_idx],
    gene_symbol = df$gene_symbol[order_idx],
    beta = df$beta[order_idx],
    highlight = FALSE
  )
  
  if (!is.null(highlight_genes)) {
    qq_df$highlight <- qq_df$gene_symbol %in% highlight_genes
  }
  
  label_df <- qq_df %>%
    dplyr::filter(highlight) %>%
    dplyr::mutate(label = paste0(gene_symbol, " (β = ", round(beta, 5), ")"))
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(data = subset(qq_df, !highlight), size = 1) +
    geom_point(data = subset(qq_df, highlight), color = highlight_color, size = 1.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      color = highlight_color,
      size = label_text_size,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = label_box_padding,
      point.padding = label_point_padding,
      segment.color = "grey40",
      segment.linetype = "dashed"
    ) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = paste0(title_prefix, " ", description, " (ID: ", phenocode, ")"),
      subtitle = paste0(
        "n genes = ", n_genes,
        ", n variants = ", n_variants,
        ", annotation: ", annotation_types
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_text_size, face = "bold"),
      plot.subtitle = element_text(size = subtitle_text_size),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )
}

qqplot_genebass_by_phenocode <- function(
    genebass_data,
    phenocode,
    pvalue_column = "pvalue",
    annotation_filter = NULL,
    title_prefix = "QQ-plot:",
    highlight_genes = NULL,
    highlight_color = "blue",
    label_box_padding = 0.7,
    label_point_padding = 0.5,
    title_text_size = 14,
    subtitle_text_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    label_text_size = 3,
    label_fields = c("gene_symbol", "beta", "pvalue", "annotation")
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Pakiet 'ggrepel' jest wymagany. Zainstaluj go przez install.packages('ggrepel')")
  }
  
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  if (!is.null(annotation_filter)) {
    df <- df %>% dplyr::filter(annotation %in% annotation_filter)
  }
  
  if (nrow(df) == 0) {
    stop(paste("Brak danych dla phenocode:", phenocode))
  }
  
  description <- df$description[1] %||% paste("Phenocode", phenocode)
  n_variants <- nrow(df)
  n_genes <- dplyr::n_distinct(df$gene_symbol)
  annotation_types <- df$annotation %>% unique() %>% sort() %>% paste(collapse = ", ")
  
  qq_df <- df %>%
    dplyr::arrange(.data[[pvalue_column]]) %>%
    dplyr::mutate(
      expected = -log10(ppoints(dplyr::n())),
      observed = -log10(.data[[pvalue_column]]),
      highlight = if (!is.null(highlight_genes)) gene_symbol %in% highlight_genes else FALSE
    )
  
  label_df <- qq_df %>%
    dplyr::filter(highlight) %>%
    dplyr::mutate(
      pvalue_fmt = format(.data[[pvalue_column]], scientific = TRUE, digits = 2),
      label_line1 = paste0(
        if ("gene_symbol" %in% label_fields) gene_symbol else NULL,
        if ("beta" %in% label_fields) paste0(" (β = ", round(beta, 5), ")") else ""
      ),
      label_line2 = paste0(
        if ("pvalue" %in% label_fields) paste0("p = ", pvalue_fmt) else "",
        if ("annotation" %in% label_fields) paste0(", ", annotation) else ""
      ),
      label = paste(label_line1, label_line2, sep = "\n")
    )
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(data = subset(qq_df, !highlight), size = 1) +
    geom_point(data = subset(qq_df, highlight), color = highlight_color, size = 1.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      color = highlight_color,
      size = label_text_size,
      max.overlaps = Inf,
      box.padding = label_box_padding,
      point.padding = label_point_padding,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.linetype = "dashed"
    ) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = paste0(title_prefix, " ", description, " (ID: ", phenocode, ")"),
      subtitle = paste0(
        "n genes = ", n_genes,
        ", n variants = ", n_variants,
        ", annotation: ", annotation_types
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_text_size, face = "bold"),
      plot.subtitle = element_text(size = subtitle_text_size),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )
}


qqplot_genebass_by_phenocode <- function(
    genebass_data,
    phenocode,
    pvalue_column = "pvalue",
    annotation_filter = NULL,
    title_prefix = "QQ-plot:",
    highlight_genes = NULL,
    highlight_color = "blue",
    label_box_padding = 0.7,
    label_point_padding = 0.5,
    title_text_size = 14,
    subtitle_text_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    label_text_size = 3,
    label_fields = c("gene_symbol", "beta", "pvalue", "annotation"),
    label_strategy = "all"  # nowy argument: "all", "lowest_p", "highest_beta"
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Pakiet 'ggrepel' jest wymagany. Zainstaluj go przez install.packages('ggrepel')")
  }
  
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  if (!is.null(annotation_filter)) {
    df <- df %>% dplyr::filter(annotation %in% annotation_filter)
  }
  
  if (nrow(df) == 0) {
    stop(paste("Brak danych dla phenocode:", phenocode))
  }
  
  description <- df$description[1] %||% paste("Phenocode", phenocode)
  n_variants <- nrow(df)
  n_genes <- dplyr::n_distinct(df$gene_symbol)
  annotation_types <- df$annotation %>% unique() %>% sort() %>% paste(collapse = ", ")
  
  qq_df <- df %>%
    dplyr::arrange(.data[[pvalue_column]]) %>%
    dplyr::mutate(
      expected = -log10(ppoints(dplyr::n())),
      observed = -log10(.data[[pvalue_column]]),
      highlight = if (!is.null(highlight_genes)) gene_symbol %in% highlight_genes else FALSE
    )
  
  if (!is.null(highlight_genes)) {
    label_df <- switch(
      label_strategy,
      "all" = qq_df %>% dplyr::filter(highlight),
      "lowest_p" = qq_df %>%
        dplyr::filter(highlight) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_min(.data[[pvalue_column]], with_ties = FALSE) %>%
        dplyr::ungroup(),
      "highest_beta" = qq_df %>%
        dplyr::filter(highlight) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_max(abs(beta), with_ties = FALSE) %>%
        dplyr::ungroup(),
      stop("Niepoprawna wartość argumentu 'label_strategy'. Użyj: 'all', 'lowest_p' lub 'highest_beta'.")
    )
    
    label_df <- label_df %>%
      dplyr::mutate(
        pvalue_fmt = format(.data[[pvalue_column]], scientific = TRUE, digits = 2),
        label_line1 = paste0(
          if ("gene_symbol" %in% label_fields) gene_symbol else NULL,
          if ("beta" %in% label_fields) paste0(" (β = ", round(beta, 5), ")") else ""
        ),
        label_line2 = paste0(
          if ("pvalue" %in% label_fields) paste0("p = ", pvalue_fmt) else "",
          if ("annotation" %in% label_fields) paste0(", ", annotation) else ""
        ),
        label = paste(label_line1, label_line2, sep = "\n")
      )
  } else {
    label_df <- NULL
  }
  
  ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(data = subset(qq_df, !highlight), size = 1) +
    geom_point(data = subset(qq_df, highlight), color = highlight_color, size = 1.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      color = highlight_color,
      size = label_text_size,
      max.overlaps = Inf,
      box.padding = label_box_padding,
      point.padding = label_point_padding,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.linetype = "dashed"
    ) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = paste0(title_prefix, " ", description, " (ID: ", phenocode, ")"),
      subtitle = paste0(
        "n genes = ", n_genes,
        ", n variants = ", n_variants,
        ", annotation: ", annotation_types
      )
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_text_size, face = "bold"),
      plot.subtitle = element_text(size = subtitle_text_size),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )
}


qqplot_genebass_by_phenocode_v2 <- function(
    genebass_data,
    phenocode,
    pvalue_column = "pvalue",
    annotation_filter = NULL,
    title_prefix = "QQ-plot:",
    highlight_genes = NULL,
    highlight_color = "blue",
    label_box_padding = 0.7,
    label_point_padding = 0.5,
    title_text_size = 14,
    subtitle_text_size = 12,
    axis_title_size = 12,
    axis_text_size = 10,
    label_text_size = 3,
    label_fields = c("gene_symbol", "beta", "pvalue", "annotation"),
    label_strategy = "all",
    lambda_gc = FALSE,
    lambda_gc_text_size = 4,
    gene_signature_name = NULL  # Nowy argument
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Pakiet 'ggrepel' jest wymagany. Zainstaluj go przez install.packages('ggrepel')")
  }
  
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  if (!is.null(annotation_filter)) {
    df <- df %>% dplyr::filter(annotation %in% annotation_filter)
  }
  
  if (nrow(df) == 0) {
    stop(paste("Brak danych dla phenocode:", phenocode))
  }
  
  description <- df$description[1] %||% paste("Phenocode", phenocode)
  n_variants <- nrow(df)
  n_genes <- dplyr::n_distinct(df$gene_symbol)
  annotation_types <- df$annotation %>% unique() %>% sort() %>% paste(collapse = ", ")
  
  lambda_value <- if (lambda_gc) {
    chisq_vals <- qchisq(1 - df[[pvalue_column]], df = 1)
    median(chisq_vals, na.rm = TRUE) / 0.455
  } else {
    NA
  }
  
  qq_df <- df %>%
    dplyr::arrange(.data[[pvalue_column]]) %>%
    dplyr::mutate(
      expected = -log10(ppoints(dplyr::n())),
      observed = -log10(.data[[pvalue_column]]),
      highlight = if (!is.null(highlight_genes)) gene_symbol %in% highlight_genes else FALSE
    )
  
  if (!is.null(highlight_genes)) {
    label_df <- switch(
      label_strategy,
      "all" = qq_df %>% dplyr::filter(highlight),
      "lowest_p" = qq_df %>%
        dplyr::filter(highlight) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_min(.data[[pvalue_column]], with_ties = FALSE) %>%
        dplyr::ungroup(),
      "highest_beta" = qq_df %>%
        dplyr::filter(highlight) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_max(abs(beta), with_ties = FALSE) %>%
        dplyr::ungroup(),
      stop("Niepoprawna wartość argumentu 'label_strategy'. Użyj: 'all', 'lowest_p' lub 'highest_beta'.")
    )
    
    label_df <- label_df %>%
      dplyr::mutate(
        pvalue_fmt = format(.data[[pvalue_column]], scientific = TRUE, digits = 2),
        label_line1 = paste0(
          if ("gene_symbol" %in% label_fields) gene_symbol else NULL,
          if ("beta" %in% label_fields) paste0(" (β = ", round(beta, 5), ")") else ""
        ),
        label_line2 = paste0(
          if ("pvalue" %in% label_fields) paste0("p = ", pvalue_fmt) else "",
          if ("annotation" %in% label_fields) paste0(", ", annotation) else ""
        ),
        label = paste(label_line1, label_line2, sep = "\n")
      )
  } else {
    label_df <- NULL
  }
  
  subtitle_text <- paste0(
    "n genes = ", n_genes,
    ", n variants = ", n_variants,
    ", annotation: ", annotation_types,
    if (!is.null(gene_signature_name)) paste0(", gene signature: ", gene_signature_name) else ""
  )
  
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(data = subset(qq_df, !highlight), size = 1) +
    geom_point(data = subset(qq_df, highlight), color = highlight_color, size = 1.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      color = highlight_color,
      size = label_text_size,
      max.overlaps = Inf,
      box.padding = label_box_padding,
      point.padding = label_point_padding,
      min.segment.length = 0,
      segment.color = "grey40",
      segment.linetype = "dashed"
    ) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = paste0(title_prefix, " ", description, " (ID: ", phenocode, ")"),
      subtitle = subtitle_text
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_text_size, face = "bold"),
      plot.subtitle = element_text(size = subtitle_text_size),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )
  
  if (lambda_gc) {
    p <- p + annotate("text", x = 0.1, y = max(qq_df$observed, na.rm = TRUE),
                      label = paste0("λGC = ", round(lambda_value, 3)),
                      hjust = 0, vjust = 1, size = lambda_gc_text_size, fontface = "plain")
  }
  
  return(p)
}

qqplot_genebass_by_phenocode_v3 <- function(
    genebass_data,
    phenocode,
    pvalue_column = "pvalue",
    annotation_filter = NULL,
    title_prefix = "QQ-plot:",
    highlight_genes = NULL,
    highlight_color = "blue",
    label_box_padding = 0.7,
    label_point_padding = 0.5,
    label_text_size = 3,
    label_force = 2,
    label_force_pull = 0.2,
    label_max_iter = 10000,
    label_min_segment_length = 0.1,
    lambda_gc = FALSE,
    lambda_gc_text_size = 4,
    gene_signature_name = NULL,
    label_strategy = "all",
    label_fields = c("gene_symbol", "beta", "pvalue", "annotation"),
    title_text_size = 14,
    subtitle_text_size = 12,
    axis_title_size = 12,
    axis_text_size = 10
) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Pakiet 'ggrepel' jest wymagany. Zainstaluj go przez install.packages('ggrepel')")
  }
  
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  if (!is.null(annotation_filter)) {
    df <- df %>% dplyr::filter(annotation %in% annotation_filter)
  }
  
  if (nrow(df) == 0) {
    stop(paste("Brak danych dla phenocode:", phenocode))
  }
  
  description <- df$description[1] %||% paste("Phenocode", phenocode)
  n_variants <- nrow(df)
  n_genes <- dplyr::n_distinct(df$gene_symbol)
  annotation_types <- df$annotation %>% unique() %>% sort() %>% paste(collapse = ", ")
  
  lambda_value <- if (lambda_gc) {
    chisq_vals <- qchisq(1 - df[[pvalue_column]], df = 1)
    median(chisq_vals, na.rm = TRUE) / 0.455
  } else {
    NA
  }
  
  qq_df <- df %>%
    dplyr::arrange(.data[[pvalue_column]]) %>%
    dplyr::mutate(
      expected = -log10(ppoints(dplyr::n())),
      observed = -log10(.data[[pvalue_column]]),
      highlight = if (!is.null(highlight_genes)) gene_symbol %in% highlight_genes else FALSE
    )
  
  if (!is.null(highlight_genes)) {
    label_df <- switch(
      label_strategy,
      "all" = qq_df %>% dplyr::filter(highlight),
      "lowest_p" = qq_df %>%
        dplyr::filter(highlight) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_min(.data[[pvalue_column]], with_ties = FALSE) %>%
        dplyr::ungroup(),
      "highest_beta" = qq_df %>%
        dplyr::filter(highlight) %>%
        dplyr::group_by(gene_symbol) %>%
        dplyr::slice_max(abs(beta), with_ties = FALSE) %>%
        dplyr::ungroup(),
      stop("Niepoprawna wartość argumentu 'label_strategy'. Użyj: 'all', 'lowest_p' lub 'highest_beta'.")
    )
    
    label_df <- label_df %>%
      dplyr::mutate(
        pvalue_fmt = format(.data[[pvalue_column]], scientific = TRUE, digits = 2),
        label_line1 = paste0(
          if ("gene_symbol" %in% label_fields) gene_symbol else NULL,
          if ("beta" %in% label_fields) paste0(" (β = ", round(beta, 5), ")") else ""
        ),
        label_line2 = paste0(
          if ("pvalue" %in% label_fields) paste0("p = ", pvalue_fmt) else "",
          if ("annotation" %in% label_fields) paste0(", ", annotation) else ""
        ),
        label = paste(label_line1, label_line2, sep = "\n")
      )
  } else {
    label_df <- NULL
  }
  
  subtitle_text <- paste0(
    "n genes = ", n_genes,
    ", n variants = ", n_variants,
    ", annotation: ", annotation_types,
    if (!is.null(gene_signature_name)) paste0(", gene signature: ", gene_signature_name) else ""
  )
  
  p <- ggplot(qq_df, aes(x = expected, y = observed)) +
    geom_point(data = subset(qq_df, !highlight), size = 1) +
    geom_point(data = subset(qq_df, highlight), color = highlight_color, size = 1.5) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = label),
      color = highlight_color,
      size = label_text_size,
      max.overlaps = Inf,
      box.padding = label_box_padding,
      point.padding = label_point_padding,
      force = label_force,
      force_pull = label_force_pull,
      max.iter = label_max_iter,
      min.segment.length = label_min_segment_length,
      segment.color = "grey40",
      segment.linetype = "dashed"
    ) +
    labs(
      x = "Expected -log10(p)",
      y = "Observed -log10(p)",
      title = paste0(title_prefix, " ", description, " (ID: ", phenocode, ")"),
      subtitle = subtitle_text
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = title_text_size, face = "bold"),
      plot.subtitle = element_text(size = subtitle_text_size),
      axis.title = element_text(size = axis_title_size),
      axis.text = element_text(size = axis_text_size)
    )
  
  if (lambda_gc) {
    p <- p + annotate("text", x = 0.1, y = max(qq_df$observed, na.rm = TRUE),
                      label = paste0("λGC = ", round(lambda_value, 3)),
                      hjust = 0, vjust = 1, size = lambda_gc_text_size, fontface = "plain")
  }
  
  return(p)
}


calculate_lambda_gc <- function(genebass_data, phenocode, pvalue_column = "pvalue") {
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0)
  
  if (nrow(df) == 0) {
    warning(paste("Brak prawidłowych p-value dla phenocode:", phenocode))
    return(NA)
  }
  
  chisq_vals <- qchisq(1 - df[[pvalue_column]], df = 1)
  lambda_gc <- median(chisq_vals, na.rm = TRUE) / 0.455
  
  return(lambda_gc)
}

calculate_lambda_gc <- function(genebass_data, phenocode, pvalue_column = "pvalue") {
  df <- genebass_data %>%
    dplyr::filter(phenocode == !!phenocode) %>%
    dplyr::filter(!is.na(.data[[pvalue_column]]) & .data[[pvalue_column]] > 0 & .data[[pvalue_column]] < 1)
  
  if (nrow(df) == 0) {
    warning(paste("Brak prawidłowych p-value dla phenocode:", phenocode))
    return(NA)
  }
  
  chisq_vals <- qchisq(1 - df[[pvalue_column]], df = 1)
  lambda_gc <- median(chisq_vals, na.rm = TRUE) / 0.455
  
  return(lambda_gc)
}

# ##############################################################################
# ---- analysis ----
# ##############################################################################

 