plotly_overlap_heatmap_interactive <- function(
    data_list,
    data_type,
    p_threshold = 0.05,
    overlap_threshold = 1,
    color_scale_range = NULL,
    palette = c("navy", "white", "firebrick3")
) {
  # 📦 Pakiety
  library(dplyr)
  library(tidyr)
  library(plotly)
  library(stringr)
  
  # 1️⃣ Pobranie danych
  data <- data_list[[data_type]]$list
  log2or <- data$log2_odds_ratio_matrix
  pval <- data$p_value_matrix
  chi2 <- data$chi2_matrix
  overlap <- data$number_overlap_matrix
  overlap_genes <- data$overlap_genes_matrix
  
  rows <- rownames(log2or)
  cols <- colnames(log2or)
  
  # 2️⃣ Przygotowanie tabeli do plotly
  df <- expand.grid(
    row = rows,
    col = cols,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      log2or = as.vector(log2or),
      pval = as.vector(pval),
      chi2 = as.vector(chi2),
      overlap = as.vector(overlap),
      genes = as.vector(overlap_genes)
    )
  
  # 3️⃣ Filtrowanie i tooltip
  df <- df %>%
    mutate(
      signif = ifelse(pval < p_threshold & overlap >= overlap_threshold, TRUE, FALSE),
      tooltip = paste0(
        "<b>", row, "</b> × <b>", col, "</b><br>",
        "log2(OR): ", sprintf("%.2f", log2or), "<br>",
        "p-value: ", signif(pval, 3), "<br>",
        "Chi²: ", sprintf("%.2f", chi2), "<br>",
        "Overlapping genes (n=", overlap, "):<br>",
        str_trunc(genes, 300, ellipsis = "…")
      )
    )
  
  # 4️⃣ Kolorystyka
  if (is.null(color_scale_range)) {
    max_abs <- max(abs(df$log2or), na.rm = TRUE)
    color_scale_range <- c(-max_abs, max_abs)
  }
  
  # 5️⃣ Interaktywna heatmapa
  p <- plot_ly(
    data = df,
    x = ~col,
    y = ~row,
    z = ~log2or,
    text = ~tooltip,
    hoverinfo = "text",
    type = "heatmap",
    colorscale = list(
      list(0, palette[1]),
      list(0.5, palette[2]),
      list(1, palette[3])
    ),
    zmin = color_scale_range[1],
    zmax = color_scale_range[2],
    colorbar = list(title = "log2(OR)")
  )
  
  # 6️⃣ Oznaczenie istotnych pól (ramka)
  signif_df <- df %>% filter(signif)
  if (nrow(signif_df) > 0) {
    p <- p %>% add_trace(
      data = signif_df,
      x = ~col,
      y = ~row,
      type = "scatter",
      mode = "markers",
      marker = list(
        symbol = "square-open",
        color = "black",
        size = 14,
        line = list(width = 2)
      ),
      hoverinfo = "none",
      showlegend = FALSE
    )
  }
  
  # 7️⃣ Ustawienia wyglądu
  p <- p %>%
    layout(
      title = paste0("Interactive overlap heatmap: ", data_type),
      xaxis = list(title = "", tickangle = -45),
      yaxis = list(title = "", autorange = "reversed"),
      margin = list(l = 100, b = 100)
    )
  
  return(p)
}
