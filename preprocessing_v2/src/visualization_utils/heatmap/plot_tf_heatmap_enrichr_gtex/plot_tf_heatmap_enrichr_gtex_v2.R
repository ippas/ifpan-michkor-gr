plot_tf_heatmap_enrichr_gtex <- function(
    enrichr_list,                         # list(signature -> df), np. res$enrichr$overlap_only
    nuclearReceptors_GO0004879,           # character vector
    signature_order = c(
      "global_up", "global_down",
      "brain_up",  "brain_down",
      "blood_up",  "blood_down",
      "lung_up",   "lung_down"
    ),
    
    ## --- selection ---
    top_n = 10,
    
    ## --- TF extraction ---
    tf_regex = "^[^ ]+",
    
    ## --- alias mapping (optional) ---
    alias_map = c(
      "LXR"      = "NR1H3",
      "SA1"      = "STAG1",
      "SMC1"     = "SMC1A",
      "CJUN"     = "JUN",
      "TCFCP2L1" = "TFCP2L1",
      "Nerf2"    = "NFE2L2",
      "NRF2"     = "NFE2L2",
      "RING1B"   = "RING1",
      "AF4"      = "AFF4"
    ),
    
    ## --- FDR binning ---
    fdr_breaks = c(0, 1e-30, 1e-20, 1e-10, 1e-6, 1e-2, Inf),
    fdr_labels = c(
      "≤1e-30",
      "(1e-30, 1e-20]",
      "(1e-20, 1e-10]",
      "(1e-10, 1e-6]",
      "(1e-6, 0.01]",
      ">0.01"
    ),
    fill_pal = c(
      ">0.01"          = "white",
      "(1e-6, 0.01]"   = "#fae8e6",
      "(1e-10, 1e-6]"  = "#eba8a1",
      "(1e-20, 1e-10]" = "#d46b66",
      "(1e-30, 1e-20]" = "#b01728",
      "≤1e-30"         = "#5a0011"
    ),
    text_pal = c(
      ">0.01"          = "black",
      "(1e-6, 0.01]"   = "black",
      "(1e-10, 1e-6]"  = "black",
      "(1e-20, 1e-10]" = "black",
      "(1e-30, 1e-20]" = "white",
      "≤1e-30"         = "white"
    ),
    
    ## --- GTEx on/off ---
    use_gtex = TRUE,                     # <--- NOWE
    
    ## --- GTEx settings ---
    gtex_tissues_in  = c("Brain", "Lung", "Whole_Blood"),
    gtex_tissue_map  = c("Brain" = "brain", "Lung" = "lung", "Whole_Blood" = "blood"),
    gtex_threshold   = 1,                 # median_max > threshold
    gtex_fun         = multi_downloadMedianSubtissuesGTEx_tissueSummaryExpression,
    gtex_verbose     = TRUE,
    gtex_show_progress = TRUE,
    gtex_keep_unmapped = FALSE,
    
    ## --- dots ---
    dot_pal = c(brain = "#2E86FF", lung = "#00A676", blood = "#E31A1C"),
    dot_size = 3,
    dot_stroke = 0.8,
    dot_hollow_fill = "white",
    dot_hollow_color = "grey40",
    dot_offsets_global = c(brain = -0.23, lung = 0.00, blood = 0.23),
    dot_y_offset = -0.30,                 # relative to tile center
    
    ## --- tiles / spacing ---
    tile_width  = 0.97,
    tile_height = 0.97,
    tile_border_color = "white",
    tile_border_lwd = 0.6,
    
    ## --- nuclear receptor strip ---
    nr_purple = "#B388FF",
    strip_h = 0.10,                       # fraction of tile height (approx)
    strip_color = NULL,                   # if NULL uses nr_purple
    
    ## --- theme / labels ---
    base_size = 12,
    tf_text_size = 3,
    tf_text_face = "bold",
    x_text_angle = 45,
    x_text_hjust = 0,
    x_text_vjust = 0,
    y_text_face = "bold",
    panel_border_color = "black",
    panel_border_lwd = 1.2,
    legend_position = "right",
    y_label = "rank",
    
    ## --- SVG saving (optional) ---
    svg_file = NULL,
    svg_width = 10,
    svg_height = 8
) {
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(stringr)
    library(tibble)
    library(tidyr)
    library(ggplot2)
  })
  
  if (is.null(strip_color)) strip_color <- nr_purple
  
  ## ------------------------------------------------------------
  ## A) Top TFs per signature (top_n)
  ## ------------------------------------------------------------
  top10_per_signature <- enrichr_list %>%
    imap_dfr(~ .x %>%
               mutate(
                 TF  = str_extract(Term, tf_regex),
                 FDR = as.numeric(FDR)
               ) %>%
               group_by(TF) %>%
               slice_min(FDR, n = 1, with_ties = FALSE) %>%
               ungroup() %>%
               arrange(FDR) %>%
               slice_head(n = top_n) %>%
               mutate(rank = row_number(),
                      signature = .y) %>%
               select(signature, rank, TF, everything())
    )
  
  ## aliases (optional)
  if (!is.null(alias_map) && length(alias_map) > 0) {
    top10_per_signature <- top10_per_signature %>%
      mutate(TF = ifelse(TF %in% names(alias_map), unname(alias_map[TF]), TF))
  }
  
  ## ------------------------------------------------------------
  ## B) (OPTIONAL) Download GTEx summary for TFs
  ## ------------------------------------------------------------
  TF_expression_GTEx <- NULL
  TF_expr_flags      <- NULL
  
  if (isTRUE(use_gtex)) {
    TF_expression_GTEx <- gtex_fun(
      gene_symbol    = top10_per_signature$TF,
      verbose        = gtex_verbose,
      show_progress  = gtex_show_progress,
      keep_unmapped  = gtex_keep_unmapped
    )
    
    ## ------------------------------------------------------------
    ## D) GTEx flags per TF (brain/lung/blood)
    ##     expects: query_gene_symbol, tissue, median_max
    ## ------------------------------------------------------------
    TF_expr_flags <- TF_expression_GTEx %>%
      filter(tissue %in% gtex_tissues_in) %>%
      mutate(
        tissue = recode(tissue, !!!gtex_tissue_map),
        expressed = median_max > gtex_threshold
      ) %>%
      group_by(query_gene_symbol, tissue) %>%
      summarise(expressed = any(expressed), .groups = "drop") %>%
      pivot_wider(
        names_from  = tissue,
        values_from = expressed,
        values_fill = FALSE
      )
  }
  
  ## ------------------------------------------------------------
  ## C) Heatmap df (complete grid + bins + coords)
  ## ------------------------------------------------------------
  hm_df <- top10_per_signature %>%
    mutate(
      signature = factor(signature, levels = signature_order),
      rank      = factor(rank, levels = 1:top_n, ordered = TRUE)
    ) %>%
    complete(signature, rank) %>%
    mutate(
      signature = factor(signature, levels = signature_order),
      rank      = factor(rank, levels = 1:top_n, ordered = TRUE),
      FDR_num   = as.numeric(FDR),
      fdr_bin   = cut(
        FDR_num,
        breaks = fdr_breaks,
        labels = fdr_labels,
        include.lowest = TRUE
      ),
      fdr_bin   = factor(fdr_bin, levels = rev(fdr_labels)),
      text_col  = text_pal[as.character(fdr_bin)],
      is_NR     = !is.na(TF) & TF %in% nuclearReceptors_GO0004879,
      x_num     = as.integer(signature),
      y_num     = as.integer(rank)
    )
  
  ## ------------------------------------------------------------
  ## E) (OPTIONAL) Dots df (global = 3 dots; tissue-specific = 1 dot)
  ## ------------------------------------------------------------
  dots_df <- NULL
  
  if (isTRUE(use_gtex)) {
    dots_df <- hm_df %>%
      mutate(signature_chr = as.character(signature)) %>%
      left_join(TF_expr_flags, by = c("TF" = "query_gene_symbol")) %>%
      mutate(
        brain = ifelse(is.na(brain), FALSE, brain),
        lung  = ifelse(is.na(lung),  FALSE, lung),
        blood = ifelse(is.na(blood), FALSE, blood),
        sig_group = case_when(
          str_detect(signature_chr, "^global_") ~ "global",
          str_detect(signature_chr, "^brain_")  ~ "brain",
          str_detect(signature_chr, "^lung_")   ~ "lung",
          str_detect(signature_chr, "^blood_")  ~ "blood",
          TRUE ~ "global"
        )
      ) %>%
      pivot_longer(
        cols = c(brain, lung, blood),
        names_to = "tissue",
        values_to = "expressed"
      ) %>%
      mutate(
        tissue = factor(tissue, levels = c("brain", "lung", "blood")),
        x_dot = x_num + case_when(
          sig_group == "global" & tissue == "brain" ~ dot_offsets_global[["brain"]],
          sig_group == "global" & tissue == "lung"  ~ dot_offsets_global[["lung"]],
          sig_group == "global" & tissue == "blood" ~ dot_offsets_global[["blood"]],
          TRUE ~ 0
        ),
        y_dot = y_num + dot_y_offset
      ) %>%
      filter(sig_group == "global" | tissue == sig_group)
  }
  
  ## ------------------------------------------------------------
  ## F) Plot
  ## ------------------------------------------------------------
  p <- ggplot(hm_df, aes(signature, y_num, fill = fdr_bin)) +
    geom_tile(
      width = tile_width,
      height = tile_height,
      color = tile_border_color,
      size = tile_border_lwd,
      na.rm = FALSE
    ) +
    geom_rect(
      data = hm_df %>% filter(is_NR),
      inherit.aes = FALSE,
      aes(
        xmin = x_num - tile_width/2,
        xmax = x_num + tile_width/2,
        ymin = y_num + tile_height/2 - strip_h,
        ymax = y_num + tile_height/2
      ),
      fill = strip_color,
      color = NA
    ) +
    ## TF text
    geom_text(
      aes(label = TF),
      color = hm_df$text_col,
      size = tf_text_size,
      fontface = tf_text_face,
      na.rm = TRUE
    ) +
    scale_y_reverse(breaks = 1:top_n, labels = 1:top_n) +
    scale_x_discrete(position = "top") +
    scale_fill_manual(values = fill_pal, drop = FALSE, name = "FDR") +
    labs(x = NULL, y = y_label) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = element_text(angle = x_text_angle, hjust = x_text_hjust, vjust = x_text_vjust),
      axis.text.y = element_text(face = y_text_face),
      legend.position = legend_position,
      panel.border = element_rect(colour = panel_border_color, fill = NA, size = panel_border_lwd)
    )
  
  ## --- add dots + dot legend only if GTEx ON ---
  if (isTRUE(use_gtex)) {
    p <- p +
      ## hollow dots (always)
      geom_point(
        data = dots_df,
        inherit.aes = FALSE,
        aes(x_dot, y_dot),
        shape  = 21,
        size   = dot_size,
        stroke = dot_stroke,
        fill   = dot_hollow_fill,
        color  = dot_hollow_color,
        na.rm  = TRUE
      ) +
      ## expressed dots (colored overlay)
      geom_point(
        data = dots_df %>% filter(expressed),
        inherit.aes = FALSE,
        aes(x_dot, y_dot, color = tissue),
        shape  = 16,
        size   = dot_size,
        na.rm  = TRUE
      ) +
      scale_color_manual(
        values = dot_pal,
        breaks = c("brain", "lung", "blood"),
        labels = c("Brain", "Lung", "Blood"),
        name = "TF expression (GTEx)"
      )
  }
  
  ## ------------------------------------------------------------
  ## G) Optional SVG save
  ## ------------------------------------------------------------
  if (!is.null(svg_file)) {
    svg(filename = svg_file, width = svg_width, height = svg_height)
    print(p)
    dev.off()
  }
  
  ## return everything useful
  return(list(
    plot = p,
    hm_df = hm_df,
    dots_df = dots_df,
    TF_expr_flags = TF_expr_flags,
    top10_per_signature = top10_per_signature,
    TF_expression_GTEx = TF_expression_GTEx
  ))
}


out <- plot_tf_heatmap_enrichr_gtex(
  enrichr_list = res$enrichr$overlap_only,
  nuclearReceptors_GO0004879 = nuclearReceptors_GO0004879,
  use_gtex = FALSE,
  tile_width = 0.97,
  tile_height = 0.97
)
out$plot
