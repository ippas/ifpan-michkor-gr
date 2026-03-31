create_customRect_patch_legend <- function(
    colors,
    labels,
    box_size = 1,
    spacing = 0.6,
    text_size = 5,
    box_linewidth = 1.2,
    fontface = "plain",
    ratio = 1,
    x_start = 0
) {
  n <- length(colors)
  if (length(labels) != n) stop("`colors` and `labels` must have the same length.")
  
  df <- data.frame(
    xmin = seq(x_start, by = (box_size + spacing), length.out = n),
    xmax = seq(x_start + box_size, by = (box_size + spacing), length.out = n),
    ymin = 0,
    ymax = 1,
    color = colors,
    label = labels
  )
  
  ggplot(df) +
    geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, color = color),
              fill = NA, linewidth = box_linewidth, show.legend = FALSE) +
    geom_text(aes(x = xmax + 0.3, y = (ymin + ymax) / 2, label = label),
              hjust = 0, size = text_size, fontface = fontface) +
    scale_color_identity() +
    coord_fixed(ratio = ratio) +
    xlim(min(df$xmin) - 0.5, max(df$xmax) + 2.5) +
    ylim(0, 1.5) +
    theme_void() +
    theme(plot.margin = margin(t = 5, b = 5))
}
