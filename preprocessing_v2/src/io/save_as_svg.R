save_as_svg <- function(plot, filename, width = 6, height = 4) {
  if (inherits(plot, "ggplot")) {
    ggsave(filename = filename, plot = plot, device = "svg",
           width = width, height = height)
  } else if (inherits(plot, c("Heatmap", "HeatmapList"))) {
    grDevices::svg(filename = filename, width = width, height = height)
    draw(plot)
    dev.off()
  } else {
    warning("Nieznany typ obiektu – nie zapisano pliku.")
  }
  
  return(plot)
}
