library(tibble)

gr_signature_long <- tribble(
  ~grSignature_derivation, ~tissue, ~direction, ~value,
  
  # systemic
  "systemic", "brain", "UP",   87,
  "systemic", "brain", "DOWN", 38,
  "systemic", "brain", "SUM", 125,
  
  "systemic", "blood", "UP",   66,
  "systemic", "blood", "DOWN", 10,
  "systemic", "blood", "SUM",  76,
  
  "systemic", "lung", "UP",    35,
  "systemic", "lung", "DOWN",  35,
  "systemic", "lung", "SUM",   70,
  
  # brainCells
  "brainCells", "brain", "UP",   63,
  "brainCells", "brain", "DOWN", 22,
  "brainCells", "brain", "SUM",  85,
  
  "brainCells", "blood", "UP",   0,
  "brainCells", "blood", "DOWN", 10,
  "brainCells", "blood", "SUM",  10,
  
  "brainCells", "lung", "UP",    0,
  "brainCells", "lung", "DOWN",  3,
  "brainCells", "lung", "SUM",   3,
  
  # bloodCells
  "bloodCells", "brain", "UP",   16,
  "bloodCells", "brain", "DOWN", 14,
  "bloodCells", "brain", "SUM",  30,
  
  "bloodCells", "blood", "UP",   10,
  "bloodCells", "blood", "DOWN", 24,
  "bloodCells", "blood", "SUM",  34,
  
  "bloodCells", "lung", "UP",    0,
  "bloodCells", "lung", "DOWN",  9,
  "bloodCells", "lung", "SUM",   9,
  
  # lungCells
  "lungCells", "brain", "UP",   14,
  "lungCells", "brain", "DOWN", 28,
  "lungCells", "brain", "SUM",  42,
  
  "lungCells", "blood", "UP",   9,
  "lungCells", "blood", "DOWN", 9,
  "lungCells", "blood", "SUM",  18,
  
  "lungCells", "lung", "UP",    43,
  "lungCells", "lung", "DOWN",  12,
  "lungCells", "lung", "SUM",   55
)

gr_signature_long


library(dplyr)
library(ggplot2)


# Dane do słupków
bars_df <- gr_signature_long %>%
  filter(direction %in% c("UP", "DOWN")) %>%
  mutate(
    tissue = factor(tissue, levels = c("brain", "blood", "lung")),
    grSignature_derivation = factor(
      grSignature_derivation,
      levels = c("systemic", "brainCells", "bloodCells", "lungCells")
    ),
    direction = factor(direction, levels = c("UP", "DOWN"))
  )

# Dane do etykiet SUM
sum_df <- gr_signature_long %>%
  filter(direction == "SUM") %>%
  mutate(
    tissue = factor(tissue, levels = c("brain", "blood", "lung")),
    grSignature_derivation = factor(
      grSignature_derivation,
      levels = c("systemic", "brainCells", "bloodCells", "lungCells")
    )
  )
# kolejność stackowania: UP na dole, DOWN na górze
bars_df <- bars_df %>%
  mutate(direction = factor(direction, levels = c("DOWN", "UP")))

svg(
  filename = "/home/mateusz/projects/ifpan-michkor-gr/results_v2/enrichr/enrichr_cellMarker_18.01.2026/barplot_enrichrCellMarkerGrouped_grSignatures_18.01.2026.svg",
  width = 10,   # mm -> cale (180 mm, publikacyjna szerokość)
  height = 6  # mm -> cale
  # pointsize = 12
)

ggplot(bars_df, aes(x = grSignature_derivation, y = value, fill = direction)) +
  # geom_col(width = 0.55) +
  geom_col(width = 0.55, colour = "black", linewidth = 2) +

  facet_wrap(~ tissue, nrow = 1) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(values = c(UP = "firebrick", DOWN = "navyblue")) +
  labs(
    x = "GR signature derivation",
    y = "Number of genes",
    fill = "Regulation"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    
    strip.background = element_rect(fill = "white", colour = "black", size = 1),
    strip.text = element_text(face = "bold", size = 14, colour = "black"),
    
    axis.title.x = element_text(size = 14, colour = "black"),
    axis.title.y = element_text(size = 14, colour = "black"),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1, colour = "black"),
    axis.text.y  = element_text(size = 14, colour = "black"),
    
    legend.title = element_text(size = 14, colour = "black"),
    legend.text  = element_text(size = 14, colour = "black"),
    
    panel.border = element_rect(colour = "black", fill = NA, size = 1.4),
    panel.grid = element_blank()
  ) +
  geom_text(
    data = sum_df,
    aes(x = grSignature_derivation, y = value, label = value),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 14/2.83,
    fontface = "bold"
  ) 

dev.off()

