# ##############################################################################
# ---- Neural ----
# ##############################################################################
df_icd10 <- tibble(
  icd10_category = c("F00-09","F10-19","F20-29","F30-39","F40-49",
                     "F50-59","F60-69","F70-79","F80-89","F90-98"),
  mean_cs        = c(1.484, 1.526, 3.756, 1.780, 2.627,
                     1.730, 0.000, 0.000, 0.000, 0.000),
  signif_n       = c(2, 3, 2, 7, 3,
                     1, 0, 0, 0, 0),
  category_n     = c(40, 44, 27, 70, 27,
                     12, 8, 11, 17, 7)
)

bar_fill  <- "#552c17ff"
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

p_mean_cs_icd10 <- ggplot(
  df_icd10 %>%
    mutate(
      total_tests = category_n * 2,
      label = paste0(icd10_category, " (", signif_n, "/", total_tests, ")")
    ) %>%
    arrange(mean_cs) %>%
    mutate(label = factor(label, levels = label)),
  aes(x = mean_cs, y = label)
) +
  geom_col(width = bar_width,
           fill = bar_fill,
           color = bar_edge,
           linewidth = bar_lwd,
           alpha = 0.7) +
  scale_x_continuous(position = "top") +
  labs(
    x = "Mean combined score",
    y = "ICD-10 category (significant / all tested associations)"
  ) +
  theme_classic() +
  theme_black_text

svg(
  filename = file.path("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplots", "barplot_meanCS_neural_26.02.2026.svg"),
  width = 5,
  height = 3
)

p_mean_cs_icd10

dev.off()

# ##############################################################################
# ---- blood ----
# ##############################################################################
# ---- dane (ręcznie) ----
df_icd10 <- tibble(
  icd10_category = c("F30-39","F20-29","F00-09","F10-19","F40-49",
                     "F80-89","F50-59","F60-69","F70-79","F90-98"),
  mean_cs        = c(2.382, 2.864, 2.099, 1.456, 1.346389,
                     0.000, 0.000, 0.000, 0.000, 0.000),
  signif_n       = c(4, 1, 2, 4, 1,
                     0, 0, 0, 0, 0),
  category_n     = c(70, 27, 40, 44, 27,
                     17, 12, 8, 11, 7)
)

# ---- styl (jak wcześniej) ----
bar_fill  <- "#552c17ff"
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# ---- wykres ----
p_mean_cs_icd10 <- ggplot(
  df_icd10 %>%
    mutate(
      total_tests = category_n * 2,
      label = paste0(icd10_category, " (", signif_n, "/", total_tests, ")")
    ) %>%
    arrange(mean_cs) %>%
    mutate(label = factor(label, levels = label)),
  aes(x = mean_cs, y = label)
) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge,
           linewidth = bar_lwd, alpha = 0.7) +
  scale_x_continuous(position = "top") +
  labs(
    x = "Mean combined score",
    y = "ICD-10 category (significant / all tested associations)"
  ) +
  theme_classic() +
  theme_black_text

svg(
  filename = file.path("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplots", "barplot_meanCS_blood_26.02.2026.svg"),
  width = 5,
  height = 3
)

p_mean_cs_icd10

dev.off()


# ##############################################################################
# ---- lung ----
# ##############################################################################\
library(tidyverse)

# ---- dane (zestaw 3) ----
df_icd10 <- tibble(
  icd10_category = c("F30-39","F20-29","F00-09","F10-19","F40-49",
                     "F80-89","F50-59","F60-69","F70-79","F90-98"),
  mean_cs        = c(2.220, 3.021, 4.736, 1.895, 1.774,
                     2.116, 0.000, 0.000, 0.000, 4.846),
  signif_n       = c(12, 2, 1, 2, 2,
                     2, 0, 0, 0, 2),
  category_n     = c(70, 27, 40, 44, 27,
                     17, 12, 8, 11, 7)
)

# ---- styl ----
bar_fill  <- "#552c17ff"
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# ---- wykres ----
p_mean_cs_icd10 <- ggplot(
  df_icd10 %>%
    mutate(
      total_tests = category_n * 2,
      label = paste0(icd10_category, " (", signif_n, "/", total_tests, ")")
    ) %>%
    arrange(mean_cs) %>%
    mutate(label = factor(label, levels = label)),
  aes(x = mean_cs, y = label)
) +
  geom_col(width = bar_width,
           fill = bar_fill,
           color = bar_edge,
           linewidth = bar_lwd,
           alpha = 0.7) +
  scale_x_continuous(position = "top") +
  labs(
    x = "Mean combined score",
    y = "ICD-10 category (significant / all tested associations)"
  ) +
  theme_classic() +
  theme_black_text


svg(
  filename = file.path("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplots", "barplot_meanCS_lung_26.02.2026.svg"),
  width = 5,
  height = 3
)

p_mean_cs_icd10

dev.off()



# ##############################################################################
# ##############################################################################

# ##############################################################################
# ---- Neural ----
# ##############################################################################
df_icd10 <- tibble(
  icd10_category = c("F00-09","F10-19","F20-29","F30-39","F40-49",
                     "F50-59","F60-69","F70-79","F80-89","F90-98"),
  mean_cs        = c(1.484, 1.526, 3.756, 1.780, 2.627,
                     1.730, 0.000, 0.000, 0.000, 0.000),
  signif_n       = c(2, 3, 2, 7, 3,
                     1, 0, 0, 0, 0),
  category_n     = c(40, 44, 27, 70, 27,
                     12, 8, 11, 17, 7)
)

bar_fill  <- "#552c17ff"
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

p_mean_cs_icd10 <- ggplot(
  df_icd10 %>%
    mutate(
      total_tests = category_n * 2,
      label = paste0(icd10_category, " (", signif_n, "/", total_tests, ")")
    ) %>%
    arrange(mean_cs) %>%
    mutate(label = factor(label, levels = label)),
  aes(x = mean_cs, y = label)
) +
  geom_col(width = bar_width,
           fill = bar_fill,
           color = bar_edge,
           linewidth = bar_lwd,
           alpha = 0.7) +
  scale_x_continuous(
    position = "top",
    limits = c(0, 5)
  ) +
  labs(
    x = "Mean combined score",
    y = "ICD-10 category (significant / all tested associations)"
  ) +
  theme_classic() +
  theme_black_text

svg(
  filename = file.path("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplots", "barplot_meanCS_neural_xAxisRangeFive_26.02.2026.svg"),
  width = 5,
  height = 3
)

p_mean_cs_icd10

dev.off()

# ##############################################################################
# ---- blood ----
# ##############################################################################
# ---- dane (ręcznie) ----
df_icd10 <- tibble(
  icd10_category = c("F30-39","F20-29","F00-09","F10-19","F40-49",
                     "F80-89","F50-59","F60-69","F70-79","F90-98"),
  mean_cs        = c(2.382, 2.864, 2.099, 1.456, 1.346389,
                     0.000, 0.000, 0.000, 0.000, 0.000),
  signif_n       = c(4, 1, 2, 4, 1,
                     0, 0, 0, 0, 0),
  category_n     = c(70, 27, 40, 44, 27,
                     17, 12, 8, 11, 7)
)

# ---- styl (jak wcześniej) ----
bar_fill  <- "#552c17ff"
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# ---- wykres ----
p_mean_cs_icd10 <- ggplot(
  df_icd10 %>%
    mutate(
      total_tests = category_n * 2,
      label = paste0(icd10_category, " (", signif_n, "/", total_tests, ")")
    ) %>%
    arrange(mean_cs) %>%
    mutate(label = factor(label, levels = label)),
  aes(x = mean_cs, y = label)
) +
  geom_col(width = bar_width, fill = bar_fill, color = bar_edge,
           linewidth = bar_lwd, alpha = 0.7) +
  scale_x_continuous(
    position = "top",
    limits = c(0, 5)
  ) +
  labs(
    x = "Mean combined score",
    y = "ICD-10 category (significant / all tested associations)"
  ) +
  theme_classic() +
  theme_black_text

svg(
  filename = file.path("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplots", "barplot_meanCS_blood_xAxisRangeFive_26.02.2026.svg"),
  width = 5,
  height = 3
)

p_mean_cs_icd10

dev.off()


# ##############################################################################
# ---- lung ----
# ##############################################################################\
library(tidyverse)

# ---- dane (zestaw 3) ----
df_icd10 <- tibble(
  icd10_category = c("F30-39","F20-29","F00-09","F10-19","F40-49",
                     "F80-89","F50-59","F60-69","F70-79","F90-98"),
  mean_cs        = c(2.220, 3.021, 4.736, 1.895, 1.774,
                     2.116, 0.000, 0.000, 0.000, 4.846),
  signif_n       = c(12, 2, 1, 2, 2,
                     2, 0, 0, 0, 2),
  category_n     = c(70, 27, 40, 44, 27,
                     17, 12, 8, 11, 7)
)

# ---- styl ----
bar_fill  <- "#552c17ff"
bar_edge  <- "black"
bar_width <- 0.7
bar_lwd   <- 1

theme_black_text <- theme(
  plot.title   = element_text(color = "black"),
  axis.title.x = element_text(color = "black"),
  axis.title.y = element_text(color = "black"),
  axis.text.x  = element_text(color = "black"),
  axis.text.y  = element_text(color = "black"),
  axis.ticks   = element_line(color = "black"),
  axis.line    = element_line(color = "black")
)

# ---- wykres ----
p_mean_cs_icd10 <- ggplot(
  df_icd10 %>%
    mutate(
      total_tests = category_n * 2,
      label = paste0(icd10_category, " (", signif_n, "/", total_tests, ")")
    ) %>%
    arrange(mean_cs) %>%
    mutate(label = factor(label, levels = label)),
  aes(x = mean_cs, y = label)
) +
  geom_col(width = bar_width,
           fill = bar_fill,
           color = bar_edge,
           linewidth = bar_lwd,
           alpha = 0.7) +
  scale_x_continuous(
    position = "top",
    limits = c(0, 5)
  ) +
  labs(
    x = "Mean combined score",
    y = "ICD-10 category (significant / all tested associations)"
  ) +
  theme_classic() +
  theme_black_text


svg(
  filename = file.path("/home/mateusz/projects/ifpan-michkor-gr/results_v2/overlap/overlap_grTissues_GWASCatalogDisGeNETgenebass/figures/barplots", "barplot_meanCS_lung_xAxisRangeFive_26.02.2026.svg"),
  width = 5,
  height = 3
)

p_mean_cs_icd10

dev.off()