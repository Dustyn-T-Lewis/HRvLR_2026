# GO Slim overview: every GO Slim theme (02_HRvLR/02_Set_Contrasts; tested but left out of the
# figure selection) in the nine contrasts, rows ordered by a clustering of their NES profiles.
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds")
saved <- readRDS(here(inputs[["set_tests"]]))
holds <- "Every GO Slim theme in the nine contrasts: fgsea NES and FDR within GO Slim."

fgsea <- filter(saved$set_tests, method == "fgsea", collection == "GO_Slim")
columns <- contrast_columns(unique(fgsea$contrast))
# Themes with similar NES across the nine contrasts sit together; a missing score counts as 0.
profiles <- fgsea |>
  select(set_id, contrast, nes) |>
  pivot_wider(names_from = contrast, values_from = nes, values_fill = 0) |>
  tibble::column_to_rownames("set_id") |>
  as.matrix()
profiles[is.na(profiles)] <- 0
tree <- hclust(dist(profiles), method = "average")
rows <- tibble(
  window = "GO Slim", set_id = rownames(profiles), order = order(tree$order)
)
cells <- transmute(fgsea,
  window = "GO Slim", set_id, pathway, collection, column = contrast,
  nes, fdr
)

panel <- pathway_grid(
  cells, rows, columns,
  title = "GO Slim themes in the nine contrasts",
  subtitle = "Every GO Slim theme; rows ordered by average-linkage clustering of their NES"
)
panel_data <- cells |>
  left_join(select(rows, set_id, order), by = "set_id") |>
  transmute(order, set_id, label = pathway_label(pathway), contrast = column, nes, fdr) |>
  arrange(order, match(contrast, columns$column))

source_values <- inner_join(panel_data, fgsea, by = c("set_id", "contrast"))
stopifnot(
  nrow(columns) == 9, all(count(panel_data, set_id)$n == 9),
  nrow(source_values) == nrow(panel_data),
  identical(source_values$nes.x, source_values$nes.y),
  identical(source_values$fdr.x, source_values$fdr.y)
)

panel_height <- 3 * nrow(rows) + 30
ggsave(
  here("04_Figures", "S3", "b_reports", "panels", "d_go_slim.pdf"), panel,
  width = 180, height = panel_height, units = "mm", device = cairo_pdf
)
