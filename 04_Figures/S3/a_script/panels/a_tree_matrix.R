# Every candidate set per block (02_HRvLR/02_Set_Contrasts) in its tree's leaf order, in all nine
# contrasts. Grey and white bands alternate between tree clusters at Jaccard height 0.85; the
# cluster's representative, the set F3 draws when it is among the first eight, is in bold.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(set_tests = "02_HRvLR/02_Set_Contrasts/c_data/set_tests.rds")
saved <- readRDS(here(inputs[["set_tests"]]))
holds <- paste(
  "Every candidate set per block in tree leaf order, with its cluster at height 0.85 and",
  "whether it represents it, in all nine contrasts: fgsea NES and FDR within the database."
)

fgsea <- filter(saved$set_tests, method == "fgsea")
columns <- contrast_columns(unique(fgsea$contrast))
rows <- set_selection(saved$selection) |>
  left_join(distinct(fgsea, set_id, pathway), by = "set_id") |>
  arrange(tree_order) |>
  mutate(
    order = tree_order, shade = consecutive_id(cut_0.85) %% 2 == 0, bold = representative_0.85,
    label = set_label(pathway, same_sign, interaction), .by = block
  )
cells <- fgsea |>
  inner_join(
    select(rows, window, block, set_id, tree_order, cut_0.85, representative_0.85),
    by = "set_id", relationship = "many-to-many"
  ) |>
  transmute(window, block, set_id, pathway, collection,
    column = contrast, nes, fdr, tree_order,
    cut_0.85, representative_0.85
  )

panel <- pathway_grid(
  cells, rows, columns,
  title = "Every candidate set, grouped by tree",
  subtitle = paste0(
    "Sets in leaf order of each block's tree (shared proteins, Jaccard, average linkage); ",
    "bands alternate between clusters at height 0.85\nbold: the cluster's most specific set; ",
    "* calls of opposite sign; \u2020 interaction FDR < 0.05"
  ),
  text_size = 5
)
panel_data <- cells |>
  left_join(select(rows, window, set_id, label), by = c("window", "set_id")) |>
  transmute(
    window = as.character(window), block = sub("\n", ": ", block), tree_order, cut_0.85,
    representative_0.85, set_id, collection, label, contrast = column, nes, fdr
  ) |>
  arrange(match(window, panel_families), block, tree_order, match(contrast, columns$column))

source_values <- inner_join(panel_data, fgsea, by = c("set_id", "contrast"))
stopifnot(
  nrow(rows) == nrow(saved$selection), all(count(panel_data, window, set_id)$n == 9),
  nrow(source_values) == nrow(panel_data),
  identical(source_values$nes.x, source_values$nes.y),
  identical(source_values$fdr.x, source_values$fdr.y),
  all(count(distinct(rows, block, cut_0.85, shade), block, cut_0.85)$n == 1)
)

panel_height <- 2.7 * nrow(rows) + 34
ggsave(
  here("04_Figures", "S3", "b_reports", "panels", "a_tree_matrix.pdf"), panel,
  width = 180, height = panel_height, units = "mm", device = cairo_pdf
)
