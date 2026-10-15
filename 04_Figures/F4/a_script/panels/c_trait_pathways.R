# Trait pathways: per window, the representative of each tree cluster at Jaccard height 0.85
# (03_Continuous/02_Trait_Pathways), up to ten by rank, drawn against every trait of the window with
# fgsea NES and FDR. Change windows (traits as % change) on the left, level windows on the right.
library(here)
library(dplyr)
library(ggplot2)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "pathway_grid.R"))

inputs <- c(trait_pathways = "03_Continuous/02_Trait_Pathways/c_data/trait_pathways.rds")
saved <- readRDS(here(inputs[["trait_pathways"]]))
cap <- 10
holds <- "Panel c: gene sets drawn in the trait pathway matrix, fgsea NES and FDR per trait."

# Window names as strips: what the proteome side of each window is.
halves <- list(
  change = c(baseline = "T1 proteome", training = "Training", acute = "Acute"),
  level = c(baseline_level = "T1", post_level = "T2", post_bout_level = "T3")
)
column_strips <- c(change = "Trait change, T1 to T2", level = "Trait level")
window_names <- unlist(unname(halves))
stopifnot(setequal(names(window_names), saved$trait_selection$window))

rows <- saved$trait_selection |>
  filter(representative_0.85) |>
  slice_min(rank, n = cap, by = window) |>
  mutate(label = set_label(pathway, same_sign), order = rank)
cells <- saved$pathways |>
  semi_join(rows, by = c("window", "set_id")) |>
  transmute(window, set_id, pathway, collection, column = trait, nes, fdr)
nes_limit <- ceiling(max(abs(cells$nes), na.rm = TRUE) * 10) / 10

half_grid <- function(half, title = NULL, subtitle = NULL) {
  strip <- halves[[half]]
  d <- filter(cells, window %in% names(strip))
  traits <- unique(d$column)
  pathway_grid(
    mutate(d, window = strip[window]),
    filter(rows, window %in% names(strip)) |>
      mutate(window = factor(strip[window], strip)),
    tibble::tibble(
      column = traits, column_window = column_strips[[half]],
      markup = trait_label(traits, change = FALSE)
    ),
    title = title, subtitle = subtitle, nes_limit = nes_limit, x_angle = 90
  )
}
panel <- half_grid(
  "change",
  title = "Gene sets that follow the traits",
  subtitle = paste0(
    "fgsea on each window-trait protein ranking, everyone pooled; one set per tree cluster ",
    "(Jaccard height 0.85), up to ten per window; * calls of opposite sign"
  )
) +
  half_grid("level") +
  plot_layout(guides = "collect", widths = c(8, 6))

panel_data <- cells |>
  inner_join(select(rows, window, set_id, rank, label), by = c("window", "set_id")) |>
  transmute(window, rank, set_id, collection, label, trait = column, nes, fdr) |>
  arrange(match(window, names(window_names)), rank, trait)

source_values <- inner_join(panel_data, saved$pathways, by = c("window", "set_id", "trait"))
traits_per_window <- count(distinct(saved$pathways, window, trait), window, name = "n_traits")
cells_per_set <- count(panel_data, window, set_id) |> inner_join(traits_per_window, by = "window")
stopifnot(
  all(count(rows, window)$n <= cap), all(rows$representative_0.85),
  nrow(source_values) == nrow(panel_data),
  identical(source_values$nes.x, source_values$nes.y),
  identical(source_values$fdr.x, source_values$fdr.y),
  nrow(cells_per_set) == nrow(rows), all(cells_per_set$n == cells_per_set$n_traits)
)

panel_height <- 100
ggsave(
  here("04_Figures", "F4", "b_reports", "panels", "c_trait_pathways.pdf"), panel,
  width = 180, height = panel_height, units = "mm", device = cairo_pdf
)
