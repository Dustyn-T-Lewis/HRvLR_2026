# F4: everyone pooled. Trait changes and links against chance, then the trait pathways across the
# page. The protein profile groups are in S4.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("F4")
stopifnot(names(panels) == c(
  "a_trait_changes", "b_links_vs_chance", "c_trait_pathways"
))
plots <- purrr::map2(lapply(panels, `[[`, "panel"), letters[seq_along(panels)], tag_panel)

# Widths in 5 mm columns of a 250 mm page: a 85 mm beside b 165 mm, then c full width.
# free() keeps a and b from aligning their axes to the long set names of panel c.
layout <- c(area(1, 1, 1, 17), area(1, 18, 1, 50), area(2, 1, 2, 50))
figure <- wrap_plots(
  free(plots$a_trait_changes), free(plots$b_links_vs_chance), plots$c_trait_pathways,
  design = layout, heights = c(70, 100)
)
save_figure(figure, "F4", 170)

write_figure_workbook(panels, "F4")
