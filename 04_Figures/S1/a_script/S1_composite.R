# S1: the technical checks behind F1. Row 1: clustering, variance shares, PC drivers. Row 2: DiffVar
# and the T2 and T3 PCAs.
library(here)
library(patchwork)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S1")
stopifnot(names(panels) == c(
  "a_no_clusters", "b_variance_shares", "c_diffvar", "d_pc_covariates", "e_pca_levels"
))
# Panels are not aligned with each other, so a wrapped panel's title stays beside its tag.
plots <- lapply(panels, \(env) free(env$panel))

# Widths in 5 mm columns of a 250 mm page: a 110, b 70, d 70 mm; then c and e 125 mm each.
layout <- c(
  area(1, 1, 1, 22), area(1, 23, 1, 36), area(2, 1, 2, 25), area(1, 37, 1, 50),
  area(2, 26, 2, 50)
)
figure <- wrap_plots(plots, design = layout, heights = c(95, 72)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, "S1", 167)

write_figure_workbook(panels, "S1")

# The PC-covariate scatter pages behind panel d, a separate PDF.
sys.source(here("04_Figures", "S1", "a_script", "extra", "pc_scatters.R"), envir = new.env())
