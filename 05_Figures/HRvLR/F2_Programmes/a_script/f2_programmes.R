# Figure 2: programmes
#
# Figure 2 asks which co-expression programmes (WGCNA modules) move with training and the bout, and
# whether HR and LR differ in them. It reads saved outputs only; each panel's rows come from the
# rule its builder states.
#
# | panel | shows |
# |---|---|
# | a | what each module is |
# | b | each module on the nine contrasts |
# | c | gene sets called per contrast, fgsea and fry |
# | d | the module with the smallest interaction p, person by person |
#
# A module eigengene is the first principal component of a module's proteins: one value per sample
# that summarises the module. fgsea ranks a gene set against all other proteins and returns a
# normalised enrichment score (NES); fry tests whether a set's proteins move together across people.
#
# Packages
#
# The panel builders and `save_figure()` come from the shared helpers in
# `00_Setup/Shared_functions/`.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(ggh4x)
  library(patchwork)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_functions", "figure_style.R"))
source(here("00_Setup", "Shared_functions", "figure_panels.R"))

# Inputs
#
# Six saved outputs: the contrast order, the module tests, the set tests and their summary, the
# modules and the module labels. The chunk creates `b_reports/` and `c_data/`, stops if any input is
# missing, and reads each table. Labels are sorted largest module first.

groups <- c("HR", "LR")
stage <- here("05_Figures", "HRvLR", "F2_Programmes")
walk(file.path(stage, c("b_reports", "c_data")), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  design = "03_Group_comparisons/HRvLR/c_data/01_design_contrasts.xlsx",
  module_tests = "03_Group_comparisons/HRvLR/c_data/04_module_contrasts.xlsx",
  set_tests = "03_Group_comparisons/HRvLR/c_data/03_set_tests.rds",
  set_summary = "03_Group_comparisons/HRvLR/c_data/03_set_contrasts.xlsx",
  modules = "01_Processing/06_Module_Building/c_data/modules.rds",
  module_labels = "01_Processing/07_Module_Characterisation/c_data/module_characterisation.xlsx"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

contrast_order <- readxl::read_excel(paths[["design"]], "contrasts")$contrast
eigengene_tests <- readxl::read_excel(paths[["module_tests"]], "eigengene_tests")
fry_tests <- readxl::read_excel(paths[["module_tests"]], "fry_tests")
set_tests <- readRDS(paths[["set_tests"]])$set_tests
modules <- readRDS(paths[["modules"]])
labels <- readxl::read_excel(paths[["module_labels"]], "labels") |>
  arrange(desc(n_proteins))

# Lead module
#
# Panel d shows one module: the one with the smallest eigengene p in either interaction contrast.
# Picking the smallest p inflates that module's difference; the panel shows how it spreads across
# people.

lead <- eigengene_tests |>
  filter(grepl("Interaction", contrast)) |>
  slice_min(p, n = 1, with_ties = FALSE)
lead

# Panels
#
# Each builder takes its tables and returns `list(plot, data, source)`: the plot, the numbers it
# draws, and one row naming the file and rule behind them.

panels <- list(
  a = panel_module_identity(labels, inputs[["module_labels"]]),
  b = panel_module_contrasts(
    eigengene_tests, fry_tests, labels$module, groups, inputs[["module_tests"]]
  ),
  c = panel_set_calls(set_tests, contrast_order, groups, inputs[["set_tests"]]),
  d = panel_module_people(
    modules$eigengenes, as_tibble(modules$meta), lead$module, eigengene_tests,
    module_identity(labels), groups, inputs[["modules"]],
    "smallest eigengene p in either interaction contrast"
  )
)

# Validation
#
# Plotted numbers against their source tables. The chunk stops if a module lacks eigengene tests,
# panel b is missing a module-contrast cell, a panel c count differs from `set_summary`, or panel d
# drops a person.

set_summary <- readxl::read_excel(paths[["set_summary"]], "set_summary")
fgsea_totals <- panels$c$data |>
  filter(startsWith(as.character(call), "fgsea")) |>
  inner_join(set_summary, by = "contrast")
fry_totals <- filter(panels$c$data, startsWith(as.character(call), "fry\nFDR")) |>
  inner_join(set_summary, by = "contrast")
stopifnot(
  setequal(labels$module, eigengene_tests$module),
  nrow(panels$b$data) == length(contrast_order) * nrow(labels),
  nrow(fgsea_totals) == length(contrast_order), all(fgsea_totals$n == fgsea_totals$n_fdr_fgsea),
  all(fry_totals$n == fry_totals$n_fdr_fry),
  n_distinct(panels$d$data$subject) == n_distinct(modules$meta$subject)
)

# Figure
#
# The chunk assembles the four panels, writes `b_reports/f2_programmes.pdf` at 300 by 175 mm, and
# prints the figure.

top <- (panels$a$plot | panels$b$plot) + plot_layout(widths = c(1, 1.35))
right <- (panels$c$plot / panels$d$plot) + plot_layout(heights = c(1.3, 1))
figure <- (free(top) | right) +
  plot_layout(widths = c(2.15, 1)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, file.path(stage, "b_reports"), "f2_programmes", width = 300, height = 175)

# Save
#
# `c_data/f2_programmes.xlsx` holds one sheet per panel's numbers, a sheet naming each panel's
# source and rule, and the input manifest with md5 sums.

write_workbook(
  c(
    list(panel_sources = list_rbind(map(panels, "source"))),
    set_names(map(panels, "data"), paste0("panel_", names(panels)))
  ),
  c(
    "Per panel: what it shows, the file and sheet it reads, and the rule that picks its rows.",
    "a: module size, top set and its FDR, three top members, STRING enrichment ratio.",
    "b: eigengene logFC, p and FDR, and fry p and FDR, per module and contrast.",
    "c: sets called up and down per contrast by fgsea (FDR), fry (FDR) and fry (p).",
    "d: each person's eigengene for the module with the smallest interaction p."
  ),
  file.path(stage, "c_data", "f2_programmes.xlsx"), inputs, paths
)
