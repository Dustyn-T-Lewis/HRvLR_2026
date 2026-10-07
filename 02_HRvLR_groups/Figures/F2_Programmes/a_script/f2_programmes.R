# Figure 2: which co-expression programmes move with training and the bout, and do HR and LR
# differ in them? Reads saved outputs only; each panel's rows come from the rule its builder states.
#   a  what each module is            b  each module on the nine contrasts
#   c  gene sets called per contrast, fgsea and fry
#   d  the module with the smallest interaction p, person by person

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
source(here("R", "helpers.R"))
source(here("R", "figure_theme.R"))
source(here("R", "figure_panels.R"))

groups <- c("HR", "LR")
stage <- here("02_HRvLR_groups", "Figures", "F2_Programmes")
walk(file.path(stage, c("b_reports", "c_data")), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  design = "02_HRvLR_groups/01_Design_Contrasts/c_data/design_contrasts.xlsx",
  module_tests = "02_HRvLR_groups/04_Module_Contrasts/c_data/module_contrasts.xlsx",
  set_tests = "02_HRvLR_groups/03_Set_Contrasts/c_data/set_tests.rds",
  set_summary = "02_HRvLR_groups/03_Set_Contrasts/c_data/set_contrasts.xlsx",
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

lead <- eigengene_tests |>
  filter(grepl("Interaction", contrast)) |>
  slice_min(p, n = 1, with_ties = FALSE)

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

# Plotted numbers against their source tables.
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

top <- (panels$a$plot | panels$b$plot) + plot_layout(widths = c(1, 1.35))
right <- (panels$c$plot / panels$d$plot) + plot_layout(heights = c(1.3, 1))
figure <- (free(top) | right) +
  plot_layout(widths = c(2.15, 1)) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, file.path(stage, "b_reports"), "f2_programmes", width = 300, height = 175)

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
