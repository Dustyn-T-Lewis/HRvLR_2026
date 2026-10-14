# Figure 3: people
#
# Figure 3 asks whether one feature can tell HR from LR, and whether a protein's link to strength
# gain differs between them. It reads saved outputs only; each panel's rows come from the rule its
# builder states.
#
# | panel | shows |
# |---|---|
# | a | the best single feature per classification task, with its AUC interval |
# | b | the group-by-trait difference against the continuous growth-by-strength interaction |
# | c | within-group Spearman rho, HR against LR, for the same proteins |
# | d | the lead group-difference protein, person by person |
#
# AUC is the area under the ROC curve: the chance that a randomly chosen member of one group scores
# higher than one of the other. Module preservation between the groups is in
# `05_Module_Preservation`.
#
# Packages
#
# The panel builders, `within_subject_change()`, `phenotype_tables()` and `save_figure()` come from
# the shared helpers in `00_Setup/Shared_functions/`.

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
# Seven saved outputs: the classification results at protein, set and module level, the HR/LR trait
# associations, the continuous trait associations, the corrected protein matrix and the phenotype
# table. The chunk creates `b_reports/` and `c_data/`, stops if any input is missing, and reads each
# sheet.

groups <- c("HR", "LR")
stage <- here("05_Figures", "HRvLR", "F3_People")
walk(file.path(stage, c("b_reports", "c_data")), dir.create, recursive = TRUE, showWarnings = FALSE)
classify <- "03_Group_comparisons/HRvLR/c_data"
inputs <- c(
  proteins = file.path(classify, "05_protein_classification.xlsx"),
  sets = file.path(classify, "05_set_classification.xlsx"),
  modules = file.path(classify, "05_module_classification.xlsx"),
  associate = "03_Group_comparisons/HRvLR/c_data/06_trait_associations.xlsx",
  continuous = file.path(
    "04_Continuous_no_group/01_Trait_Associations/01_Proteins",
    "c_data/protein_trait_associations.xlsx"
  ),
  corrected = "01_Processing/03_Quantification/c_data/DAList_corrected.rds",
  phenotype = "00_Setup/Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

read_sheet <- \(name, sheet) readxl::read_excel(paths[[name]], sheet, guess_max = 1e5)
feature_levels <- c("proteins", "sets", "modules")
auc <- map(set_names(feature_levels), \(level) read_sheet(level, "auc"))
classify_counts <- map(set_names(feature_levels), \(level) read_sheet(level, "counts"))
group_difference <- read_sheet("associate", "group_difference")
within_group <- read_sheet("associate", "within_group")
interaction <- read_sheet("continuous", "interaction")
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
group_of <- set_names(phenotype$group, phenotype$subject)
dal <- readRDS(paths[["corrected"]])

# Lead protein
#
# The lead protein is the one with the smallest group-difference p. The chunk takes its T1 to T2
# change per person from the corrected matrix and the trait it was tested against, per person.

lead <- group_difference |>
  filter(feature_type == "proteins") |>
  slice_min(p, n = 1, with_ties = FALSE)
change <- within_subject_change(dal$data, as_tibble(dal$metadata), "T1", "T2")[lead$feature, ]
trait_table <- phenotype_tables(phenotype)$pairings[[lead$window]]
trait <- set_names(trait_table[[lead$trait]], trait_table$subject)
lead

# Panels
#
# Each builder takes its tables and returns `list(plot, data, source)`: the plot, the numbers it
# draws, and one row naming the file and rule behind them.

panels <- list(
  a = panel_auc_top(auc, classify_counts, groups, inputs[["proteins"]]),
  b = panel_group_vs_continuous(
    group_difference, interaction, lead$window, lead$trait, groups, inputs[["associate"]]
  ),
  c = panel_within_rho(
    within_group, group_difference, lead$window, lead$trait, groups, inputs[["associate"]]
  ),
  d = panel_group_slope(lead, change, trait, group_of, within_group, groups, inputs[["associate"]])
)

# Validation
#
# Plotted numbers against their source tables. The chunk stops if panel a lacks a row per feature
# type and task or its tested counts differ from the classification counts, panel b drops a protein,
# or the Spearman rho recomputed from panel d's people differs from the saved within-group rho.

lead_rho <- filter(
  within_group,
  feature == lead$feature, window == lead$window, trait == lead$trait
)
tested <- imap(classify_counts, \(counts, level) {
  summarise(counts, n_tested = sum(n_tested), .by = task) |> mutate(feature_type = level)
}) |>
  list_rbind() |>
  inner_join(panels$a$data, by = c("feature_type", "task"))
person_rho <- panels$d$data |>
  summarise(rho = cor(value, trait, method = "spearman"), n = n(), .by = group) |>
  inner_join(lead_rho, by = "group")
stopifnot(
  nrow(panels$a$data) == length(feature_levels) * 5,
  nrow(tested) == nrow(panels$a$data), all(tested$n_tested.x == tested$n_tested.y),
  nrow(panels$b$data) == nrow(dal$data),
  isTRUE(all.equal(person_rho$rho.x, person_rho$rho.y)), all(person_rho$n == person_rho$n_subjects)
)

# Figure
#
# The chunk assembles the four panels, writes `b_reports/f3_people.pdf` at 300 by 175 mm, and prints
# the figure.

figure <- wrap_plots(
  (panels$a$plot | panels$b$plot) + plot_layout(widths = c(1.7, 1)),
  (panels$c$plot | panels$d$plot | plot_spacer()) + plot_layout(widths = c(1, 1.3, 0.35)),
  ncol = 1, heights = c(1, 1)
) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, file.path(stage, "b_reports"), "f3_people", width = 300, height = 175)

# Save
#
# `c_data/f3_people.xlsx` holds one sheet per panel's numbers, a sheet naming each panel's source
# and rule, and the input manifest with md5 sums.

write_workbook(
  c(
    list(panel_sources = list_rbind(map(panels, "source"))),
    set_names(map(panels, "data"), paste0("panel_", names(panels)))
  ),
  c(
    "Per panel: what it shows, the file and sheet it reads, and the rule that picks its rows.",
    "a: best feature per task and feature type: AUC, DeLong and bootstrap CIs, p, n at p < 0.05.",
    "b: per protein, group-difference t and continuous growth x strength t, with p and FDR.",
    "c: per protein, Spearman rho in each group, and the group-difference evidence.",
    "d: the lead protein's training change and the trait, per person."
  ),
  file.path(stage, "c_data", "f3_people.xlsx"), inputs, paths
)
