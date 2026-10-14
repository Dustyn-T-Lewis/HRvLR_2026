# Figure 3: can one feature tell HR from LR, and does a protein's link to strength gain differ
# between them? Reads saved outputs only; each panel's rows come from the rule its builder states.
#   a  the best single feature per classification task, with its AUC interval
#   b  the group-by-trait difference against the continuous growth-by-strength interaction
#   c  within-group Spearman rho, HR against LR, for the same proteins
#   d  the lead group-difference protein, person by person
# Module preservation between the groups is in 05_Module_Preservation.

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

groups <- c("HR", "LR")
stage <- here("02_HRvLR_groups", "Figures", "F3_People")
walk(file.path(stage, c("b_reports", "c_data")), dir.create, recursive = TRUE, showWarnings = FALSE)
classify <- "02_HRvLR_groups/06_Classification"
inputs <- c(
  proteins = file.path(classify, "01_Proteins/c_data/protein_classification.xlsx"),
  sets = file.path(classify, "02_Sets/c_data/set_classification.xlsx"),
  modules = file.path(classify, "03_Modules/c_data/module_classification.xlsx"),
  associate = "02_HRvLR_groups/07_Trait_Associations/c_data/trait_associations.xlsx",
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
group_of <- set_names(phenotype$arm, phenotype$subject)
dal <- readRDS(paths[["corrected"]])

lead <- group_difference |>
  filter(level == "proteins") |>
  slice_min(p, n = 1, with_ties = FALSE)
change <- within_subject_change(dal$data, as_tibble(dal$metadata), "T1", "T2")[lead$feature, ]
trait_table <- phenotype_tables(phenotype)$pairings[[lead$window]]
trait <- set_names(trait_table[[lead$trait]], trait_table$subject)

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

# Plotted numbers against their source tables.
lead_rho <- filter(
  within_group,
  feature == lead$feature, window == lead$window, trait == lead$trait
)
tested <- imap(classify_counts, \(counts, level) {
  summarise(counts, n_tested = sum(n_tested), .by = task) |> mutate(level = level)
}) |>
  list_rbind() |>
  inner_join(panels$a$data, by = c("level", "task"))
person_rho <- panels$d$data |>
  summarise(rho = cor(value, trait, method = "spearman"), n = n(), .by = arm) |>
  inner_join(lead_rho, by = c(arm = "group"))
stopifnot(
  nrow(panels$a$data) == length(feature_levels) * 5,
  nrow(tested) == nrow(panels$a$data), all(tested$n_tested.x == tested$n_tested.y),
  nrow(panels$b$data) == nrow(dal$data),
  isTRUE(all.equal(person_rho$rho.x, person_rho$rho.y)), all(person_rho$n == person_rho$n_subjects)
)

figure <- wrap_plots(
  (panels$a$plot | panels$b$plot) + plot_layout(widths = c(1.7, 1)),
  (panels$c$plot | panels$d$plot | plot_spacer()) + plot_layout(widths = c(1, 1.3, 0.35)),
  ncol = 1, heights = c(1, 1)
) +
  plot_annotation(tag_levels = "a") &
  figure_tags
save_figure(figure, file.path(stage, "b_reports"), "f3_people", width = 300, height = 175)

write_workbook(
  c(
    list(panel_sources = list_rbind(map(panels, "source"))),
    set_names(map(panels, "data"), paste0("panel_", names(panels)))
  ),
  c(
    "Per panel: what it shows, the file and sheet it reads, and the rule that picks its rows.",
    "a: best feature per task and level: AUC, DeLong and bootstrap CIs, p, counts at p < 0.05.",
    "b: per protein, group-difference t and continuous growth x strength t, with p and FDR.",
    "c: per protein, Spearman rho in each group, and the group-difference evidence.",
    "d: the lead protein's training change and the trait, per person."
  ),
  file.path(stage, "c_data", "f3_people.xlsx"), inputs, paths
)
