# Figure 2: where HR and LR differ. a, pathways in the two interaction contrasts; b, c, the
# OXPHOS catch-up over training and the opposite responses to the bout, per person; d, the black
# (sarcomere) module's link to strength gain in each arm; e, whether the labels can be predicted
# at all. Reads saved outputs only; per-person values come from saved scores and eigengenes.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})
source(here("R", "helpers.R"))
source(here("R", "figure_theme.R"))
source(here("R", "figure_panels.R"))

stage <- here("06_Figures", "F02_Groups")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  set_tests = "03_Pathways/03_Contrasts/c_data/set_tests.rds",
  gene_sets = "03_Pathways/01_Gene_Sets/c_data/gene_sets.rds",
  scores = "03_Pathways/02_Scores/c_data/singscore.rds",
  modules = "04_Modules/01_Build/c_data/modules.rds",
  module_associate = "04_Modules/04_Associate/c_data/04_associate.xlsx",
  panel = "05_Integration/03_Panel/c_data/03_panel.xlsx",
  protein_classify = "02_Proteins/03_Classify/c_data/03_classify.xlsx",
  set_classify = "03_Pathways/04_Classify/c_data/04_classify.xlsx",
  module_classify = "04_Modules/03_Classify/c_data/03_classify.xlsx",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Missing inputs: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
read_sheet <- \(input, sheet) readxl::read_excel(paths[[input]], sheet)
set_tests <- readRDS(paths[["set_tests"]])$set_tests
gene_sets <- readRDS(paths[["gene_sets"]])
modules <- readRDS(paths[["modules"]])
sample_meta <- modules$meta
oxphos_id <- "Hallmark|HALLMARK_OXIDATIVE_PHOSPHORYLATION"
oxphos_test <- \(contrast_, method_) {
  filter(set_tests, set_id == oxphos_id, contrast == contrast_, method == method_)
}
arm_scales <- list(
  scale_colour_manual(values = arm_colours, name = NULL),
  scale_shape_manual(values = arm_shapes, name = NULL)
)
shade_window <- \(from, to, family) {
  annotate(
    "rect",
    xmin = from, xmax = to, ymin = -Inf, ymax = Inf,
    fill = contrast_tone("fill")[[family]], alpha = 0.5
  )
}

# a: pathways at fgsea FDR < 0.05 in either interaction contrast, the six strongest per contrast
# after dropping sets that share half their measured proteins with a stronger one.
interactions <- c(Training = "Training_Interaction", Acute = "Acute_Interaction")
members <- map(gene_sets$sets, \(s) intersect(s, gene_sets$gene_universe))
fgsea <- set_tests |>
  filter(method == "fgsea", contrast %in% interactions, collection != "GO_Slim")
chosen <- map(interactions, \(cn) {
  ranked <- filter(fgsea, contrast == cn, fdr < 0.05) |> arrange(desc(abs(nes)))
  head(prune_overlapping(ranked$set_id, members), 6)
}) |>
  unlist() |>
  unique()
interaction_sets <- fgsea |>
  filter(set_id %in% chosen) |>
  mutate(
    window = factor(names(interactions)[match(contrast, interactions)], names(interactions)),
    label = pathway_label(pathway),
    log_fdr = pmin(-log10(fdr), 12)
  )
row_order <- interaction_sets |>
  filter(window == "Acute") |>
  arrange(nes) |>
  pull(label)
nes_limit <- ceiling(max(abs(interaction_sets$nes)) * 10) / 10
panel_interaction <- ggplot(
  mutate(interaction_sets, label = factor(label, row_order)), aes(window, label)
) +
  geom_point(aes(size = log_fdr, fill = nes), shape = 21, stroke = 0.2) +
  scale_fill_distiller(
    palette = "PiYG", direction = 1, limits = c(-nes_limit, nes_limit), name = "NES"
  ) +
  scale_size_area(
    max_size = 3.5, limits = c(0, 12), breaks = c(2, 6, 12),
    labels = c("2", "6", "12+"), name = "-log10 FDR"
  ) +
  labs(
    x = "interaction contrast", y = NULL, title = "Where HR and LR responses differ",
    subtitle = "six strongest fgsea sets per contrast"
  ) +
  figure_theme +
  theme(axis.text.y = element_text(size = 5.5))

# b, c: each person's Hallmark OXPHOS score (singscore) at T1, T2 and T3.
scores <- readRDS(paths[["scores"]])$scores
oxphos <- tibble(sample_id = colnames(scores), score = scores[oxphos_id, ]) |>
  inner_join(sample_meta, by = "sample_id")
oxphos_means <- summarise(oxphos, score = mean(score), .by = c(arm, timepoint))
catch_up <- oxphos_test("Baseline_HRvLR", "fgsea")
training_gain <- oxphos_test("Training_Interaction", "fgsea")
panel_catch_up <- ggplot(oxphos, aes(timepoint, score, colour = arm, shape = arm)) +
  shade_window(1, 2, "Training") +
  shade_window(2, 3, "Acute") +
  geom_line(aes(group = subject), linewidth = 0.25, alpha = 0.5) +
  geom_point(size = 1, alpha = 0.7) +
  geom_line(aes(group = arm), oxphos_means, linewidth = 1.1) +
  arm_scales +
  scale_x_discrete(expand = expansion(add = 0.05)) +
  labs(
    x = NULL, y = "OXPHOS score (Hallmark, singscore)", title = "OXPHOS catch-up",
    subtitle = sprintf(
      "HR - LR at T1 NES %.2f; training int. %+.2f", catch_up$nes, training_gain$nes
    )
  ) +
  figure_theme

acute_change <- oxphos |>
  select(subject, arm, timepoint, score) |>
  pivot_wider(names_from = timepoint, values_from = score) |>
  transmute(subject, arm, change = T3 - T2) |>
  filter(!is.na(change))
acute_fgsea <- oxphos_test("Acute_Interaction", "fgsea")
acute_fry <- oxphos_test("Acute_Interaction", "fry")
panel_acute <- ggplot(acute_change, aes(arm, change, colour = arm, shape = arm)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(position = position_jitter(width = 0.12, seed = 1), size = 1.4) +
  stat_summary(fun = median, geom = "crossbar", width = 0.5, linewidth = 0.3, show.legend = FALSE) +
  arm_scales +
  labs(
    x = NULL, y = "OXPHOS score, T3 - T2", title = "Acute OXPHOS response",
    subtitle = sprintf(
      "interaction NES %.2f; fry p = %.2f", acute_fgsea$nes, acute_fry$p
    )
  ) +
  figure_theme +
  theme(legend.position = "none")

# d: the black module's training change against leg-press gain, one line per arm.
phenotype <- phenotype_tables(readr::read_csv(paths[["phenotype"]], show_col_types = FALSE))
black <- tibble(
  sample_id = colnames(modules$eigengenes), eigengene = modules$eigengenes["black", ]
) |>
  inner_join(sample_meta, by = "sample_id") |>
  filter(timepoint %in% c("T1", "T2")) |>
  select(subject, arm, timepoint, eigengene) |>
  pivot_wider(names_from = timepoint, values_from = eigengene) |>
  transmute(subject, arm, change = T2 - T1) |>
  inner_join(select(phenotype$pairings$training, subject, pct_leg_press_1rm), by = "subject") |>
  filter(!is.na(change), !is.na(pct_leg_press_1rm))
arm_test <- read_sheet("module_associate", "arm_difference") |>
  filter(feature == "black", window == "training", outcome == "pct_leg_press_1rm")
flip <- read_sheet("module_associate", "flips") |>
  filter(feature == "black", window == "training", outcome == "pct_leg_press_1rm")
panel_slope <- ggplot(black, aes(pct_leg_press_1rm, change, colour = arm, shape = arm)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.6) +
  geom_point(size = 1.5) +
  arm_scales +
  labs(
    x = "leg-press 1RM gain (%)", y = "black (sarcomere) module change, T2 - T1",
    title = "Sarcomere-strength link by arm",
    subtitle = sprintf(
      "rho HR %.2f, LR %+.2f; FDR %.2g", flip$hr_rho, flip$lr_rho, arm_test$fdr
    )
  ) +
  figure_theme

# e: cross-validated HR vs LR classification against 200 label permutations, for the module and
# protein panels and the DIABLO models. Single features are summarised in the subtitle.
panels <- read_sheet("panel", "panels") |>
  transmute(model = paste(features, input), auc = auc_mean, null_median, null_95, permutation_p)
diablo <- read_sheet("panel", "diablo") |>
  transmute(
    model = paste("DIABLO", input, if_else(phenotype, "+ phenotype", "")),
    auc = auc_mean, null_median, null_95, permutation_p
  )
classification <- bind_rows(panels, diablo) |>
  mutate(model = factor(trimws(model), rev(trimws(model))))
single_features <- tibble(
  level = c("proteins", "pathways", "modules"),
  input = c("protein_classify", "set_classify", "module_classify")
) |>
  mutate(fdr_05 = map_int(input, \(i) sum(read_sheet(i, "counts")$n_fdr_05)))
panel_classify <- ggplot(classification, aes(y = model)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", linewidth = 0.25) +
  geom_linerange(aes(xmin = null_median, xmax = null_95), linewidth = 2, colour = "grey85") +
  geom_point(aes(auc), size = 1.4) +
  scale_x_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(
    x = "cross-validated AUC", y = NULL, title = "HR and LR cannot be told apart",
    subtitle = "grey: permutation median to 95th percentile"
  ) +
  figure_theme +
  theme(axis.text.y = element_text(size = 5.5))

# Every plotted number traces to its source table before anything is written.
stopifnot(
  nrow(interaction_sets) == 2 * length(chosen), length(chosen) >= 6,
  all(count(oxphos, arm, timepoint)$n %in% c(7, 8)), nrow(acute_change) == 15,
  nrow(catch_up) == 1, nrow(training_gain) == 1, nrow(acute_fgsea) == 1, nrow(acute_fry) == 1,
  nrow(arm_test) == 1, nrow(flip) == 1, nrow(black) == 14,
  nrow(classification) == 18
)

figure <- wrap_plots(
  panel_interaction, panel_catch_up, panel_acute, panel_slope, panel_classify,
  nrow = 1
) +
  plot_annotation(tag_levels = "a") &
  figure_tags &
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
save_figure(figure, figure_dir, "f02_groups", width = 300, height = 110)

sheets <- list(
  interaction_sets = select(interaction_sets, window, set_id, collection, label, nes, fdr),
  oxphos_scores = select(oxphos, sample_id, subject, arm, timepoint, score),
  acute_change = acute_change,
  black_change = black,
  classification = classification,
  single_features = select(single_features, level, fdr_05)
)
holds <- c(
  "Panel a: plotted sets per interaction contrast, with NES and fgsea FDR.",
  "Panel b: each sample's Hallmark OXPHOS singscore.",
  "Panel c: each person's T3 - T2 change in that score.",
  "Panel d: each person's black module training change and leg-press gain.",
  "Panel e: cross-validated AUC and permutation null per panel and DIABLO model.",
  "Single-feature HR vs LR tests at FDR < 0.05, summed over tasks, per level."
)
write_workbook(sheets, holds, file.path(out, "f02_groups.xlsx"), inputs, paths)
message("wrote f02_groups.pdf, f02_groups.png and f02_groups.xlsx")
