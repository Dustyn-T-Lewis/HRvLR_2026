# Figure 3: growth and strength as continuous traits, with no HR/LR labels. a, every continuous
# screen against its shuffle null; b-d, the one result that passes BH: sarcomere proteins whose
# training change depends on growth x leg-press gain, its robustness, the people it leans on and
# its pathway echo; e, two proteins that track strength at baseline; f, whether people cluster.
# Reads saved outputs only.

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

stage <- here("06_Figures", "F03_Continuum")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  protein_associate = "02_Proteins/04_Associate/c_data/04_associate.xlsx",
  set_associate = "03_Pathways/05_Associate/c_data/05_associate.xlsx",
  module_associate = "04_Modules/04_Associate/c_data/04_associate.xlsx",
  continuous = "05_Integration/04_Continuous/c_data/04_continuous.xlsx",
  clustering = "05_Integration/05_Clustering/c_data/05_clustering.xlsx",
  modules = "04_Modules/01_Build/c_data/modules.rds",
  corrected = "01_Preprocess/03_Quantify/c_data/DAList_corrected.rds",
  gene_sets = "03_Pathways/01_Gene_Sets/c_data/gene_sets.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Missing inputs: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
read_sheet <- \(input, sheet) readxl::read_excel(paths[[input]], sheet)
phenotype <- phenotype_tables(readr::read_csv(paths[["phenotype"]], show_col_types = FALSE))
growth <- select(phenotype$pairings$training, subject, comp_hypertrophy, pct_leg_press_1rm)
short_subject <- \(x) sub("^(HR|LR)_", "", x)

# a: every continuous screen's count against 200 shuffles of whole phenotype rows, as observed /
# shuffle 95th percentile; above 1 beats the null. Proteins count p < 0.05, pathways count sets
# at FDR < 0.05.
screens <- bind_rows(
  read_sheet("protein_associate", "chance") |>
    transmute(level = "proteins", window, outcome, observed = nominal, null_95, p = perm_p_count),
  read_sheet("continuous", "pathway_chance") |>
    transmute(level = "pathways", window, outcome, observed = sets_fdr_05, null_95, p = perm_p)
) |>
  mutate(
    ratio = observed / null_95, beats = p < 0.05,
    level = factor(level, c("proteins", "pathways")),
    screen = paste0(level, ": ", window, " vs ", outcome)
  )
panel_screens <- ggplot(screens, aes(level, ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.25) +
  geom_point(
    aes(colour = beats),
    position = position_jitter(width = 0.15, seed = 1), size = 1.2
  ) +
  ggrepel::geom_text_repel(
    aes(label = sprintf("%s vs %s\np = %.2f", window, outcome, p)), filter(screens, beats),
    size = 1.8, seed = 1, nudge_x = 0.3, min.segment.length = 0, segment.size = 0.2
  ) +
  scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "black"), guide = "none") +
  labs(
    x = NULL, y = "observed / shuffle 95th percentile", title = "Screens against chance",
    subtitle = sprintf(
      "%d of %d continuous screens beat 200 phenotype shuffles", sum(screens$beats), nrow(screens)
    )
  ) +
  figure_theme

# b: the proteins at BH in the training interaction (leg press), before and after adjusting for
# fibre composition (magenta and yellow module changes; 05_Integration, composition).
interaction <- read_sheet("protein_associate", "interaction") |>
  filter(window == "training", strength_test == "pct_leg_press_1rm", fdr < 0.05)
composition <- read_sheet("continuous", "composition") |>
  filter(
    level == "proteins", window == "training", strength_test == "pct_leg_press_1rm",
    feature %in% interaction$feature
  )
adjusted <- interaction |>
  select(feature, label, fdr, loo_max_p, loo_worst_subject) |>
  left_join(select(composition, feature, fdr_adjusted), by = "feature") |>
  mutate(label = reorder(label, -log10(fdr)), survives = fdr_adjusted < 0.05)
panel_proteins <- ggplot(adjusted, aes(y = label)) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) +
  geom_segment(aes(x = -log10(fdr), xend = -log10(fdr_adjusted)), colour = "grey70") +
  geom_point(aes(x = -log10(fdr)), size = 1.3, colour = direction_colours[["Down"]]) +
  geom_point(aes(x = -log10(fdr_adjusted), shape = survives), size = 1.3, colour = "grey30") +
  scale_shape_manual(
    values = c(`TRUE` = 16, `FALSE` = 1), labels = c("lost", "kept"),
    name = "after composition"
  ) +
  labs(
    x = "-log10 FDR", y = NULL, title = "Sarcomere proteins track growth x strength",
    subtitle = sprintf(
      "%d at FDR < 0.05 (pink), %d after composition; all negative",
      nrow(adjusted), sum(adjusted$survives)
    )
  ) +
  figure_theme +
  theme(axis.text.y = element_text(size = 5.2))

# c: the black module's training change against leg-press gain, coloured by growth, with the
# people the protein results lean on (leave-one-out) named.
modules <- readRDS(paths[["modules"]])
black <- tibble(
  sample_id = colnames(modules$eigengenes), eigengene = modules$eigengenes["black", ]
) |>
  inner_join(modules$meta, by = "sample_id") |>
  filter(timepoint %in% c("T1", "T2")) |>
  select(subject, arm, timepoint, eigengene) |>
  pivot_wider(names_from = timepoint, values_from = eigengene) |>
  transmute(subject, arm, change = T2 - T1) |>
  inner_join(growth, by = "subject") |>
  filter(!is.na(change), !is.na(pct_leg_press_1rm))
drivers <- count(adjusted, subject = loo_worst_subject, name = "proteins")
black <- left_join(black, drivers, by = "subject") |>
  mutate(name = if_else(!is.na(proteins), sprintf("%s (%d)", short_subject(subject), proteins), NA))
module_test <- read_sheet("module_associate", "interaction") |>
  filter(feature == "black", window == "training", strength_test == "pct_leg_press_1rm")
panel_people <- ggplot(black, aes(pct_leg_press_1rm, change)) +
  geom_hline(yintercept = 0, linewidth = 0.25) +
  geom_point(aes(fill = comp_hypertrophy, shape = arm), size = 2, stroke = 0.3) +
  ggrepel::geom_text_repel(
    aes(label = name), filter(black, !is.na(name)),
    size = 1.9, seed = 1, min.segment.length = 0, segment.size = 0.2
  ) +
  scale_fill_distiller(palette = "Purples", direction = 1, name = "growth\n(composite)") +
  scale_shape_manual(values = c(HR = 21, LR = 24), name = NULL) +
  labs(
    x = "leg-press 1RM gain (%)", y = "black (sarcomere) module change, T2 - T1",
    title = "Who carries it",
    subtitle = sprintf(
      "module FDR %.2g; named: left out, they weaken most proteins (n)", module_test$fdr
    )
  ) +
  figure_theme

# d: the pathways at BH in the same interaction.
gene_sets <- readRDS(paths[["gene_sets"]])
pathways <- read_sheet("set_associate", "interaction") |>
  filter(window == "training", strength_test == "pct_leg_press_1rm", fdr < 0.05) |>
  mutate(
    source = c(GOBP = "GO", Reactome = "R", Hallmark = "H", GO_Slim = "GO Slim", KEGG_Legacy = "K")[
      sub("[|].*", "", feature)
    ],
    name = paste0(
      pathway_label(gene_sets$set_catalog$pathway[match(feature, gene_sets$set_catalog$set_id)]),
      " [", source, "]"
    ),
    name = reorder(name, -log10(fdr)),
    direction = if_else(effect > 0, "Up", "Down")
  )
panel_pathways <- ggplot(pathways, aes(-log10(fdr), name, fill = direction)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", linewidth = 0.25) +
  scale_fill_manual(values = direction_colours, guide = "none") +
  labs(
    x = "-log10 FDR", y = NULL, title = "Pathway echo",
    subtitle = sprintf(
      "%d sets at FDR < 0.05; %d negative", nrow(pathways), sum(pathways$direction == "Down")
    )
  ) +
  figure_theme +
  theme(axis.text.y = element_text(size = 5.2))

# e: two proteins whose T1 abundance tracks leg-press strength at T1, across all 15 people.
baseline_hits <- read_sheet("protein_associate", "change_score") |>
  filter(window == "baseline_level", outcome == "leg_press_1rm", label %in% c("ENO2", "UQCR11"))
baseline_chance <- read_sheet("protein_associate", "chance") |>
  filter(window == "baseline_level", outcome == "leg_press_1rm")
corrected <- readRDS(paths[["corrected"]])
t1 <- corrected$metadata |> filter(timepoint == "T1")
baseline_values <- baseline_hits |>
  select(feature, label, rho, fdr) |>
  pmap(\(feature, label, rho, fdr) {
    tibble(
      label = sprintf("%s (rho %+.2f)", label, rho), subject = t1$subject,
      abundance = corrected$data[feature, t1$sample_id]
    )
  }) |>
  list_rbind() |>
  inner_join(select(phenotype$pairings$baseline_level, subject, leg_press_1rm), by = "subject") |>
  left_join(distinct(modules$meta, subject, arm), by = "subject")
panel_baseline <- ggplot(baseline_values, aes(abundance, leg_press_1rm)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.4, colour = "grey40") +
  geom_point(aes(colour = arm, shape = arm), size = 1.3) +
  facet_wrap(~label, scales = "free_x") +
  scale_colour_manual(values = arm_colours, name = NULL) +
  scale_shape_manual(values = arm_shapes, name = NULL) +
  labs(
    x = "log2 abundance at T1", y = "leg-press 1RM at T1 (kg)",
    title = "Proteins tracking baseline strength",
    subtitle = sprintf(
      "FDR %.2g; strongest hit beats shuffles, familywise p = %.3f",
      max(baseline_hits$fdr), baseline_chance$perm_p_strongest
    )
  ) +
  figure_theme

# f: person clusters per input against a max-statistic null, and the best agreement with HR/LR.
clusters <- read_sheet("clustering", "summary") |>
  filter(
    analysis == "people" | (analysis == "groups" & statistic == "best HR/LR ARI")
  ) |>
  mutate(
    facet = if_else(
      analysis == "people", "cluster strength (silhouette)", "agreement with HR/LR (ARI)"
    ),
    screen = if_else(
      analysis == "people", sub(" k [0-9]+$", "", screen), paste("pathway groups", screen)
    ),
    beats = perm_p < 0.05
  )
panel_clusters <- ggplot(clusters, aes(value, screen)) +
  geom_linerange(
    aes(xmin = coalesce(null_median, 0), xmax = null_95),
    linewidth = 2, colour = "grey85"
  ) +
  geom_point(aes(colour = beats), size = 1.4) +
  facet_wrap(~facet, ncol = 1, scales = "free") +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "black"), guide = "none") +
  labs(
    x = "observed (grey: null median to 95th percentile)", y = NULL,
    title = "People cluster, not by arm",
    subtitle = "black: beats the max-statistic null"
  ) +
  figure_theme +
  theme(axis.text.y = element_text(size = 5.2))

# Every plotted number traces to its source table before anything is written.
stopifnot(
  nrow(screens) == 72,
  nrow(adjusted) == 20, !anyNA(adjusted$fdr_adjusted), !anyNA(adjusted$loo_worst_subject),
  all(adjusted$loo_max_p < 0.05),
  nrow(black) == 14, nrow(module_test) == 1,
  nrow(pathways) == 23,
  nrow(baseline_hits) == 2, nrow(baseline_chance) == 1, nrow(baseline_values) == 2 * nrow(t1),
  nrow(clusters) == 15
)

figure <- wrap_plots(
  panel_screens, panel_proteins, panel_people,
  panel_pathways, panel_baseline, panel_clusters,
  nrow = 2
) +
  plot_annotation(tag_levels = "a") &
  figure_tags &
  theme(legend.position = "bottom", legend.key.size = unit(2.5, "mm"))
save_figure(figure, figure_dir, "f03_continuum", width = 300, height = 190)

sheets <- list(
  screens = select(screens, level, window, outcome, observed, null_95, ratio, p),
  interaction_proteins = select(
    adjusted, feature, label, fdr, fdr_adjusted, loo_max_p, loo_worst_subject
  ),
  black_change = select(black, subject, arm, change, pct_leg_press_1rm, comp_hypertrophy, proteins),
  interaction_pathways = select(pathways, feature, name, effect, fdr),
  baseline_strength = baseline_values,
  clustering = select(clusters, analysis, screen, statistic, value, null_median, null_95, perm_p)
)
holds <- c(
  "Panel a: per screen, observed count, shuffle 95th percentile, their ratio and permutation p.",
  "Panel b: proteins at BH, FDR before and after composition, leave-one-out max p and person.",
  "Panel c: each person's black module training change, leg-press gain, growth, driver count.",
  "Panel d: pathways at BH in the training x leg-press interaction.",
  "Panel e: ENO2 and UQCR11 at T1 against leg-press 1RM at T1.",
  "Panel f: clustering statistics against their max-statistic nulls."
)
write_workbook(sheets, holds, file.path(out, "f03_continuum.xlsx"), inputs, paths)
message("wrote f03_continuum.pdf, f03_continuum.png and f03_continuum.xlsx")
