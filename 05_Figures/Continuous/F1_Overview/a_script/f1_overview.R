# Figure 1: continuous overview
#
# Figure 1 for the continuous analysis: every person pooled and each trait kept as a number. It
# reads saved outputs only; the one computation is the proteome PCA (principal component analysis),
# a projection with no test.
#
# | panel | shows |
# |---|---|
# | a | each person's training response, trait by trait |
# | b | share of features linked to each trait (p < 0.05) per window and level, against chance |
# | c | the proteome PCA, coloured by growth and by strength gain |
# | d | gene sets leading each trait's protein ranking over training (fgsea) |
# | e | groups of proteins sharing a pattern of trait links (profile clusters) |
#
# The figure goes to `b_reports/f1_overview.pdf` and every number it draws to
# `c_data/f1_overview.xlsx`.
#
# Packages
#
# `analysis_reporting_helpers.R` supplies `write_workbook()` and `phenotype_tables()`;
# `figure_style.R` the theme, tags, thresholds and label helpers; `figure_panels.R`
# `prune_overlapping()` and `save_figure()`.

library(here)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(patchwork)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
source(here("00_Setup", "Shared_functions", "figure_style.R"))
source(here("00_Setup", "Shared_functions", "figure_panels.R"))

# Inputs
#
# The three trait-association workbooks (proteins, gene sets, modules), the trait pathways, the gene
# set members, the profile clusters, the blood-corrected protein matrix and the phenotype table. The
# chunk creates the output folders and stops if any input is missing.

folder <- "04_Continuous_no_group"
stage <- here("05_Figures", "Continuous", "F1_Overview")
walk(file.path(stage, c("b_reports", "c_data")), dir.create, recursive = TRUE, showWarnings = FALSE)
level_file <- \(level, name) file.path(folder, "01_Trait_Associations", level, "c_data", name)
inputs <- c(
  proteins = level_file("01_Proteins", "protein_trait_associations.xlsx"),
  sets = level_file("02_Sets", "set_trait_associations.xlsx"),
  modules = level_file("03_Modules", "module_trait_associations.xlsx"),
  pathways = file.path(folder, "02_Trait_Pathways/c_data/trait_pathways.xlsx"),
  gene_sets = "01_Processing/04_Gene_Sets/c_data/gene_sets.rds",
  profiles = file.path(folder, "03_Profile_Clusters/c_data/profile_clusters.xlsx"),
  corrected = "01_Processing/03_Quantification/c_data/DAList_corrected.rds",
  phenotype = "00_Setup/Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

# Phenotype
#
# `phenotype_tables()` returns one row per person with their training change scores. The seven
# change traits are the four size measures, the two strength measures and composite hypertrophy.

phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
training <- phenotype_tables(phenotype)$pairings$training
change_traits <- c(
  "pct_fcsa_mixed", "pct_fcsa_type1", "pct_fcsa_type2", "pct_mcsa", "pct_leg_press_1rm",
  "pct_leg_ext_1rm", "comp_hypertrophy"
)

# Panel a: each person's response
#
# Each person's % change per trait, scaled to a z score within trait, people ordered by composite
# hypertrophy.

people <- training |>
  select(subject, all_of(change_traits)) |>
  pivot_longer(-subject, names_to = "trait") |>
  mutate(z = (value - mean(value, na.rm = TRUE)) / sd(value, na.rm = TRUE), .by = trait) |>
  mutate(
    subject = factor(subject, training$subject[order(training$comp_hypertrophy)]),
    trait = factor(trait_label(trait), trait_label(change_traits))
  )
panel_people <- ggplot(people, aes(subject, trait, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_distiller(
    palette = "RdBu", limits = c(-2.5, 2.5), oob = scales::squish, na.value = "grey90",
    name = "z within\ntrait"
  ) +
  labs(
    title = "Each person's training response",
    subtitle = sprintf(
      "%d people, ordered by composite hypertrophy; grey: not measured", nrow(training)
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 5), panel.grid = element_blank())

# Panel b: how many features track each trait
#
# The share of features at p < 0.05 per change-score window and trait, at each level (proteins, gene
# sets, modules). Counts are summed within analysis and trait before the share is taken. The colour
# scale is white at 5%, the share chance gives.

level_names <- c(proteins = "Proteins", sets = "Gene sets", modules = "Modules")
windows <- c(
  training = "training change", baseline = "T1 level", acute = "acute change"
)
links <- imap(level_names, \(name, level) {
  readxl::read_excel(paths[[level]], "counts") |>
    filter(analysis %in% paste("change score:", names(windows))) |>
    summarise(across(c(n_tested, n_nominal), sum), .by = c(analysis, comparison)) |>
    mutate(level = name)
}) |>
  list_rbind() |>
  mutate(
    window = factor(windows[sub("change score: ", "", analysis)], windows),
    trait = trait_label(comparison), share = n_nominal / n_tested,
    level = factor(level, level_names)
  )
panel_links <- ggplot(links, aes(window, trait, fill = share)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  geom_text(aes(label = scales::percent(share, accuracy = 1)), size = 1.7) +
  facet_wrap(~level, nrow = 1) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = alpha_nominal,
    labels = scales::label_percent(), name = "p < 0.05"
  ) +
  labs(
    title = "How many features track each trait",
    subtitle = paste(
      "Share of features at p < 0.05 (limma, one value per person);",
      "white: 5%, the share chance gives"
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())

# Panel c: proteome PCA
#
# PCA of the complete blood-corrected matrix, proteins scaled, coloured by each person's growth and
# strength gain. The two traits sit on different scales, so each is coloured by rank within itself.

dal <- readRDS(paths[["corrected"]])
pcs <- prcomp(t(dal$data), scale. = TRUE)
variance <- summary(pcs)$importance[2, 1:2]
pca <- as_tibble(pcs$x[, 1:2], rownames = "sample_id") |>
  inner_join(select(as_tibble(dal$metadata), sample_id, subject, timepoint), by = "sample_id") |>
  left_join(select(training, subject, comp_hypertrophy, pct_leg_press_1rm), by = "subject") |>
  pivot_longer(c(comp_hypertrophy, pct_leg_press_1rm), names_to = "trait") |>
  mutate(trait = trait_label(trait), rank = percent_rank(value), .by = trait)
panel_pca <- ggplot(pca, aes(PC1, PC2, colour = rank, shape = timepoint)) +
  geom_point(size = 1.3) +
  facet_wrap(~trait, nrow = 1) +
  scale_colour_viridis_c(
    option = "magma", end = 0.9, breaks = c(0, 1), labels = c("low", "high"),
    name = "rank within trait"
  ) +
  scale_shape_manual(values = c(T1 = 16, T2 = 17, T3 = 15), name = NULL) +
  labs(
    title = "Proteome by growth and strength gain",
    subtitle = sprintf("%d biopsies, %s proteins", ncol(dal$data), scales::comma(nrow(dal$data))),
    x = sprintf("PC1 (%.1f%%)", 100 * variance[1]), y = sprintf("PC2 (%.1f%%)", 100 * variance[2])
  ) +
  figure_theme +
  theme(legend.position = "bottom", legend.key.width = unit(8, "mm"))

# Panel d: gene sets leading each trait
#
# fgsea (fast gene set enrichment analysis) ranks proteins by the limma t and gives each set an NES
# (normalised enrichment score) and an FDR. Panel d keeps training-window sets at FDR < 0.05 for at
# least two traits (GO Slim out), ranked by that count and total |NES|. `prune_overlapping()` drops
# any set sharing half its measured proteins with a stronger one, and the first 20 are shown in
# every trait.

gene_sets <- readRDS(paths[["gene_sets"]])
members <- gene_sets$sets
fgsea <- readxl::read_excel(paths[["pathways"]], "pathways") |>
  filter(window == "training", collection != "GO_Slim")
eligible <- fgsea |>
  summarise(
    hits = sum(fdr < alpha_fdr), strength = sum(abs(nes)), .by = c(set_id, pathway)
  ) |>
  filter(hits >= 2) |>
  arrange(desc(hits), desc(strength))
kept <- head(prune_overlapping(eligible$set_id, members), 20)
grid <- fgsea |>
  filter(set_id %in% kept) |>
  mutate(
    set = factor(
      pathway_label(pathway),
      rev(pathway_label(eligible$pathway[match(kept, eligible$set_id)]))
    ),
    trait = trait_label(trait), log_fdr = pmin(-log10(fdr), 8)
  )
nes_limit <- ceiling(max(abs(grid$nes)) * 10) / 10
panel_sets <- ggplot(grid, aes(trait, set)) +
  geom_point(aes(size = log_fdr, fill = nes), shape = 21, stroke = 0.2) +
  scale_fill_distiller(
    palette = "PiYG", direction = 1, limits = c(-nes_limit, nes_limit), name = "NES"
  ) +
  scale_size_area(max_size = 2.8, limits = c(0, 8), name = "-log10 FDR") +
  labs(
    title = "Gene sets leading each trait's ranking, training window",
    subtitle = "fgsea on the limma t; sets at FDR < 0.05 for two or more traits, overlaps pruned",
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.text.y = element_text(size = 5))

# Panel e: profile groups
#
# A profile group is a cluster of proteins that share a pattern of trait links
# (`03_Profile_Clusters`). Panel e draws the twelve largest groups' mean t per trait over training;
# every group and window is in the workbook. Rows read group number (proteins, main module).

profile_means <- readxl::read_excel(paths[["profiles"]], "profile_means")
profile_groups <- readxl::read_excel(paths[["profiles"]], "profile_groups")
profiles <- profile_means |>
  inner_join(select(profile_groups, group, n_proteins, main_module), by = "group") |>
  mutate(
    row = sprintf("%d (%d, %s)", group, n_proteins, main_module),
    window = sub(" .*", "", screen), trait = trait_label(sub("^[^ ]+ ", "", screen))
  )
largest <- slice_max(profile_groups, n_proteins, n = 12, with_ties = FALSE)$group
drawn <- profiles |>
  filter(group %in% largest, window == "training") |>
  mutate(row = reorder(row, n_proteins))
t_limit <- max(abs(drawn$mean_t))
panel_profiles <- ggplot(drawn, aes(trait, row, fill = mean_t)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_distiller(palette = "RdBu", limits = c(-t_limit, t_limit), name = "mean t") +
  labs(
    title = "Proteins that share a pattern of trait links",
    subtitle = paste(
      "Twelve largest profile groups (n proteins, main module):",
      "mean limma t, training window"
    ),
    x = NULL, y = NULL
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), panel.grid = element_blank())

# Figure
#
# Panels a to c on the top row, d and e below.

figure <- wrap_plots(
  (panel_people | panel_links | panel_pca) + plot_layout(widths = c(1.2, 1.5, 1.3)),
  (panel_sets | panel_profiles) + plot_layout(widths = c(1, 1.2)),
  ncol = 1, heights = c(1, 1.15)
) +
  plot_annotation(tag_levels = "a") &
  figure_tags

# Every selected set is drawn in every trait.

stopifnot(nrow(grid) == length(kept) * n_distinct(fgsea$trait))

# `save_figure()` writes `b_reports/f1_overview.pdf` at 300 by 175 mm.

save_figure(figure, file.path(stage, "b_reports"), "f1_overview", width = 300, height = 175)

# Save
#
# One sheet per panel.

write_workbook(
  list(
    people = select(people, subject, trait, value, z),
    links = select(links, level, window, trait, n_tested, n_nominal, share),
    pca = select(pca, sample_id, subject, timepoint, PC1, PC2, trait, value, rank),
    sets = select(grid, set_id, collection, pathway, trait, nes, fdr),
    profiles = select(profiles, group, n_proteins, main_module, screen, mean_t)
  ),
  c(
    "Panel a: each person's % change per trait, and its z score within the trait.",
    "Panel b: features tested and at p < 0.05 per level, change-score window and trait.",
    "Panel c: PC1 and PC2 per biopsy, with the person's growth and strength gain.",
    "Panel d: each selected set in every trait over training: NES and fgsea FDR.",
    "Panel e: each profile group's mean t per window and trait (the figure draws twelve)."
  ),
  file.path(stage, "c_data", "f1_overview.xlsx"), inputs, paths
)
