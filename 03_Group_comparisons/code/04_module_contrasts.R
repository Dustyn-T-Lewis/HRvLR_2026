#' ---
#' title: "Module contrasts"
#' params:
#'   grouping: HRvLR
#' format:
#'   typst:
#'     toc: true
#'     margin:
#'       x: 1.6cm
#' execute:
#'   echo: true
#'   warning: false
#' ---
#'
#' The co-expression modules from `06_Module_Building`, tested in each of the nine contrasts two
#' ways:
#'
#' - by eigengene: limma on each module's one-number summary per biopsy, with the protein design,
#'   person block and a within-person correlation re-estimated on the eigengenes;
#' - by fry: each module as a protein set on limpa's matrix, as step 03 tests gene sets.
#'
#' A third table reads each module's membership against protein-level change: within a module, the
#' Spearman correlation between a member's kME (how closely it follows the module) and its moderated
#' t. Members share a module, so these correlations describe; they do not test.
#'
#' ## Packages

#| label: packages
#| message: false
library(here)
library(limma)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(grid)
library(ComplexHeatmap)
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))
params$grouping

#' ## Inputs

#| label: inputs
data_dir <- file.path("03_Group_comparisons", params$grouping, "c_data")
inputs <- c(
  modules = "01_Processing/06_Module_Building/c_data/modules.rds",
  design = file.path(data_dir, "01_design.rds"),
  fit = file.path(data_dir, "02_fit.rds"),
  limpa = "01_Processing/03_Quantification/c_data/DAList_limpa.rds",
  phenotype = "00_Setup/Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
d <- readRDS(paths[["design"]])
modules <- subset_for_design(readRDS(paths[["modules"]]), d)
fit <- readRDS(paths[["fit"]])
complete <- subset_for_design(readRDS(paths[["limpa"]]), d)
me <- modules$eigengenes
design <- d$dal$design$design_matrix
contrast_matrix <- d$dal$design$contrast_matrix
contrast_names <- colnames(contrast_matrix)
subject <- d$dal$metadata$subject
abundance <- complete$data
members <- filter(modules$membership, module != "grey")
module_order <- modules$module_summary$module
stopifnot(
  identical(colnames(me), rownames(design)),
  identical(colnames(abundance), rownames(design)),
  identical(d$dal$metadata$sample_id, rownames(design)),
  identical(rownames(fit$eBayes_fit$t), rownames(abundance))
)
dim(me)

#' ## Eigengene tests
#'
#' Each eigengene goes through the same linear model as a protein: `lmFit()` with person as a block,
#' the nine contrasts, and robust empirical Bayes. The within-person correlation is estimated on the
#' eigengenes themselves (`duplicateCorrelation()`), because they are the values this model fits.

#| label: eigengene-tests
me_correlation <- duplicateCorrelation(me, design, block = subject)$consensus.correlation
me_fit <- lmFit(me, design, block = subject, correlation = me_correlation) |>
  contrasts.fit(contrast_matrix) |>
  eBayes(robust = TRUE)
eigengene_tests <- map(contrast_names, \(contrast) {
  topTable(me_fit, coef = contrast, number = Inf, sort.by = "none", confint = TRUE) |>
    rownames_to_column("module") |>
    transmute(
      contrast, module, logFC,
      logFC_lower = CI.L, logFC_upper = CI.R, t, p = P.Value, fdr = adj.P.Val,
      qvalue = storey_q(p)
    )
}) |>
  list_rbind()
me_correlation
eigengene_tests |>
  slice_min(p, n = 5) |>
  select(contrast, module, logFC, p, fdr)

#' ## Modules as protein sets

#| label: fry-tests
fry_correlation <- duplicateCorrelation(abundance, design, block = subject)$consensus.correlation
module_rows <- split(match(members$uniprot_id, rownames(abundance)), members$module)
fry_tests <- map(contrast_names, \(contrast) {
  fry(abundance, module_rows, design, contrast_matrix[, contrast],
    block = subject, correlation = fry_correlation, sort = "none"
  ) |>
    rownames_to_column("module") |>
    transmute(
      contrast, module,
      n_proteins = NGenes, direction = Direction, p = PValue, fdr = FDR, qvalue = storey_q(p)
    )
}) |>
  list_rbind()
fry_tests |>
  slice_min(p, n = 5)

#' ## Membership against protein change
#'
#' WGCNA calls this module membership against gene significance. A positive rho means the module's
#' most central proteins moved up most in that contrast. Each rho carries a Fisher-z 95% interval
#' (`DescTools::CorCI()`).

#| label: membership
protein_t <- as_tibble(fit$eBayes_fit$t, rownames = "uniprot_id") |>
  pivot_longer(-uniprot_id, names_to = "contrast", values_to = "t")
membership_contrast <- members |>
  inner_join(protein_t, by = "uniprot_id", relationship = "one-to-many") |>
  summarise(
    n_members = n(),
    test = list(cor.test(kme, t, method = "spearman", exact = FALSE)),
    .by = c(module, contrast)
  ) |>
  mutate(
    rho = map_dbl(test, "estimate"),
    rho_lower = map2_dbl(rho, n_members, \(r, n) DescTools::CorCI(r, n)[["lwr.ci"]]),
    rho_upper = map2_dbl(rho, n_members, \(r, n) DescTools::CorCI(r, n)[["upr.ci"]]),
    p = map_dbl(test, "p.value"), test = NULL
  )

#' ## Figures
#'
#' `module_tiles()` draws one module-by-column grid: fill is the effect, the label is the nominal p
#' where below 0.05, and a black border marks BH < 0.05. It is used for all three tables.

#| label: module-tiles
report_theme <- theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
module_tiles <- function(data, columns, fill_label, limits) {
  data |>
    mutate(
      module = factor(module, levels = rev(module_order)),
      column = factor(column, levels = columns)
    ) |>
    ggplot(aes(column, module, fill = effect)) +
    geom_tile(colour = "white") +
    geom_tile(
      data = \(x) filter(x, !is.na(fdr), fdr < 0.05),
      colour = "black", linewidth = 0.8, fill = NA
    ) +
    geom_text(aes(label = if_else(p < 0.05, sprintf("%.3f", p), "")), size = 2.6) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = limits) +
    labs(x = NULL, y = NULL, fill = fill_label) +
    report_theme
}

#| label: fig-eigengene
#| fig-cap: "Eigengene moderated t per module and contrast."
#| fig-width: 8
#| fig-height: 5
module_tiles(
  rename(eigengene_tests, column = contrast, effect = t), contrast_names, "moderated t", NULL
)

#| label: fig-fry
#| fig-cap: "fry per module and contrast: -log10 p, positive when the module moved up."
#| fig-width: 8
#| fig-height: 5
fry_tests |>
  mutate(effect = -log10(p) * if_else(direction == "Up", 1, -1)) |>
  rename(column = contrast) |>
  module_tiles(contrast_names, "signed\n-log10 p", NULL)

#| label: fig-membership
#| fig-cap: "Spearman rho between member kME and protein moderated t, per module and contrast."
#| fig-width: 8
#| fig-height: 5
membership_contrast |>
  rename(column = contrast, effect = rho) |>
  mutate(fdr = NA_real_) |>
  module_tiles(contrast_names, "Spearman\nrho", c(-1, 1))

#| label: fig-primary
#| fig-cap: "Each member's kME against its Training_Interaction t; rho in each header."
#| fig-width: 10
#| fig-height: 8
members |>
  inner_join(filter(protein_t, contrast == "Training_Interaction"), by = "uniprot_id") |>
  left_join(
    filter(membership_contrast, contrast == "Training_Interaction") |> select(module, rho),
    by = "module"
  ) |>
  mutate(module = factor(module, levels = module_order)) |>
  arrange(module) |>
  mutate(panel = sprintf("%s  rho %+.2f", module, rho), panel = factor(panel, unique(panel))) |>
  ggplot(aes(kme, t)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(aes(fill = module), shape = 21, colour = "grey30", size = 1.3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linewidth = 0.4) +
  facet_wrap(~panel, ncol = 4) +
  scale_fill_identity() +
  labs(x = "kME (membership in the module)", y = "moderated t, Training_Interaction") +
  theme_minimal(base_size = 10)

#' ## Every module at a glance
#'
#' Every eigengene across biopsies, then its fold change across the nine contrasts. Rows follow
#' module size.

#| label: fig-overviews
#| fig-cap: "Every module across biopsies, then across the nine contrasts."
#| fig-width: 11
#| fig-height: 8.5
feature_label <- set_names(module_order)
contrast_results <- transmute(eigengene_tests, feature = module, contrast, effect = logFC, p, fdr)
report_arm_colours <- set_names(c("#2166AC", "#B2182B"), d$arms)
trait_table <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
biopsy_traits <- trait_table |>
  select(subject, matches("_(pre|post)_")) |>
  pivot_longer(-subject, names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*") |>
  pivot_wider(names_from = trait) |>
  transmute(
    subject,
    timepoint = if_else(when == "pre", "T1", "T2"), fcsa_mixed, mcsa, leg_press_1rm
  )
columns <- as_tibble(modules$meta) |>
  transmute(sample_id, subject, arm = as.character(arm), timepoint = as.character(timepoint)) |>
  mutate(group = factor(
    paste(arm, timepoint, sep = "_"),
    paste(rep(d$arms, each = 3), c("T1", "T2", "T3"), sep = "_")
  )) |>
  left_join(biopsy_traits, by = c("subject", "timepoint")) |>
  left_join(select(trait_table, subject, comp_hypertrophy), by = "subject") |>
  arrange(group, subject)
stopifnot(!anyDuplicated(columns$sample_id), setequal(columns$sample_id, colnames(me)))
grey_ramp <- \(x) circlize::colorRamp2(range(x, na.rm = TRUE), c("#F7F7F7", "#252525"))
tracks <- list(
  "fCSA mixed (um2)" = columns$fcsa_mixed,
  "mCSA (cm2)" = columns$mcsa, "leg press 1RM (kg)" = columns$leg_press_1rm,
  "comp_hypertrophy, per subject" = columns$comp_hypertrophy
)
sample_annotation <- HeatmapAnnotation(
  df = data.frame(arm = columns$arm, timepoint = columns$timepoint, tracks, check.names = FALSE),
  col = c(
    list(arm = report_arm_colours, timepoint = c(T1 = "grey85", T2 = "grey55", T3 = "grey25")),
    map(tracks, grey_ramp)
  ),
  na_col = "grey80", annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)
effect_limit <- cap(contrast_results$effect)
m <- me[module_order, columns$sample_id, drop = FALSE]
limit <- cap(m)
labelled_heat <- \(...) make_report_heatmap(..., labels = feature_label)
biopsy_map <- labelled_heat(m, "biopsies", c(-limit, 0, limit), "eigengene",
  row_names = TRUE, column_title = character(0), top_annotation = sample_annotation,
  column_split = columns$group, cluster_column_slices = FALSE,
  column_labels = columns$subject, column_names_gp = gpar(fontsize = 5)
)
contrast_grid <- \(value) {
  results_to_matrix(contrast_results, module_order, contrast_names, "contrast", value)
}
contrast_maps <- labelled_heat(
  contrast_grid(\(d) d$effect), "eigengene log2 fold change",
  c(-effect_limit, 0, effect_limit), "eigengene log2 FC"
) +
  labelled_heat(contrast_grid(\(d) status_of(d$p, d$fdr)), "moderated t p and BH q",
    row_names = TRUE
  )
heatmap_page(biopsy_map, "Biopsy overview: modules", "Module eigengene per biopsy.",
  legends = list()
)$draw()
heatmap_page(
  contrast_maps, "Contrast overview: modules", "Eigengene fold change per contrast."
)$draw()

#' ## Validation

#| label: validate
stopifnot(
  nrow(eigengene_tests) == nrow(me) * length(contrast_names),
  nrow(fry_tests) == nrow(me) * length(contrast_names),
  all(eigengene_tests$fdr >= eigengene_tests$p),
  all(abs(membership_contrast$rho) <= 1)
)

#' ## Summary
#'
#' `r nrow(me)` modules by 9 contrasts: `r sum(eigengene_tests$p < 0.05)` eigengene tests at
#' p < 0.05 and `r sum(eigengene_tests$fdr < 0.05)` at BH < 0.05; fry gives
#' `r sum(fry_tests$p < 0.05)` at p < 0.05 and `r sum(fry_tests$fdr < 0.05)` at BH < 0.05.
#' Within-person correlation `r round(me_correlation, 3)` on the eigengenes,
#' `r round(fry_correlation, 3)` on limpa's matrix.
#'
#' ## Save

#| label: save
out <- here(data_dir)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_workbook(
  list(
    eigengene_tests = eigengene_tests, fry_tests = fry_tests,
    membership_contrast = membership_contrast,
    correlation = tibble(
      matrix = c("eigengenes", "limpa"), correlation = c(me_correlation, fry_correlation)
    )
  ),
  c(
    "limma on module eigengenes: logFC, moderated t, p, BH and Storey q per contrast.",
    "fry on each module as a protein set: direction, p and FDR per contrast.",
    paste(
      "Spearman rho (Fisher-z 95% CI) between member kME and protein moderated t, per module and",
      "contrast."
    ),
    "Within-subject correlation on the eigengenes and on limpa's matrix."
  ),
  file.path(out, "04_module_contrasts.xlsx"),
  inputs, paths
)
