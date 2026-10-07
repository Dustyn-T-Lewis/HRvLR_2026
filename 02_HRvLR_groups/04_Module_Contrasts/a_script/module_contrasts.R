# The modules tested on the nine contrasts by eigengene and by fry, and each module's membership
# read against protein-level significance.
#
# Eigengenes: lmFit on the module eigengenes with the protein design and subject block, the
# within-subject correlation re-estimated on the eigengenes. fry: each module as a protein set
# on limpa's matrix, rotating residuals under the same design, block and a correlation
# estimated on that matrix, as 03_Set_Contrasts does for gene sets.
#
# Membership against significance (WGCNA's MM against GS): within each module, Spearman between a
# member's kME and its protein-level moderated t per contrast. A positive value means the module's
# most central proteins moved up most. Members share a module, so these correlations are
# descriptive, not tests. Membership against phenotype belongs to the continuous arc.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
})
source(here("R", "helpers.R"))

stage <- here("02_HRvLR_groups", "04_Module_Contrasts")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  modules = "01_Processing/06_Module_Building/c_data/modules.rds",
  design = "02_HRvLR_groups/01_Design_Contrasts/c_data/design.rds",
  fit = "02_HRvLR_groups/02_Protein_Contrasts/c_data/fit.rds",
  limpa = "01_Processing/03_Quantification/c_data/DAList_limpa.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
modules <- readRDS(paths[["modules"]])
d <- readRDS(paths[["design"]])
fit <- readRDS(paths[["fit"]])
complete <- readRDS(paths[["limpa"]])

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
message(
  "correlation: eigengenes ", round(me_correlation, 3), ", limpa matrix ",
  round(fry_correlation, 3)
)

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

report_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    plot.caption = element_text(hjust = 0, size = 7, colour = "grey45")
  )
# The label marks nominal p and a black border BH < 0.05.
module_tiles <- function(data, columns, fill_label, limits, title, subtitle, caption) {
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
    labs(
      x = NULL, y = NULL, fill = fill_label, title = title, subtitle = subtitle,
      caption = stringr::str_wrap(width = 150, paste(
        caption, "Label: nominal p where below 0.05. Black border: BH < 0.05.",
        "Table: c_data/module_contrasts.xlsx."
      ))
    ) +
    report_theme
}

eigengene_figure <- module_tiles(
  rename(eigengene_tests, column = contrast, effect = t), contrast_names, "moderated t", NULL,
  "Module eigengenes across the nine contrasts",
  sprintf(
    "limma on eigengenes, protein design, subject block (correlation %.3f); BH within contrast",
    me_correlation
  ),
  "Fill: moderated t of each eigengene in each contrast."
)
fry_signed <- fry_tests |>
  mutate(effect = -log10(p) * if_else(direction == "Up", 1, -1))
fry_figure <- module_tiles(
  rename(fry_signed, column = contrast), contrast_names, "signed\n-log10 p", NULL,
  "Modules as protein sets across the nine contrasts",
  sprintf(
    "limma::fry on limpa's matrix, protein design, subject block (correlation %.3f)",
    fry_correlation
  ),
  paste(
    "Fill: -log10 fry p, positive when the module moved up. fry's FDR runs over the",
    nrow(me), "modules within each contrast."
  )
)
membership_contrast_figure <- module_tiles(
  rename(membership_contrast, column = contrast, effect = rho) |> mutate(fdr = NA_real_),
  contrast_names, "Spearman\nrho", c(-1, 1),
  "Module membership against protein significance",
  "Spearman between each member's kME and its protein-level moderated t, per module and contrast",
  paste(
    "Fill: rho; positive when the module's most central proteins moved up most. Members share",
    "a module, so the p describes rather than tests."
  )
)

primary <- members |>
  inner_join(filter(protein_t, contrast == "Training_Interaction"), by = "uniprot_id") |>
  left_join(
    filter(membership_contrast, contrast == "Training_Interaction") |> select(module, rho),
    by = "module"
  ) |>
  mutate(module = factor(module, levels = module_order)) |>
  arrange(module) |>
  mutate(panel = sprintf("%s  rho %+.2f", module, rho), panel = factor(panel, unique(panel)))
primary_figure <- ggplot(primary, aes(kme, t)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(aes(fill = module), shape = 21, colour = "grey30", size = 1.3, alpha = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linewidth = 0.4) +
  facet_wrap(~panel, ncol = 4) +
  scale_fill_identity() +
  labs(
    x = "kME (membership in the module)", y = "moderated t, Training_Interaction",
    title = "Membership against the primary contrast",
    subtitle = "one point per module member; Spearman rho in each header",
    caption = paste(
      "Line: least-squares fit, for reference.",
      "Table: c_data/module_contrasts.xlsx, membership_contrast."
    )
  ) +
  report_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

sample_values <- me
column_meta <- modules$meta
feature_label <- set_names(module_order)
feature_noun <- "modules"
contrast_results <- transmute(eigengene_tests, feature = module, contrast, effect = logFC, p, fdr)
row_groups <- list(modules = module_order)

# Overviews: every feature across the 45 biopsies, then across the nine contrasts. Biopsy traits
# sit on the T1 and T2 columns; T3 has none. comp_hypertrophy is one value per subject.
report_arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
trait_table <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
biopsy_traits <- trait_table |>
  select(subject, matches("_(pre|post)_")) |>
  pivot_longer(-subject, names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*") |>
  pivot_wider(names_from = trait) |>
  transmute(
    subject,
    timepoint = if_else(when == "pre", "T1", "T2"),
    fcsa_mixed, mcsa, leg_press_1rm
  )
columns <- as_tibble(column_meta) |>
  transmute(sample_id, subject, arm = as.character(arm), timepoint = as.character(timepoint)) |>
  mutate(group = factor(paste(arm, timepoint, sep = "_"), paste(
    rep(c("HR", "LR"), each = 3), c("T1", "T2", "T3"),
    sep = "_"
  ))) |>
  left_join(biopsy_traits, by = c("subject", "timepoint")) |>
  left_join(select(trait_table, subject, comp_hypertrophy), by = "subject") |>
  arrange(group, subject)
stopifnot(!anyDuplicated(columns$sample_id), setequal(columns$sample_id, colnames(sample_values)))
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
row_order_note <- "Rows follow module size."
effect_limit <- cap(contrast_results$effect)
overviews <- imap(row_groups, \(rows, group) {
  m <- sample_values[rows, columns$sample_id, drop = FALSE]
  limit <- cap(m)
  labelled_heat <- \(...) heat(..., labels = feature_label)
  biopsy_map <- labelled_heat(m, "biopsies", c(-limit, 0, limit), "eigengene",
    row_names = TRUE, column_title = character(0), top_annotation = sample_annotation,
    column_split = columns$group, cluster_column_slices = FALSE,
    column_labels = columns$subject, column_names_gp = gpar(fontsize = 5)
  )
  contrast_grid <- \(value) as_grid(contrast_results, rows, contrast_names, "contrast", value)
  contrast_maps <- labelled_heat(
    contrast_grid(\(d) d$effect), "eigengene log2 fold change",
    c(-effect_limit, 0, effect_limit), "eigengene log2 FC"
  ) +
    labelled_heat(contrast_grid(\(d) status_of(d$p, d$fdr)), "moderated t p and BH q",
      row_names = TRUE
    )
  n_rows <- format(length(rows), big.mark = ",")
  list(
    heatmap_page(biopsy_map, sprintf("Biopsy overview: %s", group), paste(
      sprintf("All %s %s by %d biopsies,", n_rows, feature_noun, ncol(m)),
      "grouped by arm and time point; column labels are subjects.",
      "Fill: module eigengene.",
      sprintf("Colours stop at the 99th percentile of |value|, %.2g.", limit),
      row_order_note, "Grey in the trait tracks: not measured (T3 has no traits)."
    ), legends = list()),
    heatmap_page(contrast_maps, sprintf("Contrast overview: %s", group), paste(
      sprintf("All %s %s by the nine contrasts.", n_rows, feature_noun),
      "Fill: eigengene fold change from limma; status: p and BH q within contrast.",
      sprintf("Colours stop at the 99th percentile of |value|, %.2g.", effect_limit),
      row_order_note, "Table: c_data/module_contrasts.xlsx, eigengene_tests."
    ))
  )
}) |>
  list_flatten()

pdf(file.path(figures, "module_contrasts_report.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(list(
  eigengene_figure, fry_figure, membership_contrast_figure,
  primary_figure
), print)
invisible(dev.off())

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
  file.path(out, "module_contrasts.xlsx"),
  inputs, paths
)
