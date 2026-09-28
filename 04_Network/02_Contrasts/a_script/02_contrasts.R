# The modules tested on the nine contrasts by eigengene and by fry, and each module's membership
# read against protein-level significance.
#
# Eigengenes: lmFit on the module eigengenes with the protein design and subject block, the
# within-subject correlation re-estimated on the eigengenes. fry: each module as a protein set
# on the imputed matrix, rotating residuals under the same design, block and a correlation
# estimated on that matrix, as 03_Pathway_Enrichment/02_Contrasts does for gene sets.
#
# Membership against significance (WGCNA's MM against GS): within each module, Spearman between a
# member's kME and its protein-level moderated t per contrast, and its Spearman rho with each
# phenotype over the training window. For a contrast, a positive value means the module's most
# central proteins moved up most; for a phenotype, that they track it most positively. Members
# share a module, so these correlations are descriptive, not tests.

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
  library(writexl)
})

stage <- here("04_Network", "02_Contrasts")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  modules = "04_Network/01_Build_Modules/c_data/modules.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  fit = "02_Differential_Expression/02_Contrasts/c_data/fit.rds",
  associate = "02_Differential_Expression/04_Associate/c_data/associate.rds",
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
modules <- readRDS(paths[["modules"]])
d <- readRDS(paths[["design"]])
fit <- readRDS(paths[["fit"]])
phenotype <- readRDS(paths[["associate"]])$change_score |> rename(uniprot_id = feature)
imputed <- readRDS(paths[["imputed"]])

me <- modules$eigengenes
design <- d$dal$design$design_matrix
contrast_matrix <- d$dal$design$contrast_matrix
contrast_names <- colnames(contrast_matrix)
subject <- d$dal$metadata$subject
abundance <- as.matrix(imputed$data)
members <- filter(modules$membership, module != "grey")
module_order <- modules$module_summary$module
outcomes <- unique(phenotype$outcome)
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
  topTable(me_fit, coef = contrast, number = Inf, sort.by = "none") |>
    rownames_to_column("module") |>
    transmute(contrast, module, logFC, t, p = P.Value, fdr = adj.P.Val)
}) |>
  list_rbind()

fry_correlation <- duplicateCorrelation(abundance, design, block = subject)$consensus.correlation
module_rows <- split(match(members$uniprot_id, rownames(abundance)), members$module)
fry_tests <- map(contrast_names, \(contrast) {
  fry(abundance, module_rows, design, contrast_matrix[, contrast],
    block = subject, correlation = fry_correlation, sort = "none"
  ) |>
    rownames_to_column("module") |>
    transmute(contrast, module, n_proteins = NGenes, direction = Direction, p = PValue, fdr = FDR)
}) |>
  list_rbind()
message(
  "correlation: eigengenes ", round(me_correlation, 3), ", imputed matrix ",
  round(fry_correlation, 3)
)

protein_t <- as_tibble(fit$eBayes_fit$t, rownames = "uniprot_id") |>
  pivot_longer(-uniprot_id, names_to = "contrast", values_to = "t")
member_rho <- function(data, value) {
  data |>
    filter(!is.na(.data[[value]])) |>
    summarise(
      n_members = n(),
      test = list(cor.test(kme, .data[[value]], method = "spearman", exact = FALSE)),
      .by = c(module, any_of(c("contrast", "outcome")))
    ) |>
    mutate(rho = map_dbl(test, "estimate"), p = map_dbl(test, "p.value"), test = NULL)
}
membership_contrast <- members |>
  inner_join(protein_t, by = "uniprot_id", relationship = "one-to-many") |>
  member_rho("t")
membership_phenotype <- members |>
  inner_join(
    filter(phenotype, window == "training") |> select(uniprot_id, outcome, association = rho),
    by = "uniprot_id", relationship = "one-to-many"
  ) |>
  member_rho("association")

figure_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    plot.caption = element_text(hjust = 0, size = 7, colour = "grey45")
  )
# One tile per module and column; the label marks nominal p and a black border BH < 0.05.
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
        "Table: c_data/02_contrasts.xlsx."
      ))
    ) +
    figure_theme
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
    "limma::fry on the imputed matrix, protein design, subject block (correlation %.3f)",
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
membership_phenotype_figure <- module_tiles(
  rename(membership_phenotype, column = outcome, effect = rho) |> mutate(fdr = NA_real_),
  outcomes, "Spearman\nrho", c(-1, 1),
  "Module membership against phenotype association",
  "Spearman between each member's kME and its protein rho with the phenotype over training",
  paste(
    "Fill: rho; positive when central proteins track the phenotype most positively.",
    "Descriptive, as above."
  )
)

primary <- members |>
  inner_join(filter(protein_t, contrast == "Training_Interaction"), by = "uniprot_id") |>
  filter(!is.na(t)) |>
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
      "Table: c_data/02_contrasts.xlsx, membership_contrast."
    )
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

sample_values <- me
column_meta <- modules$meta
feature_label <- set_names(module_order)
feature_noun <- "modules"
contrast_results <- transmute(eigengene_tests, feature = module, contrast, effect = logFC, p, fdr)
row_groups <- list(modules = module_order)
# Heatmap and workbook helpers, the same in every contrasts, classify and associate script.
# Display only: row order, colour limits and centring touch no test.
ht_opt$message <- FALSE
diverging <- c("#5E3C99", "white", "#E66101")
status_colours <- c("not significant" = "white", "p < 0.05" = "#35978F", "BH < 0.05" = "black")
status_legend <- Legend(
  labels = c(names(status_colours), "not tested"), title = "status",
  legend_gp = gpar(fill = c(status_colours, "grey80"))
)
status_of <- \(p, fdr) {
  case_when(fdr < 0.05 ~ "BH < 0.05", p < 0.05 ~ "p < 0.05", !is.na(p) ~ "not significant")
}
# Colours stop at the 99th percentile of |x| so a few extremes do not wash out the rest.
cap <- \(x) unname(quantile(abs(x), 0.99, na.rm = TRUE))
# Long results to a feature-by-column matrix; untested cells stay NA and draw grey.
as_grid <- function(d, rows, cols, col, value) {
  m <- matrix(NA, length(rows), length(cols), dimnames = list(rows, cols))
  d <- filter(d, feature %in% rows, .data[[col]] %in% cols)
  m[cbind(d$feature, d[[col]])] <- value(d)
  m
}
# An effect map on a diverging scale when breaks are given, otherwise a status map.
heat <- function(m, title, breaks = NULL, legend = title, row_names = FALSE,
                 column_title = title, column_names_gp = gpar(fontsize = 7), ...) {
  rows <- rownames(m)
  Heatmap(m,
    name = title, na_col = "grey80",
    col = if (is.null(breaks)) status_colours else circlize::colorRamp2(breaks, diverging),
    show_heatmap_legend = !is.null(breaks) && !is.null(legend),
    heatmap_legend_param = list(title = legend, at = if (!is.null(breaks)) signif(breaks, 2)),
    cluster_rows = FALSE, cluster_columns = FALSE, column_title = column_title,
    column_title_gp = gpar(fontsize = 8), column_names_gp = column_names_gp,
    show_row_names = row_names && length(rows) <= 80,
    row_labels = gsub("\n", " ", feature_label[rows]),
    row_names_gp = gpar(fontsize = if (length(rows) <= 40) 7 else 5),
    use_raster = length(rows) > 200, ...
  )
}
heatmap_page <- function(maps, title, caption, legends = list(status_legend)) {
  list(title = title, draw = function() {
    draw(maps,
      column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
      merge_legend = TRUE, heatmap_legend_list = legends, padding = unit(c(14, 4, 4, 4), "mm")
    )
    grid.text(stringr::str_wrap(caption, 190),
      x = unit(4, "mm"), y = unit(3, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 7)
    )
  })
}
write_workbook <- function(sheets, holds, file) {
  sheets <- c(sheets, list(
    input_manifest = tibble(
      input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
    ),
    package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
      as_tibble() |>
      select(package, version = loadedversion, source)
  ))
  holds <- c(holds, "Files read, with md5.", "Packages loaded at run time.")
  write_xlsx(c(list(read_me = tibble(sheet = names(sheets), holds = holds)), sheets), file)
}

# Overviews: every feature across the 45 biopsies, then across the nine contrasts. Biopsy traits
# sit on the T1 and T2 columns; T3 has none. comp_hypertrophy is one value per subject.
arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
trait_table <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
biopsy_traits <- trait_table |>
  select(subject, matches("_(pre|post)_")) |>
  pivot_longer(-subject, names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*") |>
  pivot_wider(names_from = trait) |>
  transmute(
    subject,
    timepoint = if_else(when == "pre", "T1", "T2"),
    type1_share = 100 * fibres_type1 / fibres_mixed, fcsa_mixed, mcsa, leg_press_1rm
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
  "type I share (%)" = columns$type1_share, "fCSA mixed (um2)" = columns$fcsa_mixed,
  "mCSA (cm2)" = columns$mcsa, "leg press 1RM (kg)" = columns$leg_press_1rm,
  "comp_hypertrophy, per subject" = columns$comp_hypertrophy
)
sample_annotation <- HeatmapAnnotation(
  df = data.frame(arm = columns$arm, timepoint = columns$timepoint, tracks, check.names = FALSE),
  col = c(
    list(arm = arm_colours, timepoint = c(T1 = "grey85", T2 = "grey55", T3 = "grey25")),
    map(tracks, grey_ramp)
  ),
  na_col = "grey80", annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)
row_order <- paste(
  "Rows follow module size."
)
effect_limit <- cap(contrast_results$effect)
overviews <- imap(row_groups, \(rows, group) {
  m <- sample_values[rows, columns$sample_id, drop = FALSE]
  m <- m
  limit <- cap(m)
  biopsy_map <- heat(m, "biopsies", c(-limit, 0, limit), "eigengene",
    row_names = TRUE, column_title = character(0), top_annotation = sample_annotation,
    column_split = columns$group, cluster_column_slices = FALSE,
    column_labels = columns$subject, column_names_gp = gpar(fontsize = 5)
  )
  contrast_grid <- \(value) as_grid(contrast_results, rows, contrast_names, "contrast", value)
  contrast_maps <- heat(contrast_grid(\(d) d$effect), "eigengene log2 fold change",
    c(-effect_limit, 0, effect_limit), "eigengene log2 FC"
  ) +
    heat(contrast_grid(\(d) status_of(d$p, d$fdr)), "moderated t p and BH q", row_names = TRUE)
  n_rows <- format(length(rows), big.mark = ",")
  list(
    heatmap_page(biopsy_map, sprintf("Biopsy overview: %s", group), paste(
      sprintf("All %s %s by %d biopsies,", n_rows, feature_noun, ncol(m)),
      "grouped by arm and time point; column labels are subjects.",
      "Fill: module eigengene.",
      sprintf("Colours stop at the 99th percentile of |value|, %.2g.", limit),
      row_order, "Grey: missing."
    ), legends = list()),
    heatmap_page(contrast_maps, sprintf("Contrast overview: %s", group), paste(
      sprintf("All %s %s by the nine contrasts.", n_rows, feature_noun),
      "Fill: eigengene fold change from limma; status: p and BH q within contrast.",
      sprintf("Colours stop at the 99th percentile of |value|, %.2g.", effect_limit),
      row_order, "Grey: not tested. Table: c_data/02_contrasts.xlsx, eigengene_tests."
    ))
  )
}) |>
  list_flatten()
stopifnot(length(overviews) == 2 * length(row_groups))

pdf(file.path(figures, "02_contrasts_figures.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(list(
  eigengene_figure, fry_figure, membership_contrast_figure, membership_phenotype_figure,
  primary_figure
), print)
invisible(dev.off())

write_workbook(
  list(
    eigengene_tests = eigengene_tests, fry_tests = fry_tests,
    membership_contrast = membership_contrast, membership_phenotype = membership_phenotype,
    correlation = tibble(
      matrix = c("eigengenes", "imputed"), correlation = c(me_correlation, fry_correlation)
    )
  ),
  c(
    "limma on module eigengenes: logFC, moderated t, p and BH FDR per contrast.",
    "fry on each module as a protein set: direction, p and FDR per contrast.",
    "Spearman rho between member kME and protein moderated t, per module and contrast.",
    "Spearman rho between member kME and protein rho with each phenotype over training.",
    "Within-subject correlation on the eigengenes and on the imputed matrix."
  ),
  file.path(out, "02_contrasts.xlsx")
)
