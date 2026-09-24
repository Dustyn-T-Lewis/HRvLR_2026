# The modules tested on the nine contrasts by eigengene and by fry, and each module's membership
# read against protein-level significance.
#
# Eigengenes: lmFit on the module eigengenes with the protein design and subject block, the
# within-subject correlation re-estimated on the eigengenes. fry: each module as a protein set
# on the imputed matrix, rotating residuals under the same design, block and a correlation
# estimated on that matrix, as 03_Pathway_Enrichment/01 does for gene sets.
#
# Membership against significance (WGCNA's MM against GS): within each module, Spearman between a
# member's kME and its protein-level moderated t per contrast, and its Spearman rho with each
# phenotype over the training window. A positive value means the module's most central proteins
# moved up most. Members share a module, so these correlations are descriptive, not tests.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(limma)
})

stage <- here("04_Network", "04_test_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  modules = "04_Network/01_build_modules/c_data/modules.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  fit = "02_Differential_Expression/02_Differential/c_data/fit.rds",
  phenotype = "02_Differential_Expression/03_Phenotype/c_data/phenotype.rds",
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_build_modules and stage 02 first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
modules <- readRDS(paths[["modules"]])
d <- readRDS(paths[["design"]])
fit <- readRDS(paths[["fit"]])
phenotype <- readRDS(paths[["phenotype"]])$protein_association
imputed <- readRDS(paths[["imputed"]])

me <- modules$eigengenes
design <- d$dal$design$design_matrix
contrast_matrix <- d$dal$design$contrast_matrix
contrast_names <- colnames(contrast_matrix)
subject <- d$dal$metadata$subject
abundance <- as.matrix(imputed$data)
members <- filter(modules$membership, module != "grey")
module_order <- modules$module_summary$module
outcomes <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_fcsa_mixed", "d_nfibre_mixed",
  "d_nfibre_I", "d_mcsa", "d_1rm_legpress", "d_1rm_ext", "volume_load"
)
stopifnot(
  identical(colnames(me), rownames(design)),
  identical(colnames(abundance), rownames(design)),
  identical(d$dal$metadata$sample_id, rownames(design)),
  identical(rownames(fit$eBayes_fit$t), rownames(abundance))
)


# ---- eigengene contrasts -------------------------------------------------------------------

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


# ---- fry on module sets --------------------------------------------------------------------

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


# ---- membership against significance -------------------------------------------------------

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


# ---- figures -------------------------------------------------------------------------------

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
        "Table: c_data/04_test_modules.xlsx."
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
      "Table: c_data/04_test_modules.xlsx, membership_contrast."
    )
  ) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

pdf(file.path(figure_dir, "04_test_modules_figures.pdf"), width = 11, height = 8.5)
walk(list(
  eigengene_figure, fry_figure, membership_contrast_figure, membership_phenotype_figure,
  primary_figure
), print)
invisible(dev.off())

sheets <- list(
  eigengene_tests = eigengene_tests, fry_tests = fry_tests,
  membership_contrast = membership_contrast, membership_phenotype = membership_phenotype,
  correlation = tibble(
    matrix = c("eigengenes", "imputed"), correlation = c(me_correlation, fry_correlation)
  ),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "limma on module eigengenes: logFC, moderated t, p and BH FDR per contrast.",
  "fry on each module as a protein set: direction, p and FDR per contrast.",
  "Spearman rho between member kME and protein moderated t, per module and contrast.",
  "Spearman rho between member kME and protein rho with each phenotype over training.",
  "Within-subject correlation on the eigengenes and on the imputed matrix.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
writexl::write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "04_test_modules.xlsx"))
message("wrote 04_test_modules.xlsx and 5 figures")
