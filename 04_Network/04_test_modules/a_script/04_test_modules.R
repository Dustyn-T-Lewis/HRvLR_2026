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
# are the ones moving. Members share a module, so these correlations are descriptive, not tests.

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
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

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
phenotype <- readRDS(paths[["phenotype"]])$results
imputed <- readRDS(paths[["imputed"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

me <- modules$eigengenes
design <- d$dal$design$design_matrix
contrast_matrix <- d$dal$design$contrast_matrix
contrast_names <- colnames(contrast_matrix)
subject <- d$dal$metadata$subject
abundance <- as.matrix(imputed$data)
members <- filter(modules$membership, module != "grey")
module_order <- modules$module_summary$module
stopifnot(
  identical(colnames(me), rownames(design)),
  identical(colnames(abundance), rownames(design)),
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
    transmute(contrast, module, n = NGenes, direction = Direction, p = PValue, fdr = FDR)
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
      n = n(),
      rho = cor(kme, .data[[value]], method = "spearman"),
      p = cor.test(kme, .data[[value]], method = "spearman", exact = FALSE)$p.value,
      .by = c(module, any_of(c("contrast", "outcome")))
    )
}
membership_contrast <- members |>
  inner_join(protein_t, by = "uniprot_id", relationship = "one-to-many") |>
  member_rho("t")
membership_phenotype <- members |>
  inner_join(
    filter(phenotype, window == "training") |> select(uniprot_id = protein, outcome, r),
    by = "uniprot_id", relationship = "one-to-many"
  ) |>
  member_rho("r")


# ---- figures -------------------------------------------------------------------------------

save_figure <- function(figure, name, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 200, bg = "white"
    )
  })
}
figure_theme <- theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40")
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

save_figure(
  module_tiles(
    rename(eigengene_tests, column = contrast, effect = t), contrast_names, "moderated t", NULL,
    "Module eigengenes across the nine contrasts",
    sprintf(
      "limma on eigengenes, protein design, subject block (correlation %.3f); BH within contrast",
      me_correlation
    ),
    "Fill: moderated t of each eigengene in each contrast."
  ),
  "01_eigengene_contrasts", 9, 5.5
)
fry_signed <- fry_tests |>
  mutate(effect = -log10(p) * if_else(direction == "Up", 1, -1))
save_figure(
  module_tiles(
    rename(fry_signed, column = contrast), contrast_names, "signed\n-log10 p", NULL,
    "Modules as protein sets across the nine contrasts",
    sprintf(
      "limma::fry on the imputed matrix, protein design, subject block (correlation %.3f)",
      fry_correlation
    ),
    paste(
      "Fill: -log10 fry p, positive when the module moved up. fry's FDR runs over the twelve",
      "modules within each contrast."
    )
  ),
  "02_fry_contrasts", 9, 5.5
)
save_figure(
  module_tiles(
    rename(membership_contrast, column = contrast, effect = rho) |> mutate(fdr = NA_real_),
    contrast_names, "Spearman\nrho", c(-1, 1),
    "Module membership against protein significance",
    "Spearman between each member's kME and its protein-level moderated t, per module and contrast",
    paste(
      "Fill: rho; positive when the module's most central proteins moved up most. Members share",
      "a module, so the p describes rather than tests."
    )
  ),
  "03_membership_contrasts", 9, 5.5
)
save_figure(
  module_tiles(
    rename(membership_phenotype, column = outcome, effect = rho) |> mutate(fdr = NA_real_),
    unique(membership_phenotype$outcome), "Spearman\nrho", c(-1, 1),
    "Module membership against phenotype association",
    "Spearman between each member's kME and its protein rho with the phenotype over training",
    paste(
      "Fill: rho; positive when central proteins track the phenotype most positively.",
      "Descriptive, as above."
    )
  ),
  "04_membership_phenotype", 9, 5.5
)

primary <- members |>
  inner_join(filter(protein_t, contrast == "Training_Interaction"), by = "uniprot_id") |>
  filter(!is.na(t)) |>
  left_join(
    filter(membership_contrast, contrast == "Training_Interaction") |> select(module, rho),
    by = "module"
  ) |>
  mutate(
    module = factor(module, levels = module_order),
    panel = factor(sprintf("%s  rho %+.2f", module, rho), levels = unique(sprintf(
      "%s  rho %+.2f", module, rho
    )[order(module)]))
  )
save_figure(
  ggplot(primary, aes(kme, t)) +
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
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40")
    ),
  "05_membership_primary", 10, 8
)

packages <- c("here", "limma", "dplyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    eigengene_tests = eigengene_tests, fry_tests = fry_tests,
    membership_contrast = membership_contrast, membership_phenotype = membership_phenotype,
    correlation = c(eigengenes = me_correlation, imputed = fry_correlation),
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "module_tests.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    eigengene_tests = eigengene_tests, fry_tests = fry_tests,
    membership_contrast = membership_contrast, membership_phenotype = membership_phenotype,
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "04_test_modules.xlsx")
)
combined <- file.path(figure_dir, "04_test_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote module_tests.rds, 04_test_modules.xlsx and ", length(pages),
  " figures bundled into ", qpdf::pdf_length(combined), " pages"
)
