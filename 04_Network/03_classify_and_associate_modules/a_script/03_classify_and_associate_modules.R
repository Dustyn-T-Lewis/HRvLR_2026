# The module eigengenes put through the same three questions the proteins and sets answered: the
# nine contrasts, the eight classification tasks, and association with the ten phenotypes.
#
# The contrasts reuse the protein fit's design and subject block, with the within-subject
# correlation re-estimated on the eigengenes, so a module is tested exactly as a protein was.
# Twelve features make BH a far weaker filter than 1,900 proteins; each table reports the count
# chance would give beside it.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(limma)
})

stage <- here("04_Network", "03_classify_and_associate_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  modules = "04_Network/01_build_modules/c_data/modules.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_build_modules first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
modules <- readRDS(paths[["modules"]])
d <- readRDS(paths[["design"]])
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

me <- modules$eigengenes
design <- d$dal$design$design_matrix
contrast_matrix <- d$dal$design$contrast_matrix
subject <- d$dal$metadata$subject
stopifnot(identical(colnames(me), rownames(design)))


# ---- contrasts -----------------------------------------------------------------------------

correlation <- duplicateCorrelation(me, design, block = subject)$consensus
module_fit <- lmFit(me, design, block = subject, correlation = correlation) |>
  contrasts.fit(contrast_matrix) |>
  eBayes(robust = TRUE)
module_contrasts <- map(colnames(contrast_matrix), \(contrast) {
  topTable(module_fit, coef = contrast, number = Inf, sort.by = "none") |>
    rownames_to_column("module") |>
    transmute(contrast, module, logFC, t, p = P.Value, fdr = adj.P.Val)
}) |>
  list_rbind()
message("eigengene within-subject correlation: ", round(correlation, 3))


# ---- the eight classification tasks --------------------------------------------------------

meta <- modules$meta
arm_of <- deframe(distinct(meta, subject, arm))
by_subject <- function(mat, timepoint) {
  rows <- filter(meta, timepoint == .env$timepoint)
  set_names(as.data.frame(mat[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
}
change <- function(mat, from, to) {
  a <- by_subject(mat, from)
  b <- by_subject(mat, to)
  both <- intersect(colnames(a), colnames(b))
  b[, both, drop = FALSE] - a[, both, drop = FALSE]
}
windows <- list(
  training = change(me, "T1", "T2"), baseline = by_subject(me, "T1"), acute = change(me, "T2", "T3")
)
within_arm <- function(from, to, arm) {
  a <- by_subject(me, from)
  b <- by_subject(me, to)
  both <- intersect(colnames(a), colnames(b))
  both <- both[arm_of[both] == arm]
  list(positive = b[, both, drop = FALSE], negative = a[, both, drop = FALSE], paired = TRUE)
}
between_arm <- function(values) {
  list(
    positive = values[, arm_of[colnames(values)] == "HR", drop = FALSE],
    negative = values[, arm_of[colnames(values)] == "LR", drop = FALSE], paired = FALSE
  )
}
tasks <- list(
  Training_HR = within_arm("T1", "T2", "HR"), Training_LR = within_arm("T1", "T2", "LR"),
  Acute_HR = within_arm("T2", "T3", "HR"), Acute_LR = within_arm("T2", "T3", "LR"),
  Baseline_HRvLR = between_arm(windows$baseline),
  Trained_HRvLR = between_arm(by_subject(me, "T2")),
  Training_change_HRvLR = between_arm(windows$training),
  Acute_change_HRvLR = between_arm(windows$acute)
)
# direction = "<" pins pROC's orientation, so AUC above 0.5 means higher in the later timepoint
# or in HR.
module_auc <- imap(tasks, \(spec, name) {
  tibble(module = rownames(me), task = name, paired = spec$paired) |>
    mutate(
      auc = map_dbl(module, \(m) {
        as.numeric(pROC::auc(pROC::roc(
          controls = spec$negative[m, ], cases = spec$positive[m, ],
          direction = "<", quiet = TRUE
        )))
      }),
      p = map_dbl(module, \(m) {
        wilcox.test(spec$positive[m, ], spec$negative[m, ], paired = spec$paired)$p.value
      })
    )
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = task)


# ---- phenotype association -----------------------------------------------------------------

outcomes <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_fcsa_mixed", "d_nfibre_mixed",
  "d_nfibre_I", "d_mcsa", "d_1rm_legpress", "d_1rm_ext", "volume_load"
)
outcome_of <- function(name, subjects) deframe(select(phenotype, subject, all_of(name)))[subjects]
# The t approximation: above nine observations cor.test's "exact" Spearman p is an Edgeworth
# series that returns 0 in the far tail.
spearman_by_row <- function(values, outcome) {
  usable <- !is.na(outcome)
  fits <- apply(values[, usable, drop = FALSE], 1, \(row) {
    test <- cor.test(row, outcome[usable], method = "spearman", exact = FALSE)
    c(r = unname(test$estimate), p = test$p.value)
  })
  tibble(module = rownames(values), n = sum(usable), r = fits["r", ], p = fits["p", ])
}
module_association <- imap(windows, \(values, window) {
  map(set_names(outcomes), \(name) spearman_by_row(values, outcome_of(name, colnames(values)))) |>
    list_rbind(names_to = "outcome") |>
    mutate(window = window, .before = 1)
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(window, outcome))

chance_expectation <- bind_rows(
  summarise(module_contrasts,
    tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05), .by = contrast
  ) |>
    transmute(analysis = "contrast", comparison = contrast, tested, nominal, fdr_sig),
  summarise(module_auc, tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05),
    .by = task
  ) |>
    transmute(analysis = "classification", comparison = task, tested, nominal, fdr_sig),
  summarise(module_association,
    tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05), .by = c(window, outcome)
  ) |>
    transmute(analysis = paste0("association: ", window), comparison = outcome, tested,
      nominal, fdr_sig
    )
) |>
  mutate(expected = 0.05 * tested, .after = tested)
print(as.data.frame(
  chance_expectation |>
    summarise(across(c(tested, expected, nominal, fdr_sig), sum), .by = analysis)
))


# ---- figures -------------------------------------------------------------------------------

save_figure <- function(figure, name, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 200, bg = "white"
    )
  })
}
module_order <- rev(modules$module_summary$module)
# One tile per module and column. Every module is drawn, since there are twelve; the label marks
# nominal p and a black border BH < 0.05.
module_tiles <- function(data, columns, fill_label, midpoint, limits, title, subtitle, caption) {
  data |>
    mutate(
      module = factor(module, levels = module_order), column = factor(column, levels = columns)
    ) |>
    ggplot(aes(column, module, fill = effect)) +
    geom_tile(colour = "white") +
    geom_tile(data = \(x) filter(x, fdr < 0.05), colour = "black", linewidth = 0.8, fill = NA) +
    geom_text(aes(label = if_else(p < 0.05, sprintf("%.3f", p), "")), size = 2.6) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = midpoint, limits = limits
    ) +
    labs(
      x = NULL, y = NULL, fill = fill_label, title = title, subtitle = subtitle,
      caption = paste(
        caption, "Label: nominal p where below 0.05. Black border: BH < 0.05.",
        "Table: c_data/03_classify_and_associate_modules.xlsx."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)
    )
}

save_figure(
  module_tiles(
    rename(module_contrasts, column = contrast, effect = t), colnames(contrast_matrix),
    "moderated t", 0, NULL, "Module eigengenes across the nine contrasts",
    sprintf("limma, protein design, subject block (correlation %.3f), BH within contrast",
      correlation
    ),
    "Fill: moderated t of each eigengene in each contrast."
  ),
  "01_contrasts", 9, 5.5
)
save_figure(
  module_tiles(
    rename(module_auc, column = task, effect = auc), names(tasks), "AUC", 0.5, c(0, 1),
    "Module classification", "pROC AUC, Wilcoxon p (signed-rank within arm), BH within task",
    "Fill: AUC, above 0.5 higher at the later timepoint or in HR."
  ),
  "02_classification", 9, 5.5
)
window_titles <- c(
  training = "training change (T2 - T1)", baseline = "level at T1",
  acute = "acute change (T3 - T2)"
)
iwalk(window_titles, \(title, window) {
  save_figure(
    module_tiles(
      module_association |>
        filter(window == .env$window) |>
        rename(column = outcome, effect = r),
      outcomes, "rho", 0, c(-1, 1), paste("Module association,", title),
      "Spearman, subjects pooled, BH within outcome",
      "Fill: Spearman rho between eigengene and phenotype."
    ),
    sprintf("03_association_%02d_%s", match(window, names(window_titles)), window), 9, 5.5
  )
})

packages <- c("here", "limma", "pROC", "dplyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    contrasts = module_contrasts, auc = module_auc, association = module_association,
    chance_expectation = chance_expectation, correlation = correlation,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "module_results.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    chance_expectation = chance_expectation, contrasts = module_contrasts,
    classification = module_auc, association = arrange(module_association, p),
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "03_classify_and_associate_modules.xlsx")
)
combined <- file.path(figure_dir, "03_classify_and_associate_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote module_results.rds, 03_classify_and_associate_modules.xlsx and a ",
  length(pages), "-page figure PDF"
)
