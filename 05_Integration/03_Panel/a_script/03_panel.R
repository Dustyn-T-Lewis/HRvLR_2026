# Exploratory multi-feature panels: can an elastic-net model tell HR from LR in subjects it did not
# train on? Five inputs (T1, T2 and T3 levels, training and acute change) by two feature sets
# (module eigengenes, and proteins observed in every subject of the input). Each panel is scored
# by repeated stratified 5-fold CV and judged against 1,000 label permutations of the whole
# pipeline. Nothing here is validated outside this cohort.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
})
source(here("R", "helpers.R"))

stage <- here("05_Integration", "03_Panel")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  limpa = "01_Preprocess/03_Quantify/c_data/DAList_limpa.rds",
  modules = "04_Modules/01_Build/c_data/modules.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
n_permutations <- 1000
cores <- min(8, parallel::detectCores())

# Exploratory elastic-net panel: repeated stratified k-fold CV around nestedcv's glmnet, alpha 0.5
# fixed and lambda.1se from an inner leave-one-out. Stratified outer folds keep the HR:LR balance
# of every training set; pooled leave-one-out does not, and an empty model then scores AUC 0
# (Airola et al. 2011). x has one row per subject. Returns each repeat's out-of-fold AUC, how
# often each feature was kept across the repeat-by-fold models, and the share that kept none.
panel_cv <- function(x, y, repeats = 10, k = 5, seed = 1) {
  inner <- length(y) - ceiling(length(y) / k)
  runs <- map(seq_len(repeats), \(r) {
    set.seed(seed + r)
    fit <- suppressWarnings(nestedcv::nestcv.glmnet(y, x,
      family = "binomial", outer_folds = caret::createFolds(y, k = k), n_inner_folds = inner,
      alphaSet = 0.5, min_1se = 1, finalCV = FALSE
    ))
    list(
      auc = pROC::roc(fit$output$testy, fit$output$predyp,
        levels = levels(y), direction = "<", quiet = TRUE
      ) |>
        pROC::auc() |>
        as.numeric(),
      kept = map(fit$outer_result, \(fold) setdiff(names(fold$coef), "(Intercept)"))
    )
  })
  kept <- list_flatten(map(runs, "kept"))
  list(
    auc = map_dbl(runs, "auc"),
    frequency = table(factor(unlist(kept), colnames(x))) / length(kept),
    empty = mean(lengths(kept) == 0)
  )
}

dal <- readRDS(paths[["limpa"]])
meta <- as_tibble(dal$metadata)
modules <- readRDS(paths[["modules"]])
stopifnot(identical(colnames(modules$eigengenes), meta$sample_id))
feature_label <- c(
  with(dal$annotation, set_names(gene, uniprot_id)),
  set_names(rownames(modules$eigengenes), rownames(modules$eigengenes))
)
arm_of <- set_names(meta$arm, meta$subject)

views <- function(values) {
  list(
    T1 = subject_values(values, meta, "T1"),
    T2 = subject_values(values, meta, "T2"),
    T3 = subject_values(values, meta, "T3"),
    training = subject_change(values, meta, "T1", "T2"),
    acute = subject_change(values, meta, "T2", "T3")
  )
}
# The blood correction here sees time point only, never arm, so the labels cannot inform it. The
# eigengenes are rescored on that matrix with the module assignments from 01_Build. glmnet
# takes no missing value; keeping only proteins observed in every subject of the input avoids
# imputing outside the folds. The modules were built on all 16 subjects, a leak the protein
# panels do not share.
label_free <- blood_corrected(dal$data, meta$blood, model.matrix(~ 0 + timepoint, meta))
colours <- set_names(modules$membership$module, modules$membership$uniprot_id)[rownames(label_free)]
eigengenes <- WGCNA::moduleEigengenes(t(label_free), colors = colours, excludeGrey = TRUE)
eigengenes <- t(as.matrix(eigengenes$eigengenes))
dimnames(eigengenes) <- list(sub("^ME", "", rownames(eigengenes)), colnames(label_free))
matrices <- list(
  modules = views(eigengenes),
  proteins = views(label_free)
) |>
  map_depth(2, t)

panels <- expand_grid(features = names(matrices), input = names(matrices$modules)) |>
  mutate(result = map2(features, input, \(features, input) {
    x <- matrices[[features]][[input]]
    y <- factor(arm_of[rownames(x)], c("LR", "HR"))
    message(features, " ", input, ": ", nrow(x), " subjects, ", ncol(x), " features")
    observed <- panel_cv(x, y)
    null <- parallel::mclapply(seq_len(n_permutations), \(i) {
      set.seed(1e4 + i)
      mean(panel_cv(x, sample(y), seed = i)$auc)
    }, mc.cores = cores) |>
      unlist()
    stopifnot(length(null) == n_permutations, !anyNA(null))
    list(x = x, y = y, observed = observed, null = null)
  }))

summary_table <- panels |>
  transmute(
    features, input,
    n_subjects = map_int(result, \(r) length(r$y)),
    n_hr = map_int(result, \(r) sum(r$y == "HR")),
    n_lr = map_int(result, \(r) sum(r$y == "LR")),
    n_features = map_int(result, \(r) ncol(r$x)),
    auc_mean = map_dbl(result, \(r) mean(r$observed$auc)),
    auc_min = map_dbl(result, \(r) min(r$observed$auc)),
    auc_max = map_dbl(result, \(r) max(r$observed$auc)),
    empty_models = map_dbl(result, \(r) r$observed$empty),
    null_median = map_dbl(result, \(r) median(r$null)),
    null_95 = map_dbl(result, \(r) unname(quantile(r$null, 0.95))),
    permutation_p = map_dbl(result, \(r) permutation_p(mean(r$observed$auc), r$null))
  )
selection <- panels |>
  mutate(frequency = map(result, \(r) enframe(c(r$observed$frequency), "feature", "frequency"))) |>
  select(features, input, frequency) |>
  unnest(frequency) |>
  filter(frequency > 0) |>
  mutate(label = feature_label[feature], .after = feature) |>
  arrange(features, input, desc(frequency))

null_page <- panels |>
  transmute(features, input, null = map(result, "null")) |>
  unnest(null) |>
  ggplot(aes(null)) +
  geom_histogram(binwidth = 0.025, boundary = 0, fill = "grey70") +
  geom_vline(aes(xintercept = auc_mean), summary_table, colour = "firebrick") +
  geom_text(
    aes(-Inf, Inf, label = sprintf("AUC %.2f\npermutation p %.3f", auc_mean, permutation_p)),
    summary_table,
    hjust = -0.1, vjust = 1.2, size = 2.6
  ) +
  facet_grid(features ~ factor(input, names(matrices$modules))) +
  labs(
    title = "Exploratory panels against their permutation nulls",
    x = "mean cross-validated AUC", y = "permutations",
    caption = stringr::str_wrap(width = 170, paste(
      "Elastic net (alpha 0.5, lambda.1se from inner leave-one-out) scored by 10 repeats of",
      "stratified 5-fold CV. Grey: mean AUC over the same pipeline with HR and LR labels",
      sprintf("shuffled, %d times; red: observed.", n_permutations),
      "p = (b + 1) / (m + 1). The null sits below 0.5: when a model keeps no feature its scores",
      "follow each training fold's class mix, so compare the red line with the grey, not with 0.5.",
      "Exploratory; no external validation. The module panels carry a leak:",
      "modules were built on all 16 subjects. Table: c_data/03_panel.xlsx, panels."
    ))
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.caption = element_text(hjust = 0))
selection_page <- selection |>
  slice_max(frequency, n = 15, by = c(features, input), with_ties = FALSE) |>
  mutate(key = paste(label, input, features, sep = "___")) |>
  ggplot(aes(frequency, reorder(key, frequency))) +
  geom_col(fill = "grey40") +
  scale_y_discrete(labels = \(x) sub("___.*", "", x)) +
  facet_wrap(~ paste(features, input), scales = "free_y", ncol = 5) +
  labs(
    title = "How often each feature entered the model",
    x = "share of the 50 repeat-by-fold models", y = NULL,
    caption = paste(
      "Up to 15 features per panel. A feature kept in most models is the panel; one kept in a few",
      "is noise. Table: c_data/03_panel.xlsx, selection."
    )
  ) +
  theme_minimal(base_size = 8) +
  theme(plot.caption = element_text(hjust = 0))

pdf(file.path(figures, "03_panel_figures.pdf"), width = 11, height = 8.5)
print(null_page)
if (nrow(selection) > 0) print(selection_page)
invisible(dev.off())

write_workbook(
  list(panels = summary_table, selection = selection),
  c(
    paste(
      "Per panel: subjects, features, mean and range of the 10 repeat AUCs, share of models that",
      "kept no feature, permutation null median and 95th percentile, permutation p."
    ),
    "Per panel and feature: share of the 50 repeat-by-fold models that kept it."
  ),
  file.path(out, "03_panel.xlsx"),
  inputs, paths
)
