# Categorical dream fit: feature ~ group * timepoint + ACCUM_VL + (1|subject) on
# both matrices. The group term is the HR-vs-LR baseline (T1) difference; the
# group x timepoint terms are the divergence in trajectory. Composite hypertrophy
# is deliberately absent, since the groups were defined from it. BH within each
# matrix x coefficient block.
if (!exists("scores")) source(here::here("05_test", "association_dream", "a_script", "setup.R"))
pacman::p_load(variancePartition, BiocParallel)
set.seed(42)

fit_categorical <- function(feature_mat, matrix_label) {
  keep <- !is.na(meta$ACCUM_VL)
  design <- meta[keep, ]
  design$accum_vl_z <- as.numeric(scale(design$ACCUM_VL))
  fit <- dream(
    feature_mat[, keep],
    ~ group * timepoint + accum_vl_z + (1 | subject),
    design,
    BPPARAM = SerialParam(), suppressWarnings = TRUE, quiet = TRUE
  )
  fit <- eBayes(fit)

  coefs <- c(
    groupHR = "group",
    `groupHR:timepointT2` = "group:T2",
    `groupHR:timepointT3` = "group:T3"
  )
  bind_rows(lapply(names(coefs), function(cf) {
    tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
    tibble(
      matrix = matrix_label, feature = rownames(tt),
      coefficient = coefs[[cf]],
      estimate = tt$logFC, t = tt$t, p_value = tt$P.Value
    )
  }))
}

matrices <- list(protein = expr, pathway = scores)

cat_res <- bind_rows(lapply(names(matrices), function(ml) {
  fit_categorical(matrices[[ml]], ml)
})) |>
  group_by(matrix, coefficient) |>
  mutate(fdr = p.adjust(p_value, "BH")) |>
  ungroup() |>
  arrange(p_value)

write_csv(cat_res, file.path(dat, "dream_categorical.csv"))

message(sprintf(
  "categorical dream: %d rows; FDR<0.05 hits: %d (proteins %d, pathways %d)",
  nrow(cat_res), sum(cat_res$fdr < 0.05),
  sum(cat_res$fdr < 0.05 & cat_res$matrix == "protein"),
  sum(cat_res$fdr < 0.05 & cat_res$matrix == "pathway")
))
