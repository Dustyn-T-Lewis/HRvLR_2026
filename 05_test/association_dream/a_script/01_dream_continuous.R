# Continuous-phenotype dream fits: feature ~ phenotype * timepoint + ACCUM_VL +
# (1|subject), one fit per phenotype on both the protein and the singscore matrix.
# The phenotype main effect is the baseline (T1) association; the phenotype x
# timepoint terms are the trajectory-tracks-phenotype signal. BH is applied within
# each matrix x phenotype x coefficient block.
if (!exists("scores")) source(here::here("05_test", "association_dream", "a_script", "setup.R"))
pacman::p_load(variancePartition, BiocParallel)
set.seed(42)

phenotype_coefs <- function(trait) {
  c(trait, sprintf("%s:timepointT2", trait), sprintf("%s:timepointT3", trait))
}

fit_continuous <- function(feature_mat, matrix_label, trait) {
  keep <- !is.na(meta$ACCUM_VL) & !is.na(meta[[trait]])
  design <- meta[keep, ]
  design$accum_vl_z <- as.numeric(scale(design$ACCUM_VL))
  design[[trait]] <- as.numeric(scale(design[[trait]]))
  form <- as.formula(sprintf(
    "~ %s * timepoint + accum_vl_z + (1 | subject)", trait
  ))
  fit <- dream(
    feature_mat[, keep], form, design,
    BPPARAM = SerialParam(), suppressWarnings = TRUE, quiet = TRUE
  )
  fit <- eBayes(fit)

  bind_rows(lapply(phenotype_coefs(trait), function(cf) {
    tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
    tibble(
      matrix = matrix_label,
      feature = rownames(tt),
      phenotype = trait,
      coefficient = dplyr::recode(
        cf,
        !!trait := "phenotype",
        !!sprintf("%s:timepointT2", trait) := "phenotype:T2",
        !!sprintf("%s:timepointT3", trait) := "phenotype:T3"
      ),
      estimate = tt$logFC, t = tt$t, p_value = tt$P.Value
    )
  }))
}

matrices <- list(protein = expr, pathway = scores)

cont_res <- bind_rows(lapply(names(matrices), function(ml) {
  bind_rows(lapply(names(cont_traits), function(tr) {
    fit_continuous(matrices[[ml]], ml, tr)
  }))
}))

cont_res <- cont_res |>
  group_by(matrix, phenotype, coefficient) |>
  mutate(fdr = p.adjust(p_value, "BH")) |>
  ungroup() |>
  arrange(p_value)

write_csv(cont_res, file.path(dat, "dream_continuous.csv"))

message(sprintf(
  "continuous dream: %d rows; FDR<0.05 hits: %d (proteins %d, pathways %d)",
  nrow(cont_res), sum(cont_res$fdr < 0.05),
  sum(cont_res$fdr < 0.05 & cont_res$matrix == "protein"),
  sum(cont_res$fdr < 0.05 & cont_res$matrix == "pathway")
))
