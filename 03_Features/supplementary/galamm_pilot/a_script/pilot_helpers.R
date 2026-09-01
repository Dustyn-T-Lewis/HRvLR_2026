# Shared assembly and statistics for the galamm pilot (PREREG.md). Every
# script in the pilot reads its inputs through pilot_data(), so the
# complete-case gate and the blood-index join happen one way, once.

pilot_data <- function() {
  dal <- readRDS(here::here(
    "02_Normalization", "c_data", "DAList_normalized.rds"
  ))
  meta <- hlm_meta(dal)
  meta <- meta[match(colnames(dal$data), meta$sample), , drop = FALSE]
  blood <- blood_index_data()
  meta$blood_index <- blood$blood_index[match(meta$sample, blood$Col_ID)]
  stopifnot(!anyNA(meta$blood_index))
  mat <- dal$data[stats::complete.cases(dal$data), , drop = FALSE]
  stopifnot(nrow(mat) == 931L)
  list(mat = mat, meta = meta)
}

PHENO_ITEMS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)

# Six indicators z-scored across the 16 subjects, long format; the one
# d_1rm_ext NA drops its single row and nothing else.
pheno_long <- function() {
  pheno <- readr::read_csv(
    here::here("00_input", "c_data", "phenotype.csv"),
    show_col_types = FALSE
  )
  pheno |>
    dplyr::mutate(dplyr::across(
      dplyr::all_of(PHENO_ITEMS), \(x) as.numeric(scale(x))
    )) |>
    tidyr::pivot_longer(
      dplyr::all_of(PHENO_ITEMS),
      names_to = "item", values_to = "value"
    ) |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::mutate(item = factor(.data$item, levels = PHENO_ITEMS))
}

# Relative residual SDs per timepoint from a varIdent lme fit, T1 = 1.
varident_ratios <- function(fit) {
  rel <- coef(fit$modelStruct$varStruct,
    unconstrained = FALSE, allCoef = TRUE
  )
  c(T2 = unname(rel["T2"] / rel["T1"]), T3 = unname(rel["T3"] / rel["T1"]))
}

# Share of an indicator's variance the factor explains under a common
# residual variance: lambda^2 V / (lambda^2 V + sigma^2).
variance_explained <- function(lambda, var_eta, sigma2) {
  lambda^2 * var_eta / (lambda^2 * var_eta + sigma2)
}

# Empirical permutation p with the +1 correction, as in pi_permutation.R.
emp_p <- function(n_ge, n_perm) {
  (n_ge + 1) / (n_perm + 1)
}
