#!/usr/bin/env Rscript
# Q2, step one: can a latent hypertrophic-response factor be estimated
# from the six phenotype indicators at all, and is it just d_mcsa?
#
# The HR/LR label is a median split of a continuous response, and six
# noisy indicators plausibly measure one trait. galamm fits the
# one-factor measurement model on the z-scored indicators (95 rows:
# 16 subjects x 6 items minus the one d_1rm_ext NA), anchor
# comp_hypertrophy fixed to 1. A per-item residual variance
# (dispformula = ~ (1 | item)) is attempted first and the homoscedastic
# model kept if it fails, with the attempt recorded either way.
#
# Read 02_q2_measurement.csv for loadings, Wald SEs and the share of each
# indicator's variance the factor explains, and 02_q2_scores.csv for the
# empirical-Bayes factor scores beside d_mcsa. The PREREG stop rules:
# no convergence or all free loadings |lambda|/SE < 2 closes Q2 as
# unidentifiable; Spearman rho > 0.9 against d_mcsa closes it as "the
# factor is d_mcsa". A null here is a factor nobody can estimate.

pacman::p_load(here, dplyr, tidyr, readr, tibble, galamm)
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here("03_Features", "04_galamm_pilot", "a_script", "pilot_helpers.R"))

set.seed(42)
OUT <- here("03_Features", "04_galamm_pilot", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

long <- pheno_long()
lambda <- matrix(c(1, rep(NA, 5)), ncol = 1)

fit_measurement <- function(disp) {
  args <- list(
    formula = value ~ 0 + item + (0 + eta | subject), data = long,
    load_var = "item", lambda = lambda, factor = "eta"
  )
  if (disp) args$dispformula <- ~ (1 | item)
  tryCatch(do.call(galamm, args), error = function(e) NULL)
}

mod <- fit_measurement(disp = TRUE)
disp_used <- !is.null(mod)
if (!disp_used) mod <- fit_measurement(disp = FALSE)
stopifnot("measurement model did not converge" = !is.null(mod))
message(
  "per-item dispersion: ",
  if (disp_used) "converged" else "failed at start values, common sigma kept"
)

load_tbl <- as.data.frame(factor_loadings(mod))
vc <- as.data.frame(VarCorr(mod))
var_eta <- vc$vcov[vc$grp == "subject"]
sigma2 <- vc$vcov[vc$grp == "Residual"]

loadings <- tibble(
  item = PHENO_ITEMS,
  loading = load_tbl$eta,
  se = load_tbl$SE,
  z = load_tbl$eta / load_tbl$SE,
  var_explained = variance_explained(load_tbl$eta, var_eta, sigma2)
)

pheno <- read_csv(
  here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)
scores <- tibble(
  subject = rownames(ranef(mod)$subject),
  eta = ranef(mod)$subject$eta
) |>
  left_join(
    pheno |> select("subject", "group_arm", "d_mcsa", "comp_hypertrophy"),
    by = "subject"
  )
rho_mcsa <- cor(scores$eta, scores$d_mcsa, method = "spearman")
rho_comp <- cor(scores$eta, scores$comp_hypertrophy, method = "spearman")

write_csv(loadings, file.path(OUT, "02_q2_measurement.csv"))
write_csv(
  mutate(scores, rho_d_mcsa = rho_mcsa, rho_comp = rho_comp),
  file.path(OUT, "02_q2_scores.csv")
)

n_identified <- sum(abs(loadings$z) >= 2, na.rm = TRUE)
message(paste(capture.output(print(as.data.frame(loadings))), collapse = "\n"))
message(sprintf(
  "free loadings with |z| >= 2: %d of 5 | rho %.3f vs d_mcsa, %.3f vs comp",
  n_identified, rho_mcsa, rho_comp
))
message(sprintf(
  "PREREG gate: %s",
  if (n_identified == 0) {
    "Q2 closes as unidentifiable"
  } else if (rho_mcsa > 0.9) {
    "rho > 0.9, the factor is d_mcsa; protein fits are confirmation only"
  } else {
    "measurement model passes, proceed to per-protein fits"
  }
))
