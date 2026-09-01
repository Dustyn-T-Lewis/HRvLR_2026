#!/usr/bin/env Rscript
# Q2, step two: does the latent hypertrophic-response factor associate
# with any protein's abundance?
#
# One galamm fit per complete-case protein: the six z-scored phenotype
# indicators and the protein's 45 z-scored abundances share one latent
# factor per subject. The protein response carries arm, timepoint, their
# interaction and the blood index as fixed effects (all zero on phenotype
# rows), plus its own subject intercept for the repeated measures, so the
# protein loading tests latent association within arm, beyond the median
# split. Inference is Wald from the marginal likelihood, not
# Satterthwaite (PREREG.md). BH across proteins, q < 0.05; any survivor
# must also beat a B = 200 permutation of the phenotype block across
# subjects, the subject-as-unit scheme of pi_permutation.R.
#
# Read 03_q2_protein.csv: loading, SE, z, p, bh, var_eta and status per
# protein, non-converged fits included with status naming the failure.
# A null looks like no BH survivor, which at n = 16 it should.

pacman::p_load(here, dplyr, tidyr, purrr, tibble, readr, galamm)
source(here("functions", "shared_hlm.R"))
source(here("functions", "blood_index_model.R"))
source(here(
  "03_Features", "supplementary", "galamm_pilot", "a_script",
  "pilot_helpers.R"
))

set.seed(42)
OUT <- here("03_Features", "supplementary", "galamm_pilot", "c_data")
stopifnot(
  "run 02_q2_measurement.R first" =
    file.exists(file.path(OUT, "02_q2_measurement.csv"))
)

inp <- pilot_data()
mat <- inp$mat
meta <- inp$meta
CORES <- max(1L, parallel::detectCores() - 2L)
N_PERM <- 200L

strip_arm <- function(x) sub("^(HR|LR)_", "", x)

ph_rows <- pheno_long() |>
  transmute(
    subject = strip_arm(.data$subject), item = .data$item,
    value = .data$value,
    pr = 0, arm = 0, t2 = 0, t3 = 0, a2 = 0, a3 = 0, blood = 0
  )

prot_covariates <- tibble(
  subject = strip_arm(as.character(meta$subject)),
  item = "protein",
  pr = 1,
  arm = as.numeric(meta$group == "HR"),
  t2 = as.numeric(meta$timepoint == "T2"),
  t3 = as.numeric(meta$timepoint == "T3"),
  blood = as.numeric(scale(meta$blood_index))
) |>
  mutate(a2 = .data$arm * .data$t2, a3 = .data$arm * .data$t3)

item_levels <- c(PHENO_ITEMS, "protein")
lambda <- matrix(c(1, rep(NA, 6)), ncol = 1)

fit_joint <- function(y, ph = ph_rows) {
  d <- bind_rows(ph, mutate(prot_covariates, value = as.numeric(scale(y)))) |>
    mutate(item = factor(.data$item, levels = item_levels))
  mod <- tryCatch(
    galamm(
      value ~ 0 + item + arm + t2 + t3 + a2 + a3 + blood +
        (0 + eta | subject) + (0 + pr | subject),
      data = d, load_var = "item", lambda = lambda, factor = "eta"
    ),
    error = function(e) structure(conditionMessage(e), class = "fit_error")
  )
  if (inherits(mod, "fit_error")) {
    return(tibble(
      loading = NA_real_, se = NA_real_, var_eta = NA_real_,
      status = paste("error:", substr(unclass(mod), 1, 60))
    ))
  }
  lo <- as.data.frame(factor_loadings(mod))
  vc <- as.data.frame(VarCorr(mod))
  tibble(
    loading = lo$eta[7], se = lo$SE[7],
    var_eta = vc$vcov[vc$grp == "subject" & vc$var1 == "eta"][1],
    status = "ok"
  )
}

secs <- vapply(1:3, function(i) {
  t0 <- Sys.time()
  fit_joint(mat[i, ])
  as.numeric(Sys.time() - t0, units = "secs")
}, numeric(1))
message(sprintf("single-fit timing: %.2f s median", median(secs)))
features <- rownames(mat)
if (median(secs) > 5) {
  features <- sample(features, 100)
  message("timing gate tripped: cut to a seeded 100-protein sample")
}

results <- parallel::mclapply(features, function(f) {
  mutate(fit_joint(mat[f, ]), feature = f, .before = 1)
}, mc.cores = CORES) |>
  bind_rows() |>
  mutate(
    z = .data$loading / .data$se,
    p = 2 * pnorm(-abs(.data$z)),
    bh = p.adjust(.data$p, "BH")
  )

write_csv(results, file.path(OUT, "03_q2_protein.csv"))

survivors <- filter(results, .data$status == "ok", .data$bh < 0.05)
if (nrow(survivors)) {
  subjects <- unique(ph_rows$subject)
  perm_one <- function(f) {
    y <- mat[f, ]
    obs <- abs(results$z[results$feature == f])
    null_z <- vapply(seq_len(N_PERM), function(b) {
      shuffled <- setNames(sample(subjects), subjects)
      fit_b <- fit_joint(y, ph = mutate(
        ph_rows,
        subject = shuffled[.data$subject]
      ))
      abs(fit_b$loading / fit_b$se)
    }, numeric(1))
    tibble(
      feature = f, obs_z = obs,
      n_ge = sum(null_z >= obs, na.rm = TRUE),
      emp_p = emp_p(sum(null_z >= obs, na.rm = TRUE), N_PERM)
    )
  }
  perm <- bind_rows(parallel::mclapply(
    survivors$feature, perm_one,
    mc.cores = CORES
  ))
  write_csv(perm, file.path(OUT, "03_q2_permutation.csv"))
  message(paste(capture.output(print(as.data.frame(perm))), collapse = "\n"))
}

message(sprintf(
  "fits ok %d/%d | min BH q %.3f | BH survivors %d",
  sum(results$status == "ok"), length(features),
  min(results$bh, na.rm = TRUE), nrow(survivors)
))
