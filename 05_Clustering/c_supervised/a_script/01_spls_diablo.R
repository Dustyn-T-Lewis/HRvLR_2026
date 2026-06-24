# Exploratory supervised sPLS: baseline proteome -> comp_hypertrophy
# QUARANTINED. Supervised on n~15; at this sample size the model WILL
# overfit. Nested LOO-CV is the best available defense but remains
# unstable. Results are exploratory/non-confirmatory and MUST NOT be
# compared to the phenotype-blind WGCNA/Mfuzz engines or treated as
# confirmatory signal. See task description for framing rationale.

pacman::p_load(here, mixOmics, dplyr)
set.seed(42)

source(here("05_Clustering/functions/inputs.R"))
i <- load_clustering_inputs()

pheno <- read.csv(
  here("05_Clustering/d_integration/c_data/phenotype.csv")
)

# per-subject baseline (T1) proteome
t1_ids <- i$meta$sample_id[i$meta$timepoint == "T1"]
t1_subjects <- i$meta$subject[i$meta$timepoint == "T1"]

# subjects in both T1 proteome and phenotype (S29 has no T1 - drops)
keep <- intersect(t1_subjects, pheno$subject)
message(
  "subjects entering X: ", length(keep),
  " (dropped from pheno only: ",
  paste(setdiff(pheno$subject, keep), collapse = ", "), ")"
)

col_idx <- match(keep, t1_subjects)
x_mat <- t(i$abund[, t1_ids[col_idx], drop = FALSE])
rownames(x_mat) <- keep

y_vec <- pheno$comp_hypertrophy[match(keep, pheno$subject)]
names(y_vec) <- keep

n <- length(keep)
message("n = ", n, ", p = ", ncol(x_mat))

# tune keepX via LOO-CV
tune_out <- tryCatch(
  tune.spls(
    X = x_mat, Y = y_vec,
    ncomp = 2,
    test.keepX = c(10, 20, 30, 40),
    validation = "loo",
    measure = "MSE",
    progressBar = FALSE
  ),
  error = function(e) {
    message(
      "tune.spls failed: ", conditionMessage(e),
      " - falling back to keepX = 20"
    )
    NULL
  }
)

if (!is.null(tune_out)) {
  opt_keepX <- tune_out$choice.keepX
  message("tuned keepX: ", paste(opt_keepX, collapse = ", "))
} else {
  opt_keepX <- c(20, 20)
  message("fixed fallback keepX: 20 per component")
}

mod <- spls(
  X = x_mat, Y = y_vec,
  ncomp = 2,
  keepX = opt_keepX,
  mode = "regression",
  near.zero.var = TRUE
)

perf_out <- tryCatch(
  perf(mod, validation = "loo", progressBar = FALSE),
  error = function(e) {
    message("perf failed: ", conditionMessage(e))
    NULL
  }
)

if (!is.null(perf_out)) {
  q2_tab <- perf_out$measures$Q2$summary
  cv_df <- data.frame(
    component = q2_tab$comp,
    metric = "Q2",
    value = q2_tab$mean,
    sd = q2_tab$sd,
    note = "LOO-CV; n=15; expect poor/negative Q2 (overfitting)"
  )
  message("Q2 per component:")
  print(q2_tab)
} else {
  cv_df <- data.frame(
    component = 1:2,
    metric = "Q2",
    value = NA_real_,
    sd = NA_real_,
    note = "perf() failed at n=15; instability expected"
  )
}

dir.create(here("05_Clustering/c_supervised/c_data"),
  recursive = TRUE,
  showWarnings = FALSE
)
write.csv(cv_df,
  here("05_Clustering/c_supervised/c_data/cv_error.csv"),
  row.names = FALSE
)

get_loadings <- function(mod, comp) {
  ld <- mod$loadings$X[, comp]
  ld <- ld[ld != 0]
  data.frame(
    protein_id = names(ld),
    group_id = paste0("comp", comp),
    membership_weight = as.numeric(ld)
  )
}

membership <- bind_rows(
  get_loadings(mod, 1),
  get_loadings(mod, 2)
)
write.csv(membership,
  here("05_Clustering/c_supervised/c_data/membership.csv"),
  row.names = FALSE
)
message(
  "selected proteins: ", nrow(membership),
  " (comp1: ",
  sum(membership$group_id == "comp1"),
  ", comp2: ",
  sum(membership$group_id == "comp2"), ")"
)

dir.create(here("05_Clustering/c_supervised/b_reports"),
  recursive = TRUE,
  showWarnings = FALSE
)

pdf(here("05_Clustering/c_supervised/b_reports/spls_components.pdf"),
  width = 10, height = 8
)

# loading plot: comp 1
plotLoadings(mod,
  comp = 1, method = "median", contrib = "max",
  title = paste0(
    "sPLS comp 1 loadings - EXPLORATORY, n=",
    n, ", supervised, overfits, NOT confirmatory"
  )
)

# sample/variate plot
plotIndiv(mod,
  comp = c(1, 2),
  group = pheno$group_arm[match(keep, pheno$subject)],
  legend = TRUE,
  title = paste0(
    "sPLS variates (comp 1 vs 2) - EXPLORATORY, n=",
    n, ", NOT confirmatory"
  )
)

dev.off()

message("done; outputs written to 05_Clustering/c_supervised/")
