# Decomposes the module-eigengene -> Delta mCSA claim along the two axes that
# could be inflating it: whole-cohort imputation, and whole-cohort WGCNA.
#
#   matrix : missforest (1896 proteins, 12.2% of cells imputed)
#            complete   (the proteins with zero missing values, nothing imputed)
#   wgcna  : cohort     (fit once on all samples, as shipped)
#            foldrefit  (refit on the 15 training subjects inside every LOSO fold)
#
# WGCNA never sees the outcome, so the per-fold module fits are y-independent and
# are cached once and reused across the 200 permutations.

pacman::p_load(here, dplyr, tidyr, readr, tibble, glmnet, WGCNA, withr, parallel)

source(here("05_test", "prediction_shared", "_harness.R"))

WGCNA::disableWGCNAThreads()
OUT <- here("05_test", "complete_case_check", "c_data")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

N_PERM <- 200L
PERM_SEED <- 42L
CORES <- max(1L, parallel::detectCores() - 2L)

# Same knobs as F06/clustering.R::run_wgcna so the two are comparable.
fit_modules <- function(expr) {
  powers <- 1:20
  sft <- WGCNA::pickSoftThreshold(
    expr,
    powerVector = powers, networkType = "signed", verbose = 0
  )
  r2 <- -sign(sft$fitIndices$slope) * sft$fitIndices$SFT.R.sq
  hit <- which(r2 > 0.87)
  power <- if (length(hit)) powers[hit[1]] else if (nrow(expr) < 20) 18L else 16L

  net <- WGCNA::blockwiseModules(
    expr,
    networkType = "signed", power = power,
    minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
    randomSeed = 12345, verbose = 0
  )
  net$colors
}

# PC1 per module with the centre, scale and loadings kept, so a held-out sample
# can be projected through the training fit rather than re-estimated with it.
learn_eigengenes <- function(expr_train, colors) {
  mods <- setdiff(sort(unique(colors)), "grey")
  lapply(setNames(mods, mods), function(m) {
    x <- expr_train[, colors == m, drop = FALSE]
    ctr <- colMeans(x)
    scl <- apply(x, 2, stats::sd)
    scl[scl == 0 | !is.finite(scl)] <- 1
    z <- sweep(sweep(x, 2, ctr, "-"), 2, scl, "/")
    v <- svd(z, nu = 1, nv = 1)$v[, 1]
    # positive eigengene = higher mean module abundance
    if (stats::cor(as.numeric(z %*% v), rowMeans(z)) < 0) v <- -v
    list(genes = colnames(x), ctr = ctr, scl = scl, v = v)
  })
}

project_eigengenes <- function(fits, expr) {
  out <- vapply(fits, function(f) {
    z <- sweep(sweep(expr[, f$genes, drop = FALSE], 2, f$ctr, "-"), 2, f$scl, "/")
    as.numeric(z %*% f$v)
  }, numeric(nrow(expr)))
  rownames(out) <- rownames(expr)
  out
}

# Subject-level phase feature from a sample x module eigengene matrix.
phase_rows <- function(me, meta, phase, subjects) {
  tps <- switch(phase,
    baseline = "T1",
    training = c("T1", "T2"),
    acute = c("T2", "T3")
  )
  rows <- lapply(subjects, function(s) {
    ids <- vapply(tps, function(tp) {
      hit <- meta$sample[meta$subject == s & meta$timepoint == tp]
      if (length(hit) == 1) hit else NA_character_
    }, character(1))
    if (anyNA(ids) || !all(ids %in% rownames(me))) {
      return(NULL)
    }
    if (phase == "baseline") me[ids[[1]], ] else me[ids[[2]], ] - me[ids[[1]], ]
  })
  names(rows) <- subjects
  rows <- rows[!vapply(rows, is.null, logical(1))]
  x <- do.call(rbind, rows)
  rownames(x) <- names(rows)
  x
}

# Out-of-fold predictions where each fold carries its own feature matrix. The
# cohort variant passes the same matrix every fold; foldrefit passes its own.
loso_from_folds <- function(folds, y) {
  vapply(seq_along(folds), function(o) {
    sc <- scale_train_apply(folds[[o]]$x_train, folds[[o]]$x_test)
    fit_predict_glmnet(sc$train, y[-o], sc$test, "gaussian")
  }, numeric(1))
}

run_cell <- function(folds, y, label) {
  obs <- stat_q2(y, loso_from_folds(folds, y))
  pm <- withr::with_seed(PERM_SEED, replicate(N_PERM, sample.int(length(y))))
  null <- unlist(parallel::mclapply(seq_len(N_PERM), function(b) {
    yp <- y[pm[, b]]
    stat_q2(yp, loso_from_folds(folds, yp))
  }, mc.cores = CORES))
  null <- null[is.finite(null)]
  tibble(
    cell = label, n = length(y), n_modules = ncol(folds[[1]]$x_train),
    q2 = obs, perm_p = (sum(null >= obs) + 1) / (length(null) + 1),
    null_q2_mean = mean(null)
  )
}

norm_da <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
imp_da <- readRDS(here(
  "02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"
))

meta <- imp_da$metadata |>
  transmute(sample = Col_ID, subject = Subject_ID, timepoint = Timepoint)

complete <- rowSums(is.na(as.matrix(norm_da$data))) == 0
matrices <- list(
  missforest = t(as.matrix(imp_da$data)),
  complete   = t(as.matrix(norm_da$data)[complete, , drop = FALSE])
)
cat(sprintf(
  "missforest: %d proteins | complete-case: %d proteins (%.1f%% of the imputed set)\n",
  ncol(matrices$missforest), ncol(matrices$complete),
  100 * ncol(matrices$complete) / ncol(matrices$missforest)
))

pheno <- read_csv(here("04_Figures", "F06", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)
y_all <- setNames(pheno$d_mcsa, pheno$subject)

subjects <- sort(unique(meta$subject))

# One WGCNA fit per (matrix, held-out subject), reused across all three phases.
cat("fitting per-fold WGCNA (this is the slow part)\n")
refit <- lapply(matrices, function(expr) {
  cohort_colors <- fit_modules(expr)
  per_fold <- lapply(setNames(subjects, subjects), function(s) {
    keep <- rownames(expr) %in% meta$sample[meta$subject != s]
    colors <- fit_modules(expr[keep, , drop = FALSE])
    fits <- learn_eigengenes(expr[keep, , drop = FALSE], colors)
    project_eigengenes(fits, expr)
  })
  cohort_fits <- learn_eigengenes(expr, cohort_colors)
  list(cohort = project_eigengenes(cohort_fits, expr), per_fold = per_fold)
})

results <- list()
for (mat in names(matrices)) {
  for (wg in c("cohort", "foldrefit")) {
    for (phase in c("baseline", "training", "acute")) {
      me_cohort <- refit[[mat]]$cohort
      x_cohort <- phase_rows(me_cohort, meta, phase, subjects)
      subj <- intersect(rownames(x_cohort), names(y_all)[!is.na(y_all)])
      y <- as.numeric(y_all[subj])

      folds <- lapply(seq_along(subj), function(o) {
        me <- if (wg == "cohort") {
          me_cohort
        } else {
          refit[[mat]]$per_fold[[subj[[o]]]]
        }
        x <- phase_rows(me, meta, phase, subjects)[subj, , drop = FALSE]
        vary <- apply(x[-o, , drop = FALSE], 2, function(c) stats::sd(c) > 0)
        list(
          x_train = x[-o, vary, drop = FALSE],
          x_test  = x[o, vary, drop = FALSE]
        )
      })

      label <- paste(mat, wg, phase, sep = " | ")
      cat("running:", label, "\n")
      results[[label]] <- run_cell(folds, y, label)
    }
  }
}

out <- bind_rows(results) |>
  tidyr::separate(cell, c("matrix", "wgcna", "phase"), sep = " \\| ") |>
  arrange(desc(q2))

write_csv(out, file.path(OUT, "complete_case_metrics.csv"))
print(as.data.frame(out), digits = 4)
