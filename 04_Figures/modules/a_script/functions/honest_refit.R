pacman::p_load(dplyr, tibble, purrr)

jaccard <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

# Scaled PC1 loading, sign-oriented so the score tracks the row mean.
fit_pc1 <- function(train_mat) {
  ctr <- colMeans(train_mat)
  scl <- apply(train_mat, 2, sd)
  scl[scl == 0 | !is.finite(scl)] <- 1
  z <- sweep(sweep(train_mat, 2, ctr, "-"), 2, scl, "/")
  v <- svd(z, nu = 0, nv = 1)$v[, 1]
  score <- as.numeric(z %*% v)
  if (cor(score, rowMeans(z)) < 0) v <- -v
  list(center = ctr, scale = scl, loading = v)
}

project_pc1 <- function(fit, new_mat) {
  z <- sweep(sweep(new_mat, 2, fit$center, "-"), 2, fit$scale, "/")
  as.numeric(z %*% fit$loading)
}

match_modules <- function(colors_full, colors_train) {
  full_levels <- setdiff(unique(colors_full), "grey")
  train_levels <- setdiff(unique(colors_train), "grey")
  purrr::map_dfr(full_levels, function(f) {
    prot_f <- names(colors_full)[colors_full == f]
    j <- vapply(train_levels, function(t) {
      jaccard(prot_f, names(colors_train)[colors_train == t])
    }, numeric(1))
    tibble(full = f, train = train_levels[which.max(j)], jaccard = max(j))
  })
}

run_honest_refit <- function(abund, sample_meta, outcome_named, phase, power = NULL) {
  pacman::p_load(WGCNA)
  WGCNA::disableWGCNAThreads()
  full <- fit_modules(abund, sample_meta)
  if (is.null(power)) power <- full$chosen_power
  subjects <- unique(sample_meta$subject)

  fold <- function(hold) {
    train_ids <- sample_meta$sample_id[sample_meta$subject != hold]
    expr_tr <- t(abund[, train_ids])
    set.seed(42)
    net <- WGCNA::blockwiseModules(
      expr_tr,
      networkType = "signed", TOMType = "signed", power = power,
      minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
      pamRespectsDendro = FALSE, saveTOMs = FALSE, randomSeed = 42, verbose = 0
    )
    ct <- WGCNA::labels2colors(net$colors)
    names(ct) <- colnames(expr_tr)
    mm <- match_modules(full$colors, ct) |> filter(jaccard >= 0.5)
    hold_ids <- sample_meta$sample_id[sample_meta$subject == hold]
    purrr::map_dfr(seq_len(nrow(mm)), function(i) {
      prot <- names(full$colors)[full$colors == mm$full[i]]
      fitp <- fit_pc1(t(abund[prot, train_ids]))
      me_hold <- project_pc1(fitp, t(abund[prot, hold_ids, drop = FALSE]))
      tibble(
        module = paste0("ME", mm$full[i]), subject = hold,
        sample_id = hold_ids, me = me_hold
      )
    })
  }

  held <- purrr::map_dfr(subjects, fold)
  me_phase <- .phase_eigengene(
    held |> tidyr::pivot_wider(id_cols = sample_id, names_from = module, values_from = me) |>
      tibble::column_to_rownames("sample_id"),
    sample_meta, phase
  )

  ins <- .phase_eigengene(full$eigengenes, sample_meta, phase)
  modules <- intersect(unique(ins$module), unique(me_phase$module))
  purrr::map_dfr(modules, function(m) {
    xi <- ins |> filter(module == m)
    yi <- outcome_named[xi$subject]
    xo <- me_phase |> filter(module == m)
    yo <- outcome_named[xo$subject]
    tibble(
      module = m,
      auc_insample = module_auc(yi, xi$value),
      auc_loso = module_auc(yo, xo$value)
    ) |>
      mutate(drop = auc_insample - auc_loso)
  })
}
