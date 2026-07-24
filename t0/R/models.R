# Model specs, leave-one-subject-out resampling, the shared metric set, and the nested evaluator.
# Engines expose tuneable hyperparameters; ranger keeps impurity importance so vip is not needed.

# The three feature levels as three recipes over the same protein frame. protein keeps proteins;
# pathway and module swap in the fold-aware steps. subject is an id, .target the outcome.
build_level_recipes <- function(frame, gene_sets, gene_map, cfg) {
  base <- recipes::recipe(.target ~ ., data = frame) |>
    recipes::update_role(subject, new_role = "id") |>
    recipes::step_zv(recipes::all_predictors()) |>
    recipes::step_impute_median(recipes::all_predictors())
  list(
    protein = base |>
      recipes::step_nzv(recipes::all_predictors()) |>
      recipes::step_normalize(recipes::all_predictors()),
    pathway = base |>
      step_singscore(recipes::all_predictors(),
        gene_sets = gene_sets, gene_map = gene_map, min_size = cfg$gene_sets$min_size
      ) |>
      recipes::step_normalize(recipes::all_predictors()),
    module = base |>
      step_wgcna_me(recipes::all_predictors(), wgcna = cfg$wgcna) |>
      recipes::step_normalize(recipes::all_predictors())
  )
}

spec_glmnet <- function() {
  parsnip::logistic_reg(penalty = tune::tune(), mixture = tune::tune()) |>
    parsnip::set_engine("glmnet") |>
    parsnip::set_mode("classification")
}

spec_ranger <- function(cfg) {
  parsnip::rand_forest(mtry = tune::tune(), min_n = tune::tune(), trees = cfg$models$ranger$trees) |>
    parsnip::set_engine("ranger", importance = "impurity") |>
    parsnip::set_mode("classification")
}

loso_folds <- function(data, cfg) {
  rsample::group_vfold_cv(data,
    group = dplyr::all_of(cfg$cv$group),
    v = dplyr::n_distinct(data[[cfg$cv$group]])
  )
}

metric_set_t0 <- function() {
  yardstick::metric_set(yardstick::roc_auc, yardstick::bal_accuracy, yardstick::accuracy)
}

# Nested LOSO: tune on the inner resamples of each training fold, refit, predict the held-out subject.
# Every recipe transform (singscore, WGCNA) is re-prepped per fold, so no target sees its own network.
eval_nested_loso <- function(wflow, data, target, cfg) {
  outer <- loso_folds(data, cfg)
  mset <- metric_set_t0()
  preds <- purrr::map(outer$splits, function(split) {
    tr <- rsample::analysis(split)
    te <- rsample::assessment(split)
    inner <- rsample::vfold_cv(tr, v = min(cfg$cv$inner_v, nrow(tr) - 1))
    tuned <- tune::tune_grid(wflow,
      resamples = inner, grid = cfg$cv$tune_grid,
      metrics = mset, control = tune::control_grid(save_pred = FALSE)
    )
    best <- tune::select_best(tuned, metric = "roc_auc")
    fit <- tune::finalize_workflow(wflow, best) |> parsnip::fit(tr)
    dplyr::bind_cols(
      dplyr::select(te, dplyr::all_of(c(cfg$cv$group, target))),
      predict(fit, te, type = "prob"),
      predict(fit, te)
    )
  }) |> purrr::list_rbind()

  truth <- factor(preds[[target]])
  prob_col <- grep("^\\.pred_", names(preds), value = TRUE)[1]
  tibble::tibble(
    roc_auc = yardstick::roc_auc_vec(truth, preds[[prob_col]]),
    bal_accuracy = yardstick::bal_accuracy_vec(truth, preds$.pred_class),
    accuracy = yardstick::accuracy_vec(truth, preds$.pred_class)
  )
}

# Permutation null: shuffle the target within the modeling frame and re-run the whole nested LOSO.
# Seeds vary by index so the shuffles differ but the run stays reproducible.
perm_null <- function(wflow, data, target, cfg, n = cfg$permutation$reps) {
  base_seed <- cfg$seed
  purrr::map(seq_len(n), function(i) {
    withr::with_seed(base_seed + i, {
      shuffled <- data
      shuffled[[target]] <- sample(shuffled[[target]])
      eval_nested_loso(wflow, shuffled, target, cfg)
    }) |> dplyr::mutate(perm = i)
  }) |> purrr::list_rbind()
}

perm_pvalue <- function(observed, null_dist, metric = "roc_auc") {
  (sum(null_dist[[metric]] >= observed[[metric]]) + 1) / (nrow(null_dist) + 1)
}

# The full supervised sweep: every (contrast x target x feature level x engine) evaluated by nested
# LOSO against its own permutation null. The grid is data, so a new level or engine is one more row.
run_sweep <- function(proteins, meta, gene_sets, gmap, cfg) {
  purrr::map(cfg$scope$contrasts, function(cn) {
    delta <- build_delta(proteins, meta, cfg$design$contrasts[[cn]], cfg)
    purrr::map(cfg$scope$targets, function(tg) {
      frame <- make_modeling_frame(delta, meta, tg, cfg)
      recs <- build_level_recipes(frame, gene_sets, gmap, cfg)
      grid <- tidyr::expand_grid(
        level = intersect(cfg$scope$levels, names(recs)),
        model = c("glmnet", "ranger")
      )
      purrr::pmap(grid, function(level, model) {
        spec <- if (model == "glmnet") spec_glmnet() else spec_ranger(cfg)
        wf <- workflows::workflow(recs[[level]], spec)
        obs <- eval_nested_loso(wf, frame, ".target", cfg)
        null <- perm_null(wf, frame, ".target", cfg)
        tibble::tibble(
          contrast = cn, target = tg, level = level, model = model,
          roc_auc = obs$roc_auc, bal_accuracy = obs$bal_accuracy,
          p_perm = perm_pvalue(obs, null)
        )
      }) |> purrr::list_rbind()
    }) |> purrr::list_rbind()
  }) |>
    purrr::list_rbind() |>
    dplyr::mutate(
      wflow_id = paste(level, model, sep = "_"),
      q_perm = stats::p.adjust(p_perm, "BH")
    )
}

# Unsupervised counterpart: cluster each feature representation whole-cohort and ask whether the
# clusters recover either label (adjusted Rand index). Descriptive, so whole-cohort features are fine.
run_unsupervised <- function(proteins, meta, gene_sets, gmap, cfg) {
  labels <- meta |> dplyr::distinct(subject, responder, supplement)
  purrr::map(cfg$scope$contrasts, function(cn) {
    delta <- build_delta(proteins, meta, cfg$design$contrasts[[cn]], cfg)
    builders <- list(
      protein = function() delta,
      pathway = function() score_pathways(delta, proteins, gene_sets, cfg),
      module = function() module_eigengenes(delta, cfg)
    )
    levels <- lapply(builders[intersect(cfg$scope$levels, names(builders))], function(f) f())
    purrr::imap(levels, function(tbl, level) {
      mat <- tbl |>
        dplyr::inner_join(labels, "subject") |>
        dplyr::select(-subject, -responder, -supplement) |>
        as.matrix()
      mat <- scale(mat[, matrixStats::colVars(mat, na.rm = TRUE) > 0, drop = FALSE])
      mat[is.na(mat)] <- 0
      km <- withr::with_seed(cfg$seed, stats::kmeans(mat, centers = 2, nstart = 20))
      lab <- labels[match(tbl$subject, labels$subject), ]
      tibble::tibble(
        contrast = cn, level = level, k = 2,
        ari_responder = mclust_ari(km$cluster, lab$responder),
        ari_supplement = mclust_ari(km$cluster, lab$supplement)
      )
    }) |> purrr::list_rbind()
  }) |> purrr::list_rbind()
}

# Adjusted Rand index without a hard dependency on mclust.
mclust_ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  sum_comb <- function(x) sum(choose(x, 2))
  ai <- sum_comb(rowSums(tab))
  bi <- sum_comb(colSums(tab))
  index <- sum_comb(tab)
  expected <- ai * bi / choose(n, 2)
  max_index <- (ai + bi) / 2
  if (max_index == expected) {
    return(0)
  }
  (index - expected) / (max_index - expected)
}
