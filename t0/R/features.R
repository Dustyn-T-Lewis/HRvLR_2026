# Feature construction at three levels. Each returns a subjects x features tibble keyed by subject,
# so any level plugs into the same modeling frame. Deltas require a complete timepoint pair; subjects
# missing an endpoint drop out of that contrast rather than being imputed.

# Join a protein-level feature table to one per-subject target, dropping subjects without it.
make_modeling_frame <- function(delta_tbl, meta, target, cfg) {
  targets <- meta |>
    dplyr::distinct(subject, .data[[target]]) |>
    dplyr::rename(.target = dplyr::all_of(target))
  delta_tbl |>
    dplyr::inner_join(targets, by = "subject") |>
    dplyr::filter(!is.na(.target)) |>
    dplyr::relocate(subject, .target)
}

build_delta <- function(proteins, meta, contrast, cfg) {
  expr <- proteins$expr
  wide <- meta |>
    dplyr::select(sample, subject, timepoint) |>
    dplyr::filter(sample %in% colnames(expr))

  pick <- function(tp) {
    s <- wide$sample[wide$timepoint == tp]
    m <- expr[, s, drop = FALSE]
    colnames(m) <- wide$subject[match(s, wide$sample)]
    m
  }

  if (identical(contrast$type, "baseline")) {
    mat <- pick(contrast$at)
  } else {
    a <- pick(contrast$from)
    b <- pick(contrast$to)
    common <- intersect(colnames(a), colnames(b))
    mat <- b[, common, drop = FALSE] - a[, common, drop = FALSE]
  }
  tibble::as_tibble(t(mat), rownames = "subject")
}

# singscore is scored per sample against a fixed gene set, so it never learns from other subjects -
# leakage-safe by construction. Kept here for QC and reused inside step_singscore().
score_pathways <- function(feature_tbl, proteins, gene_sets, cfg) {
  mat <- feature_tbl |>
    dplyr::select(subject, dplyr::where(is.numeric)) |>
    tibble::column_to_rownames("subject") |>
    as.matrix() |>
    t()
  key <- cfg$ids$protein_key
  rownames(mat) <- proteins$annotation[[cfg$ids$gene_key]][
    match(rownames(mat), proteins$annotation[[key]])
  ]
  mat <- mat[!is.na(rownames(mat)) & !duplicated(rownames(mat)), , drop = FALSE]
  mat[] <- t(apply(mat, 1, \(r) {
    r[is.na(r)] <- stats::median(r, na.rm = TRUE)
    r
  }))

  ranked <- singscore::rankGenes(mat)
  scores <- vapply(gene_sets, function(genes) {
    genes <- intersect(genes, rownames(mat))
    if (length(genes) < cfg$gene_sets$min_size) {
      return(rep(NA_real_, ncol(mat)))
    }
    singscore::simpleScore(ranked, upSet = genes)$TotalScore
  }, numeric(ncol(mat)))

  tibble::as_tibble(scores, rownames = NULL) |>
    dplyr::mutate(subject = colnames(mat), .before = 1) |>
    dplyr::select(subject, dplyr::where(~ !anyNA(.x)))
}

# WGCNA modules learned on one set of subjects, projected onto another by stored PC loadings.
# Standalone version for QC; the fold-aware form lives in step_wgcna_me().
module_eigengenes <- function(feature_tbl, cfg) {
  mods <- fit_wgcna_modules(feature_tbl, cfg$wgcna)
  project_eigengenes(mods, feature_tbl)
}

fit_wgcna_modules <- function(feature_tbl, wcfg) {
  # blockwiseModules calls an unqualified cor(); without WGCNA attached it resolves to stats::cor
  # and errors on WGCNA's extra args. Bind WGCNA::cor for the call, restore on exit.
  old_cor <- mget("cor", envir = .GlobalEnv, ifnotfound = list(NULL))$cor
  assign("cor", WGCNA::cor, envir = .GlobalEnv)
  on.exit(
    if (is.null(old_cor)) {
      suppressWarnings(rm(list = "cor", envir = .GlobalEnv))
    } else {
      assign("cor", old_cor, envir = .GlobalEnv)
    },
    add = TRUE
  )

  datExpr <- feature_tbl |>
    tibble::column_to_rownames("subject") |>
    as.matrix()
  datExpr <- datExpr[, matrixStats::colVars(datExpr) > 0, drop = FALSE]

  power <- wcfg$power
  if (is.null(power)) {
    fit <- WGCNA::pickSoftThreshold(datExpr,
      powerVector = wcfg$power_candidates, RsquaredCut = wcfg$fit_cut,
      networkType = wcfg$network_type, verbose = 0
    )
    power <- fit$powerEstimate
    if (is.na(power)) power <- max(wcfg$power_candidates)
  }
  net <- WGCNA::blockwiseModules(datExpr,
    power = power, networkType = wcfg$network_type,
    minModuleSize = wcfg$min_module_size,
    mergeCutHeight = wcfg$merge_cut_height,
    numericLabels = TRUE, verbose = 0
  )
  keep <- net$colors != 0
  loadings <- tapply(colnames(datExpr)[keep], net$colors[keep], function(px) {
    stats::prcomp(datExpr[, px, drop = FALSE], center = TRUE, scale. = TRUE, rank. = 1)
  })
  list(loadings = loadings, power = power)
}

project_eigengenes <- function(mods, feature_tbl) {
  datExpr <- feature_tbl |>
    tibble::column_to_rownames("subject") |>
    as.matrix()
  me <- vapply(mods$loadings, function(pc) {
    px <- rownames(pc$rotation)
    stats::predict(pc, datExpr[, px, drop = FALSE])[, 1]
  }, numeric(nrow(datExpr)))
  me <- matrix(me,
    nrow = nrow(datExpr),
    dimnames = list(rownames(datExpr), paste0("ME", names(mods$loadings)))
  )
  tibble::as_tibble(me) |>
    dplyr::mutate(subject = rownames(datExpr), .before = 1)
}
