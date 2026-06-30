# Clustering compute for F06 - the figure builds its WGCNA modules and the
# module-phenotype linkages self-contained. WGCNA runs unsupervised on the full
# imputed proteome. Methods and seeds are preserved from the former standalone
# clustering stage; only print-only diagnostics (WGCNA power / mergeCut sweeps)
# and QC PDFs are dropped.
# Heavy stats packages are called namespace-qualified (not attached) so they do
# not collide in one process (e.g. WGCNA masking stats::cor).

pacman::p_load(here, dplyr, tidyr, readr, tibble)

# pi-gated input set (Xiao pi-value < 0.05 in any contrast) from the canonical
# non-imputed limma fit; abundance from the missForest-imputed matrix; gap =
# HR-LR logFC at T1/T2/T3.
load_clustering_inputs <- function() {
  primary <- readRDS(here("03_DEP", "a_non_imputed", "c_data", "01_limma_DAList.rds"))
  imputed <- readRDS(here(
    "02_Normalization", "imputation", "c_data",
    "DAList_imputed_missforest.rds"
  ))
  stopifnot(identical(imputed$metadata$Col_ID, colnames(imputed$data)))

  pi_mat <- vapply(
    primary$results, `[[`, numeric(nrow(primary$data)), "pi_score"
  )
  rownames(pi_mat) <- rownames(primary$data)
  pi_keep <- apply(pi_mat, 1, function(x) any(x < 0.05, na.rm = TRUE))
  pi_set <- rownames(pi_mat)[pi_keep]
  pi_set <- intersect(pi_set, rownames(imputed$data))

  full_abund <- imputed$data

  md <- imputed$metadata
  meta <- data.frame(
    sample_id = md$Col_ID,
    subject = md$Subject_ID,
    group_arm = md$Group,
    timepoint = md$Timepoint,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  stopifnot(
    ncol(full_abund) == 45, !anyNA(full_abund),
    identical(meta$sample_id, colnames(full_abund))
  )

  list(full_abund = full_abund, meta = meta, pi_set = pi_set)
}

# Per-subject phenotype table: T2 composite hypertrophy + the six T1->T2 deltas.
build_phenotype_table <- function(meta_path) {
  meta <- readr::read_csv(meta_path, show_col_types = FALSE)

  trait_cols <- c(
    "fCSA_Type_I_Pre", "fCSA_Type_II_Pre", "MyoVision_fCSA_Type_I__Pre",
    "mCSA_Pre", "X1RM_Leg_Pre", "X1RM._Ext_Pre"
  )
  delta_names <- c(
    "d_fcsa_I", "d_fcsa_II", "d_myovision_fcsa_I",
    "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
  )

  arm <- meta |>
    distinct(Subject_ID, Group) |>
    rename(subject = Subject_ID, group_arm = Group)

  comp <- meta |>
    filter(Timepoint == "T2") |>
    transmute(
      subject = Subject_ID,
      comp_hypertrophy = readr::parse_number(COMP.HYPERTROPHY)
    )

  t1 <- meta |>
    filter(Timepoint == "T1") |>
    select(subject = Subject_ID, all_of(trait_cols))
  t2 <- meta |>
    filter(Timepoint == "T2") |>
    select(subject = Subject_ID, all_of(trait_cols))

  deltas <- t2 |>
    left_join(t1, by = "subject", suffix = c("_t2", "_t1")) |>
    mutate(
      d_fcsa_I = fCSA_Type_I_Pre_t2 - fCSA_Type_I_Pre_t1,
      d_fcsa_II = fCSA_Type_II_Pre_t2 - fCSA_Type_II_Pre_t1,
      d_myovision_fcsa_I = MyoVision_fCSA_Type_I__Pre_t2 -
        MyoVision_fCSA_Type_I__Pre_t1,
      d_mcsa = mCSA_Pre_t2 - mCSA_Pre_t1,
      d_1rm_legpress = X1RM_Leg_Pre_t2 - X1RM_Leg_Pre_t1,
      d_1rm_ext = `X1RM._Ext_Pre_t2` - `X1RM._Ext_Pre_t1`
    ) |>
    select(subject, all_of(delta_names))

  arm |>
    left_join(comp, by = "subject") |>
    left_join(deltas, by = "subject")
}

# Unsupervised WGCNA on the full imputed proteome (no DE pre-filtering, per the
# WGCNA FAQ; mirrors YvO F06). Signed network; soft power = first with signed
# scale-free R^2 > 0.87, else the sample-size-derived signed fallback (Horvath
# FAQ table: <20->18, 20-30->16, 30-40->14, >40->12). The downstream linkages
# (run_module_prediction, run_module_responder) reduce each subject to one row
# per timepoint, so the paired design needs no random effect.
# Deterministic: single-threaded, fixed blockwiseModules randomSeed.
run_wgcna <- function(abund, meta) {
  # WGCNA must be attached: moduleEigengenes resolves cor() via do.call from the
  # search path, so namespace-only loading would hit stats::cor. This masks
  # stats::cor for the session; downstream cor() calls are stats::-qualified.
  pacman::p_load(WGCNA)
  WGCNA::disableWGCNAThreads()

  expr <- t(abund)
  gsg <- WGCNA::goodSamplesGenes(expr, verbose = 0)
  if (!gsg$allOK) expr <- expr[gsg$goodSamples, gsg$goodGenes]
  meta <- meta[match(rownames(expr), meta$sample_id), ]

  powers <- 1:20
  sft <- WGCNA::pickSoftThreshold(
    expr,
    powerVector = powers, networkType = "signed", verbose = 0
  )
  fit <- sft$fitIndices
  r2 <- -sign(fit$slope) * fit$SFT.R.sq
  n_samp <- nrow(expr)
  fallback_power <- if (n_samp < 20) 18L else if (n_samp <= 30) 16L else if (n_samp <= 40) 14L else 12L
  candidates <- which(r2 > 0.87)
  chosen_power <- if (length(candidates)) powers[candidates[1]] else fallback_power

  net <- WGCNA::blockwiseModules(
    expr,
    networkType = "signed", power = chosen_power,
    minModuleSize = 30, mergeCutHeight = 0.25, deepSplit = 2,
    randomSeed = 12345, verbose = 0
  )

  mes_all <- WGCNA::moduleEigengenes(expr, colors = net$colors)$eigengenes
  stopifnot(rownames(mes_all) == meta$sample_id)
  grey_col <- grep("^MEgrey$", names(mes_all), value = TRUE)
  mes_real <- mes_all[, setdiff(names(mes_all), grey_col), drop = FALSE]

  # sign-fix: positive ME -> higher mean abundance for the module
  for (mod_col in names(mes_real)) {
    in_module <- net$colors == sub("^ME", "", mod_col)
    mod_mean <- rowMeans(expr[, in_module, drop = FALSE])
    if (stats::cor(mes_real[[mod_col]], mod_mean) < 0) {
      mes_real[[mod_col]] <- -mes_real[[mod_col]]
    }
  }
  module_labels <- sub("^ME", "", names(mes_real))

  membership <- data.frame(
    protein_id = colnames(expr), group_id = net$colors,
    membership_weight = 1L, stringsAsFactors = FALSE
  )
  eigengene <- do.call(rbind, lapply(seq_along(module_labels), function(j) {
    data.frame(
      sample_id = meta$sample_id, subject = meta$subject,
      group_arm = meta$group_arm, timepoint = meta$timepoint,
      group_id = module_labels[j], ME = mes_real[[j]],
      stringsAsFactors = FALSE
    )
  }))

  list(membership = membership, eigengene = eigengene)
}

# Per-module baseline prediction: does each module's baseline (T1) eigengene
# predict each training-adaptation trait? In-sample correlation + leave-one-out
# cross-validated Q^2 (Q2 > 0 = predicts out of sample), BH across the grid.
run_module_prediction <- function(wgcna_eig, pheno) {
  traits <- c(
    "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_myovision_fcsa_I",
    "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
  )
  t1 <- wgcna_eig[wgcna_eig$timepoint == "T1", c("subject", "group_id", "ME")]
  loo <- function(d) {
    d <- d[stats::complete.cases(d), ]
    if (nrow(d) < 6) {
      return(NA_real_)
    }
    p <- vapply(seq_len(nrow(d)), function(i) predict(lm(y ~ x, d[-i, ]), d[i, ]), numeric(1))
    1 - sum((d$y - p)^2) / sum((d$y - mean(d$y))^2)
  }
  rows <- list()
  for (m in unique(t1$group_id)) {
    for (tr in traits) {
      d <- merge(t1[t1$group_id == m, c("subject", "ME")], pheno[, c("subject", tr)])
      d <- stats::setNames(d[, c("ME", tr)], c("x", "y"))
      d <- d[stats::complete.cases(d), ]
      ct <- suppressWarnings(stats::cor.test(d$x, d$y))
      rows[[length(rows) + 1L]] <- data.frame(
        module = m, trait = tr, n = nrow(d),
        r = unname(ct$estimate), p = ct$p.value, q2 = loo(d)
      )
    }
  }
  res <- do.call(rbind, rows)
  res$p_bh <- p.adjust(res$p, "BH")
  res
}

# Per-module responder signal: point-biserial correlation of each module
# eigengene with responder group (HR = 1) at each timepoint.
run_module_responder <- function(wgcna_eig) {
  d0 <- wgcna_eig
  d0$g <- as.integer(d0$group_arm == "HR")
  parts <- split(d0, list(d0$group_id, d0$timepoint), drop = TRUE)
  res <- do.call(rbind, lapply(parts, function(s) {
    ct <- suppressWarnings(stats::cor.test(s$ME, s$g))
    data.frame(
      module = s$group_id[1], timepoint = s$timepoint[1],
      r = unname(ct$estimate), p = ct$p.value
    )
  }))
  res$p_bh <- p.adjust(res$p, "BH")
  res
}
