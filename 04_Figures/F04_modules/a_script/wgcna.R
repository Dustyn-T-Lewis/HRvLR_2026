# Clustering compute for F04 - the figure builds its WGCNA modules and the
# module-phenotype linkages self-contained. WGCNA runs unsupervised on the full
# imputed proteome. Methods and seeds are preserved from the former standalone
# clustering stage; only print-only diagnostics (WGCNA power / mergeCut sweeps)
# and QC PDFs are dropped.
# Heavy stats packages are called namespace-qualified (not attached) so they do
# not collide in one process (e.g. WGCNA masking stats::cor).

pacman::p_load(here, dplyr, tidyr, readr, tibble)

# The trait table is built once at stage 00 (00_input/a_script/01_build_phenotype.R) and read
# here. MyoVision is not among them: its source column is a fibre COUNT, not an area.
ADAPTATION_TRAITS <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II",
  "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
)

# WGCNA runs unsupervised on the full imputed proteome, never on a DE-gated subset - filtering
# by differential expression before WGCNA collapses the network and invalidates the scale-free
# fit used to pick the soft power (Langfelder & Horvath, WGCNA FAQ).
load_clustering_inputs <- function() {
  imputed <- readRDS(here(
    "02_Normalization", "imputation", "c_data",
    "DAList_imputed_missforest.rds"
  ))
  stopifnot(identical(imputed$metadata$Col_ID, colnames(imputed$data)))

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
    ncol(imputed$data) == 45, !anyNA(imputed$data),
    identical(meta$sample_id, colnames(imputed$data))
  )

  list(full_abund = imputed$data, meta = meta)
}

# Unsupervised WGCNA on the full imputed proteome (no DE pre-filtering, per the
# WGCNA FAQ; mirrors YvO F04). Signed network; soft power = first with signed
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

  list(
    membership = membership, eigengene = eigengene,
    sft = sft, net = net, chosen_power = chosen_power,
    module_colors = net$colors
  )
}

# Per-module baseline prediction: does each module's baseline (T1) eigengene
# predict each training-adaptation trait? In-sample correlation + leave-one-out
# cross-validated Q^2 (Q2 > 0 = predicts out of sample), BH across the grid.
run_module_prediction <- function(wgcna_eig, pheno) {
  traits <- ADAPTATION_TRAITS
  t1 <- wgcna_eig[wgcna_eig$timepoint == "T1", c("subject", "group_id", "ME")]
  loo <- function(d) {
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

# Per-module phenotype LMM: each module eigengene modelled over time with a
# continuous trait as moderator and a subject random intercept,
#   ME ~ trait * timepoint + (1 | subject).
# The trait:timepoint term is the trajectory moderator (does the module's
# time-course bend with the size of the response); the trait term is the overall
# between-subject association. BH within trait, per term. This keeps the full
# repeated-measures structure that run_module_prediction (T1 only) and Panel C
# (T1->T3 delta) collapse away. Singular fits are flagged, not dropped: a
# time-invariant trait can drive the random-intercept variance to zero, which
# inflates the trait main effect, so the interaction is the term to trust.
run_module_phenotype_lmm <- function(wgcna_eig, pheno) {
  if (!requireNamespace("lmerTest", quietly = TRUE) ||
    !requireNamespace("lme4", quietly = TRUE)) {
    stop("run_module_phenotype_lmm needs the lmerTest and lme4 packages")
  }
  traits <- ADAPTATION_TRAITS
  modules <- setdiff(unique(wgcna_eig$group_id), "grey")

  fit_cell <- function(module, trait) {
    d <- merge(
      wgcna_eig[wgcna_eig$group_id == module, c("subject", "timepoint", "ME")],
      pheno[, c("subject", trait)]
    )
    names(d)[names(d) == trait] <- "trait"
    d$timepoint <- factor(d$timepoint)
    d <- d[stats::complete.cases(d[, c("ME", "trait")]), ]
    out <- data.frame(
      module = module, trait = trait, n = nrow(d),
      p_trait = NA_real_, p_trait_time = NA_real_, singular = NA
    )
    fit <- tryCatch(
      suppressWarnings(
        lmerTest::lmer(ME ~ trait * timepoint + (1 | subject), data = d)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(out)
    }
    a <- anova(fit)
    out$p_trait <- a["trait", "Pr(>F)"]
    out$p_trait_time <- a["trait:timepoint", "Pr(>F)"]
    out$singular <- lme4::isSingular(fit)
    out
  }

  grid <- expand.grid(module = modules, trait = traits, stringsAsFactors = FALSE)
  res <- do.call(rbind, Map(fit_cell, grid$module, grid$trait))
  res <- do.call(rbind, lapply(split(res, res$trait), function(s) {
    s$fdr_trait <- p.adjust(s$p_trait, "BH")
    s$fdr_trait_time <- p.adjust(s$p_trait_time, "BH")
    s
  }))
  rownames(res) <- NULL
  res
}
