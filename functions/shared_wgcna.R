# The one WGCNA engine.
#
# Modules are DEFINED on within-subject-centred abundance and SCORED on raw abundance.
# That split is the whole point, and it is not cosmetic:
#
# WGCNA correlates proteins across samples treating them as independent. Our 45 samples are
# 16 subjects x 3 timepoints, so subject identity is a dominant variance driver - the FAQ's
# "strong driver that makes a subset of samples globally different from the rest". Built on
# raw abundance the modules substantially encode WHICH SUBJECT a sample came from rather
# than co-regulation: per-module subject ICC reaches 0.51, and the module assignment barely
# survives removing subject identity (adjusted Rand index 0.288). That artifact produced the
# earlier baseline-prediction q2 of 0.713, and it is why leave-one-subject-out destroyed it:
# LOSO removes exactly the thing the module was tracking.
#
# Centring each protein within subject removes that driver, so the modules capture how
# proteins move together WITHIN a subject - the training and acute response, which is the
# question. Eigengenes are then computed on raw abundance, so the between-subject HR/LR
# contrast is still available to test; it just no longer DEFINES the modules.
#
# Parameters follow the WGCNA FAQ's two recommended deviations from default (signed network,
# bicor) and are otherwise package defaults. The qmd's methods note carries the sources.
pacman::p_load(here, dplyr, tidyr, tibble)

# bicor is the FAQ's recommended robust alternative to Pearson. maxPOutliers = 0.05 is its
# explicit "strongly recommend" - the package default of 1 disables the outlier guard
# entirely. pearsonFallback = "individual" covers zero-MAD proteins, which missForest can
# produce by imputing tied values.
BICOR_ARGS <- list(
  corType = "bicor", maxPOutliers = 0.05, pearsonFallback = "individual"
)

centre_within_subject <- function(abund, subject) {
  abund - t(apply(t(abund), 2, function(x) ave(x, subject)))
}

# The soft power is the lowest whose signed scale-free fit clears RsquaredCut, which is what
# pickSoftThreshold returns at its documented default (0.85). The FAQ's sample-size table
# (signed, n > 40 -> 12) is a FALLBACK for when the fit never clears the threshold. Our fit
# clears it, so the table does not apply.
choose_power <- function(expr) {
  # pickSoftThreshold resolves corFnc by NAME from the search path, so WGCNA must be
  # attached, not merely namespaced.
  pacman::p_load(WGCNA)
  sft <- WGCNA::pickSoftThreshold(
    expr,
    powerVector = 1:20, networkType = "signed",
    corFnc = "bicor", corOptions = list(maxPOutliers = 0.05), verbose = 0
  )
  power <- sft$powerEstimate
  if (is.na(power)) power <- if (nrow(expr) > 40) 12L else 16L
  fit <- sft$fitIndices
  list(
    power = power, sft = sft,
    r2 = (-sign(fit$slope) * fit$SFT.R.sq)[power],
    mean_k = fit$mean.k.[power]
  )
}

fit_modules <- function(abund, subject) {
  pacman::p_load(WGCNA)
  WGCNA::disableWGCNAThreads()

  expr <- t(centre_within_subject(abund, subject))
  pw <- choose_power(expr)

  net <- do.call(WGCNA::blockwiseModules, c(
    list(
      expr,
      power = pw$power,
      networkType = "signed", TOMType = "signed",
      deepSplit = 2, minModuleSize = 30, mergeCutHeight = 0.15,
      pamStage = TRUE, pamRespectsDendro = TRUE,
      maxBlockSize = 2500, numericLabels = FALSE,
      randomSeed = 54321, verbose = 0
    ),
    BICOR_ARGS
  ))

  colors <- setNames(net$colors, colnames(expr))

  # Scored on RAW abundance with the centred-derived assignments, so the between-subject
  # HR/LR contrast survives to be tested. moduleEigengenes aligns each ME with its module's
  # average expression by default (align = "along average"); no hand-rolled sign fix.
  me <- WGCNA::moduleEigengenes(t(abund), colors = colors, excludeGrey = TRUE)$eigengenes
  rownames(me) <- colnames(abund)

  list(
    colors = colors, eigengenes = me, net = net,
    power = pw$power, r2 = pw$r2, mean_k = pw$mean_k, sft = pw$sft
  )
}

eigengene_long <- function(me, meta) {
  tibble::as_tibble(me, rownames = "sample_id") |>
    tidyr::pivot_longer(-sample_id, names_to = "module", values_to = "ME") |>
    dplyr::left_join(meta, by = "sample_id")
}

# Does subject identity still drive a module? Reported, not corrected: it is the diagnostic
# that justifies the centred definition. Raw R2 is inflated by 15 subject df on 45 samples
# (expectation 0.34 under pure noise), so the adjusted R2 and the ICC are the honest numbers.
subject_variance <- function(me_long) {
  pacman::p_load(lme4)
  me_long |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      adj_r2 = summary(stats::lm(ME ~ factor(subject)))$adj.r.squared,
      icc = {
        v <- as.data.frame(lme4::VarCorr(
          suppressMessages(lme4::lmer(ME ~ 1 + (1 | subject)))
        ))$vcov
        v[1] / (v[1] + v[2])
      },
      .groups = "drop"
    )
}
