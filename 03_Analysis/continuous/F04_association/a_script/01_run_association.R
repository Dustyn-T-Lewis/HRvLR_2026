# Does the proteome track how much someone actually adapted? Two tiers.
#
# Tier 1 (cheap, no permutation, the actual screen): three feature levels
# (protein abundance, singscore, module eigengene) x three timepoint
# windows (T1 baseline, training = T2-T1, acute = T3-T2) x six phenotypes
# = 54 cells. Protein/pathway are fit with limma (vectorised across all
# features against the phenotype, no duplicateCorrelation needed since
# pred_contrast_matrix() already collapses every config to one row per
# subject); modules (p < n, only 12) get a plain per-eigengene lm. BH is
# applied within each of the 54 cells, never across them -- pooling would
# imply an independence the shared subjects and overlapping windows don't
# have.
#
# Tier 2 (LOSO + permutation, expensive): only for Tier-1 cells that clear
# a modest bar (>=1 feature at within-cell BH < 0.10). Full nested-LOSO
# elastic net with a subject-label-shuffle permutation null (B = 200),
# same harness as F04_classification. Reserving the expensive machinery
# for a small, named, reported subset is the concrete form of "don't
# permute by default."
#
# No BH across cells at either tier boundary, and no BH across Tier 2's
# promoted cells: same argument as F04_classification -- raw permutation p
# plus the exact promoted count is the honest alternative to a q that
# implies independence this design doesn't have. A lead needs
# perm_p < .05 AND Q2 > 0 (beats predicting the group mean).

pacman::p_load(here, dplyr, purrr, tidyr, tibble, limma, openxlsx)
source(here("functions", "pred_features.R"))
source(here("functions", "shared_prediction.R"))

OUT_DIR <- here("03_Analysis", "continuous", "F04_association", "c_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

LEVELS <- c("proteins", "singscore", "eigengenes")
CONFIGS <- c("T1", "training", "acute")
PHENOTYPES <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_mcsa",
  "d_1rm_legpress", "d_1rm_ext"
)
TIER1_BH_BAR <- 0.10

bundle <- pred_load("continuous")

# --- Tier 1 ------------------------------------------------------------

tier1_limma <- function(feat, y) {
  design <- stats::model.matrix(~y)
  fit <- limma::eBayes(limma::lmFit(t(feat), design))
  limma::topTable(fit, coef = 2, number = Inf, sort.by = "none") |>
    tibble::rownames_to_column("feature") |>
    transmute(feature, estimate = logFC, t, p = P.Value, bh = adj.P.Val)
}

tier1_lm <- function(feat, y) {
  purrr::map_dfr(colnames(feat), function(f) {
    fit <- stats::lm(feat[, f] ~ y)
    s <- summary(fit)$coefficients
    tibble(
      feature = f, estimate = s[2, 1], t = s[2, 3],
      p = s[2, 4]
    )
  }) |>
    mutate(bh = p.adjust(p, method = "BH"))
}

tier1_cell <- function(level, config, phenotype) {
  feat <- pred_contrast_matrix(
    bundle$feature_sets[[level]], bundle$meta, config
  )
  y_named <- pred_outcome(bundle, phenotype)
  al <- align_xy(feat, y_named)
  if (ncol(al$x) == 0 || nrow(al$x) < 6) {
    return(tibble())
  }
  res <- if (level == "eigengenes") {
    tier1_lm(al$x, al$y)
  } else {
    tier1_limma(al$x, al$y)
  }
  res |>
    mutate(
      level = level, config = config, phenotype = phenotype, n = nrow(al$x),
      .before = 1
    )
}

grid1 <- tidyr::expand_grid(
  level = LEVELS, config = CONFIGS, phenotype = PHENOTYPES
)
cat(sprintf("Tier 1: running %d cheap association cells...\n", nrow(grid1)))

tier1 <- purrr::pmap_dfr(grid1, tier1_cell)

tier1_summary <- tier1 |>
  summarise(
    n_tested = dplyr::n(), n_nominal = sum(p < 0.05),
    n_bh10 = sum(bh < TIER1_BH_BAR),
    .by = c(level, config, phenotype)
  )
print(as.data.frame(tier1_summary), digits = 3)

write.xlsx(
  list(tier1_cells = tier1, tier1_summary = tier1_summary),
  file.path(OUT_DIR, "F04_association_tier1.xlsx")
)

# --- Tier 2 --------------------------------------------------------------

promoted <- tier1_summary |> filter(n_bh10 >= 1)
cat(sprintf(
  "\nTier 2: %d of %d Tier-1 cells promoted (>=1 feature at BH < %.2f)\n",
  nrow(promoted), nrow(grid1), TIER1_BH_BAR
))

if (nrow(promoted) == 0) {
  cat("Nothing promoted. Tier 2 skipped: nothing to confirm.\n")
  tier2 <- tibble()
} else {
  print(as.data.frame(promoted), digits = 3)

  tier2_cell <- function(level, config, phenotype) {
    feat <- pred_contrast_matrix(
      bundle$feature_sets[[level]], bundle$meta, config
    )
    y_named <- pred_outcome(bundle, phenotype)
    al <- align_xy(feat, y_named)
    res <- suppressWarnings(run_cont_cell(al$x, al$y, "enet", phenotype))
    res$summary |>
      mutate(
        level = level, config = config, .before = 1,
        lead = perm_p_q2 < 0.05 & q2 > 0
      )
  }

  tier2 <- purrr::pmap_dfr(
    promoted |> select(level, config, phenotype), tier2_cell
  )
  print(as.data.frame(tier2 |> select(-outcome)), digits = 3)
  cat(sprintf(
    "\n%d of %d promoted cells are leads (perm_p_q2 < .05 and Q2 > 0)\n",
    sum(tier2$lead), nrow(tier2)
  ))
}

write.xlsx(
  list(
    tier1_cells = tier1, tier1_summary = tier1_summary,
    tier2_promoted = promoted, tier2_cells = tier2
  ),
  file.path(OUT_DIR, "F04_association_source_data.xlsx")
)
