# Does the proteome split the subjects on its own, with no label in sight?
#
# Model-based clustering selects the component count by BIC and can return one,
# which is the answer that says there is nothing to find (Scrucca et al. 2016).
# BIC alone is not enough here. At 15 baseline subjects a mixture model will
# report clusters in data that has none, which is the failure Senbabaoglu et
# al. (2014, Sci Rep 4:6207) documented for consensus clustering and which
# applies to any stability or fit criterion read without a null.
#
# So every cell is calibrated: draw from a single multivariate Gaussian with
# the observed mean and covariance, refit, and ask how often that structureless
# data produces a BIC gain as large as the observed one. Same n, same
# dimension, same correlation structure, no clusters.
#
# Dimension is reduced by PCA before fitting because 12 modules against 15
# subjects leaves an EEE covariance alone costing 78 parameters. The component
# count is swept over PC_GRID rather than fixed, and every value is reported;
# picking the k that gave the best answer is the failure mode this sweep
# exists to prevent.

# MASS is called qualified, never attached: MASS::select() masks dplyr::select()
# for every script sourced after this one in the same session.
pacman::p_load(here, dplyr, tidyr, purrr, tibble, readr, mclust, openxlsx)

source(here("functions", "feature_levels.R"))

OUT_DIR <- here("03_Features", "03_Subtypes", "c_data")

PC_GRID <- 2:4
G_MAX <- 4
N_NULL <- 999
N_PROTEINS <- 500
GATE_ALPHA <- 0.05

sample_meta <- function() {
  dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
  as.data.frame(dal$metadata) |>
    transmute(
      sample_id = Col_ID, subject = Subject_ID, timepoint = Timepoint
    )
}

meta <- sample_meta()
baseline <- meta |> filter(.data$timepoint == "T1")

eigengene_baseline <- function() {
  wide <- read_csv(
    here("03_Features", "02_WGCNA", "c_data", "wgcna_eigengene.csv"),
    show_col_types = FALSE
  ) |>
    filter(.data$sample_id %in% baseline$sample_id) |>
    pivot_wider(names_from = "group_id", values_from = "ME") |>
    left_join(baseline, by = "sample_id")
  mat <- as.matrix(dplyr::select(wide, -sample_id, -subject, -timepoint))
  rownames(mat) <- wide$subject
  mat
}

# The high-dimensional sensitivity space. Ranked by variance because an
# unsupervised filter cannot use a label to choose features without leaking it.
protein_baseline <- function(n_top = N_PROTEINS) {
  mat <- protein_matrix()[, baseline$sample_id, drop = FALSE]
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  keep <- order(matrixStats::rowVars(mat), decreasing = TRUE)[seq_len(n_top)]
  out <- t(mat[keep, , drop = FALSE])
  rownames(out) <- baseline$subject[match(rownames(out), baseline$sample_id)]
  out
}

# The statistic the null is built for: how much better BIC does the selected
# component count do than a single component. Returns NA when no model in the
# family is estimable, which is a legitimate outcome at this sample size.
bic_gain <- function(z) {
  fit <- try(mclust::Mclust(z, G = 1:G_MAX, verbose = FALSE), silent = TRUE)
  if (inherits(fit, "try-error") || is.null(fit)) {
    return(list(g = NA_integer_, gain = NA_real_, fit = NULL))
  }
  best <- apply(fit$BIC, 1, function(r) {
    if (all(is.na(r))) NA_real_ else max(r, na.rm = TRUE)
  })
  list(g = fit$G, gain = unname(best[fit$G] - best[1]), fit = fit)
}

cluster_cell <- function(space, mat, k) {
  scores <- stats::prcomp(mat, center = TRUE, scale. = TRUE)$x[, seq_len(k),
    drop = FALSE
  ]
  obs <- bic_gain(scores)
  null <- vapply(
    seq_len(N_NULL),
    function(i) {
      draw <- MASS::mvrnorm(nrow(scores), colMeans(scores), stats::cov(scores))
      unlist(bic_gain(draw)[c("g", "gain")])
    },
    numeric(2)
  )
  usable <- !is.na(null["gain", ])
  list(
    summary = tibble(
      space = space, n_pc = k, n_subjects = nrow(scores),
      best_g = obs$g, bic_gain = obs$gain,
      null_g_above_1 = mean(null["g", usable] > 1),
      null_gain_median = stats::median(null["gain", usable]),
      null_gain_q95 = stats::quantile(null["gain", usable], 0.95),
      p_empirical =
        (sum(null["gain", usable] >= obs$gain) + 1) / (sum(usable) + 1)
    ),
    draws = tibble(
      space = space, n_pc = k,
      null_g = null["g", usable], null_gain = null["gain", usable]
    )
  )
}

set.seed(42)
spaces <- list(eigengenes = eigengene_baseline(), proteins = protein_baseline())

fitted <- purrr::map(names(spaces), function(sp) {
  purrr::map(PC_GRID, function(k) cluster_cell(sp, spaces[[sp]], k))
}) |> purrr::list_flatten()

cells <- purrr::map_dfr(fitted, "summary")
null_draws <- purrr::map_dfr(fitted, "draws")

gate_open <- any(cells$p_empirical < GATE_ALPHA, na.rm = TRUE)

# Reported whatever the gate does, and labelled by it. A forced two-component
# split scored against the labels is informative context when the gate is shut
# and would be the headline if it opened; it is never promoted silently.
two_group <- purrr::map_dfr(names(spaces), function(sp) {
  scores <- stats::prcomp(spaces[[sp]], center = TRUE, scale. = TRUE)$x[, 1:2]
  fit <- mclust::Mclust(scores, G = 2, verbose = FALSE)
  assign_tbl <- tibble(
    subject = rownames(spaces[[sp]]),
    cluster = paste0("c", fit$classification)
  )
  labels <- read_csv(
    here(
      "03_Features", "01_Responsiveness", "c_data",
      "02_candidate_labels.csv"
    ),
    show_col_types = FALSE
  )
  purrr::map_dfr(unique(labels$label), function(l) {
    j <- inner_join(
      assign_tbl, filter(labels, .data$label == l),
      by = "subject"
    )
    tibble(
      space = sp, label = l, n = nrow(j),
      ari = mclust::adjustedRandIndex(j$cluster, j$level)
    )
  })
})

verdict <- if (gate_open) {
  "gate open: clustering beat its null"
} else {
  "descriptive only: gate shut, no structure to interpret"
}
two_group$gate_open <- gate_open
two_group$interpretation <- verdict

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(null_draws, file.path(OUT_DIR, "01_null_draws.csv"))
write.xlsx(
  list(cluster_cells = cells, forced_two_group = two_group),
  file.path(OUT_DIR, "01_subtypes.xlsx")
)

print(as.data.frame(cells), digits = 3)
message(
  "\ngate ", if (gate_open) "OPEN" else "SHUT",
  ": ", sum(cells$p_empirical < GATE_ALPHA, na.rm = TRUE), " of ", nrow(cells),
  " cells beat their null at p < ", GATE_ALPHA
)
if (!gate_open) {
  message(
    "no seventh candidate label is written; the sweep stays at six"
  )
}
