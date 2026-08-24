# Does a survivor survive its own subjects, and its own biopsies?
#
# A regression at n = 14 or 15 can clear BH on the strength of two or three
# points. Three checks a p-value cannot give: refit dropping each subject in
# turn, compare the parametric fit against a rank correlation, and adjust for
# what the needle collected. A hit that needs a particular subject, that a rank
# test cannot see, or that dissolves once tissue composition is in the model is
# a hit about something other than adaptation.
#
# The composition check matters most for the level windows. At T2 the myofibre
# fraction correlates -0.81 with change in whole-muscle CSA and the blood
# fraction +0.71, so an association between a T2 level and that phenotype is
# partly a statement about biopsy content. Adjustment is reported at three
# depths because spending five covariates on fifteen subjects leaves nine
# residual degrees of freedom, and a hit lost at that price has not been shown
# to be composition.
#
# V1 found the same pattern from the other side: its Delta mCSA axis result
# rested on LR_S14, and dropping that one subject moved the correlations it was
# built on by more than half.
#
# This runs only against cells that already cleared BH, so it never becomes a
# fishing pass over the whole sweep.

pacman::p_load(here, dplyr, purrr, tibble, readr, openxlsx)

source(here("functions", "association.R"))

OUT_DIR <- here("03_Features", "c_data")
BH_ALPHA <- 0.05

survivors <- openxlsx::read.xlsx(
  file.path(OUT_DIR, "02_association.xlsx"), "survivors"
)

if (!nrow(survivors)) {
  message("no survivor to test; nothing written")
  quit(save = "no")
}

pheno <- phenotype_table()
features <- feature_matrices()

robustness <- pmap_dfr(
  survivors |> dplyr::select(level, window, phenotype, feature),
  function(level, window, phenotype, feature) {
    delta <- subject_window(features[[level]], window)
    y <- phenotype_vector(pheno, phenotype)
    shared <- intersect(colnames(delta), names(y))
    shared <- shared[!is.na(y[shared])]
    delta <- delta[, shared, drop = FALSE]
    y <- y[shared]

    full <- associate(delta, y) |> filter(.data$feature == !!feature)

    comp <- composition_matrix(window)[colnames(delta), , drop = FALSE]
    confounded <- names(which(vapply(
      colnames(comp),
      function(p) stats::cor(comp[, p], y, method = "spearman"),
      numeric(1)
    ) |> abs() > 0.5))
    adj_bh <- function(cols) {
      if (!length(cols)) {
        return(NA_real_)
      }
      r <- associate(delta, y, adjust = comp[, cols, drop = FALSE])
      r$bh[r$feature == feature]
    }
    rho <- stats::cor.test(
      delta[feature, ], y,
      method = "spearman", exact = FALSE
    )

    loso <- map_dfr(shared, function(drop) {
      keep <- setdiff(shared, drop)
      res <- associate(delta[, keep, drop = FALSE], y[keep]) |>
        filter(.data$feature == !!feature)
      tibble(dropped = drop, bh = res$bh, p = res$p, slope = res$slope)
    })

    tibble(
      level = level, window = window, phenotype = phenotype, feature = feature,
      n = length(y),
      bh_full = full$bh, p_full = full$p,
      confounders = paste(confounded, collapse = ", "),
      bh_adj_confounders = adj_bh(confounded),
      bh_adj_all_panels = adj_bh(colnames(comp)),
      spearman_rho = unname(rho$estimate), spearman_p = rho$p.value,
      loso_min_bh = min(loso$bh), loso_max_bh = max(loso$bh),
      loso_folds_below_alpha = sum(loso$bh < BH_ALPHA),
      loso_folds = nrow(loso),
      most_influential = loso$dropped[which.max(loso$bh)],
      robust = sum(loso$bh < BH_ALPHA) >= nrow(loso) - 1 &&
        rho$p.value < BH_ALPHA &&
        (is.na(adj_bh(confounded)) || adj_bh(confounded) < BH_ALPHA)
    )
  }
)

loso_detail <- pmap_dfr(
  survivors |> dplyr::select(level, window, phenotype, feature),
  function(level, window, phenotype, feature) {
    delta <- subject_window(features[[level]], window)
    y <- phenotype_vector(pheno, phenotype)
    shared <- intersect(colnames(delta), names(y))
    shared <- shared[!is.na(y[shared])]
    map_dfr(shared, function(drop) {
      keep <- setdiff(shared, drop)
      res <- associate(delta[, keep, drop = FALSE], y[keep]) |>
        filter(.data$feature == !!feature)
      tibble(
        level = level, window = window, phenotype = phenotype,
        feature = feature, dropped = drop, bh = res$bh, slope = res$slope
      )
    })
  }
)

write.xlsx(
  list(robustness = robustness, loso_detail = loso_detail),
  file.path(OUT_DIR, "04_hit_robustness.xlsx")
)

print(as.data.frame(robustness), digits = 3)
message(
  "\n", sum(robustness$robust), " of ", nrow(robustness),
  " survivors hold by rank correlation, in all but at most one ",
  "leave-one-subject-out fold, and after adjusting for the composition ",
  "panels that confound their phenotype"
)
