# Does a survivor survive its own subjects?
#
# A regression at n = 14 or 15 can clear BH on the strength of two or three
# points. Two checks that a p-value cannot give: refit dropping each subject in
# turn, and compare the parametric fit against a rank correlation. A hit that
# needs a particular subject, or that a rank test does not see, is a hit about
# that subject rather than about the cohort.
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
    delta <- subject_change(features[[level]], window)
    y <- phenotype_vector(pheno, phenotype)
    shared <- intersect(colnames(delta), names(y))
    shared <- shared[!is.na(y[shared])]
    delta <- delta[, shared, drop = FALSE]
    y <- y[shared]

    full <- associate(delta, y) |> filter(.data$feature == !!feature)
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
      spearman_rho = unname(rho$estimate), spearman_p = rho$p.value,
      loso_min_bh = min(loso$bh), loso_max_bh = max(loso$bh),
      loso_folds_below_alpha = sum(loso$bh < BH_ALPHA),
      loso_folds = nrow(loso),
      most_influential = loso$dropped[which.max(loso$bh)],
      robust = sum(loso$bh < BH_ALPHA) == nrow(loso) &&
        rho$p.value < BH_ALPHA
    )
  }
)

loso_detail <- pmap_dfr(
  survivors |> dplyr::select(level, window, phenotype, feature),
  function(level, window, phenotype, feature) {
    delta <- subject_change(features[[level]], window)
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
  " survivors hold in every leave-one-subject-out fold and by rank correlation"
)
