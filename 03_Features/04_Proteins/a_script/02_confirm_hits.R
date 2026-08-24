# What a random split of these subjects produces, for every partition.
#
# With 24 cells the question is no longer only "is this hit real" but "how many
# hits should a sweep this size produce when nothing is there". BH controls the
# false discovery rate inside a cell given its own null; it says nothing about
# how often an arbitrary split of 16 subjects yields a survivor. So every
# partition is permuted, not only the two that hit, which gives both the
# per-cell p-value and the expected number of hit cells under the null.
#
# Subject is the unit of randomisation because the label is a subject property.
# Timepoint travels with the sample untouched, so the repeated-measures
# structure the contrasts rely on is left intact. The consensus correlation is
# estimated once on the observed data and reused, following V1's
# pi_permutation.R: re-estimating it per permutation costs far more than it
# changes, and the permutation does not disturb the structure it describes.

pacman::p_load(here, dplyr, purrr, tibble, readr, limma, openxlsx)

source(here("functions", "label_contrasts.R"))

OUT_DIR <- here("03_Features", "04_Proteins", "c_data")
BH_ALPHA <- 0.05
N_PERM <- 999

summary_tbl <- openxlsx::read.xlsx(
  file.path(OUT_DIR, "01_sweep.xlsx"), "sweep_summary"
)
labels_long <- read_csv(
  here("03_Features", "01_Responsiveness", "c_data", "02_candidate_labels.csv"),
  show_col_types = FALSE
)

reps <- summary_tbl |>
  distinct(.data$partition, .data$label, .keep_all = TRUE) |>
  slice_max(.data$r2_alone, n = 1, by = "partition", with_ties = FALSE)

mat <- protein_matrix()

# One permuted fit. The design is rebuilt because the cell assignment moved;
# the correlation and block are not, because the subject structure did not.
# `samples` is fixed across permutations: shuffling which subject is hi does
# not change which subjects the label covers, and a label built on d_1rm_ext
# covers one fewer.
fit_once <- function(label_vec, correlation, samples) {
  meta <- label_cells(label_vec)
  meta <- meta[match(samples, meta$sample_id), ]
  design <- stats::model.matrix(~ 0 + cell, meta)
  colnames(design) <- levels(meta$cell)
  cm <- limma::makeContrasts(contrasts = LABEL_CONTRASTS, levels = design)
  colnames(cm) <- LABEL_CONTRAST_NAMES
  fit <- limma::lmFit(mat[, samples, drop = FALSE], design,
    block = meta$subject, correlation = correlation
  )
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))
  map_dfr(LABEL_CONTRAST_NAMES, function(ct) {
    tt <- limma::topTable(fit2,
      coef = ct, number = Inf, adjust.method = "BH", sort.by = "none"
    )
    tibble(
      contrast = ct,
      n_bh = sum(tt$adj.P.Val < BH_ALPHA, na.rm = TRUE),
      min_p = min(tt$P.Value, na.rm = TRUE)
    )
  })
}

set.seed(42)
fitted <- map(reps$label, function(l) {
  lab <- labels_long |> filter(.data$label == l)
  vec <- setNames(lab$level, lab$subject)
  parts <- label_design(mat, vec)
  samples <- colnames(parts$mat)
  observed <- fit_once(vec, parts$correlation, samples)

  null <- map_dfr(seq_len(N_PERM), function(i) {
    shuffled <- setNames(sample(unname(vec)), names(vec))
    fit_once(shuffled, parts$correlation, samples) |> mutate(perm = i)
  })

  list(
    label = l,
    observed = observed |> mutate(label = l, .before = 1),
    null = null |> mutate(label = l, .before = 1)
  )
})

null_draws <- map_dfr(fitted, "null")

confirmation <- map_dfr(fitted, function(x) {
  x$observed |>
    left_join(
      x$null |>
        summarise(
          null_any_hit = mean(.data$n_bh >= 1),
          null_hits_mean = mean(.data$n_bh),
          null_min_p_q05 = stats::quantile(.data$min_p, 0.05),
          .by = "contrast"
        ),
      by = "contrast"
    ) |>
    left_join(
      x$null |>
        dplyr::select(contrast, perm, null_n = n_bh) |>
        left_join(x$observed, by = "contrast") |>
        summarise(
          p_count = (sum(.data$null_n >= .data$n_bh) + 1) / (N_PERM + 1),
          .by = "contrast"
        ),
      by = "contrast"
    )
})

# The sweep-level read: how many cells cleared BH, against how many an
# all-null sweep of the same shape would be expected to produce. The per-cell
# null rate is what supplies the expectation.
sweep_calibration <- tibble(
  n_cells = nrow(confirmation),
  observed_hit_cells = sum(confirmation$n_bh >= 1),
  expected_hit_cells = sum(confirmation$null_any_hit),
  p_binomial = stats::pbinom(
    sum(confirmation$n_bh >= 1) - 1,
    size = nrow(confirmation),
    prob = mean(confirmation$null_any_hit),
    lower.tail = FALSE
  )
)

confirmed <- confirmation |> filter(.data$n_bh >= 1, .data$p_count < 0.05)

write_csv(null_draws, file.path(OUT_DIR, "02_perm_null_draws.csv"))
write.xlsx(
  list(confirmation = confirmation, sweep_calibration = sweep_calibration),
  file.path(OUT_DIR, "02_confirmation.xlsx")
)

print(as.data.frame(confirmation), digits = 3)
print(as.data.frame(sweep_calibration), digits = 3)
message(
  "\n", nrow(confirmed), " of ", sum(confirmation$n_bh >= 1),
  " hit cells survive a ", N_PERM, "-permutation subject-label null; ",
  "a sweep this size expects ", round(sweep_calibration$expected_hit_cells, 1),
  " hit cells with nothing there"
)
