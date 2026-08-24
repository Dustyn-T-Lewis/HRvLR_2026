# The same three contrasts under every candidate label: whether the groups
# differed before training, whether they diverged over the training block, and
# whether they diverged over the acute bout.
#
# Ten labels resolve to eight distinct partitions, because the three fibre-area
# measures cut the cohort identically. Each partition is fitted once and every
# label sharing it inherits the result, so the sweep is 24 tests rather than the
# 30 the label count suggests.
#
# BH within each partition and contrast, never across the sweep: the partitions
# are splits of correlated outcomes on the same 16 subjects, so a q spanning
# them would claim an independence the design does not have. The cell count is
# reported instead, and 02_confirm_hits.R supplies what it would take to read
# that count: how often an arbitrary split produces a survivor.
#
# No blood-index covariate. T3 biopsies are roughly twice as bloody as T1 and
# T2, which is why the acute contrast is written as an interaction: a shift
# affecting both arms alike cancels in the difference of their acute changes.

pacman::p_load(here, dplyr, purrr, tibble, readr, openxlsx)

source(here("functions", "label_contrasts.R"))

OUT_DIR <- here("03_Features", "04_Proteins", "c_data")
BH_ALPHA <- 0.05

labels_long <- read_csv(
  here("03_Features", "01_Responsiveness", "c_data", "02_candidate_labels.csv"),
  show_col_types = FALSE
)
labels_meta <- openxlsx::read.xlsx(
  here(
    "03_Features", "01_Responsiveness", "c_data", "02_candidate_labels.xlsx"
  ),
  "labels_meta"
)

mat <- protein_matrix()
annotation <- read_csv(
  here("02_Normalization", "c_data", "normalized.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(feature = uniprot_id, gene, protein, description)

reps <- labels_meta |>
  slice_max(.data$r2_alone, n = 1, by = "partition", with_ties = FALSE)

fits <- map(reps$label, function(l) {
  lab <- labels_long |> filter(.data$label == l)
  vec <- setNames(lab$level, lab$subject)
  res <- fit_label_contrasts(mat, vec)
  res |>
    mutate(
      partition = reps$partition[reps$label == l],
      n_samples = attr(res, "n_samples"),
      within_cor = attr(res, "within_cor"), .before = 1
    )
})

sweep <- bind_rows(fits) |>
  left_join(
    labels_meta |> dplyr::select(label, partition),
    by = "partition", relationship = "many-to-many"
  )

summary_tbl <- sweep |>
  summarise(
    n_tested = sum(!is.na(.data$p)),
    n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
    n_bh = sum(.data$bh < BH_ALPHA, na.rm = TRUE),
    min_bh = min(.data$bh, na.rm = TRUE),
    .by = c(partition, label, contrast)
  ) |>
  left_join(
    labels_meta |> dplyr::select(label, internal, r2_alone, n_hi, n_lo),
    by = "label"
  ) |>
  arrange(.data$partition, .data$contrast)

survivors <- sweep |>
  filter(.data$bh < BH_ALPHA) |>
  left_join(annotation, by = "feature") |>
  arrange(.data$bh)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(sweep, file.path(OUT_DIR, "01_sweep_full.csv"))
write.xlsx(
  list(
    sweep_summary = summary_tbl,
    survivors = survivors,
    equivalence_v1 = verify_v1_equivalence()
  ),
  file.path(OUT_DIR, "01_sweep.xlsx")
)

print(as.data.frame(summary_tbl), digits = 3)
message(
  "\n", nrow(survivors), " protein-contrast survivors at BH < ", BH_ALPHA,
  " across ", dplyr::n_distinct(summary_tbl$partition) *
    dplyr::n_distinct(summary_tbl$contrast),
  " partition-contrast cells (", nrow(summary_tbl), " label-contrast rows)"
)
