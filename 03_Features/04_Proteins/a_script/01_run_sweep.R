# The same two contrasts under every candidate label. Baseline asks whether the
# groups differ before training started, Training_Interaction whether their
# training response diverged.
#
# Six labels, twelve tests. BH within each label and contrast, never across the
# sweep: five of the six labels are splits of correlated outcomes on the same 16
# subjects, so a q spanning them would claim an independence the design does
# not have. The cell count is reported instead.
#
# No blood-index covariate here, unlike V1's acute contrasts. The T3 biopsies
# are the bloodier ones and neither of these two contrasts touches T3.

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

sweep <- map_dfr(unique(labels_long$label), function(l) {
  lab <- labels_long |> filter(.data$label == l)
  vec <- setNames(lab$level, lab$subject)
  res <- fit_label_contrasts(mat, vec)
  res |>
    mutate(
      label = l, n_samples = attr(res, "n_samples"),
      within_cor = attr(res, "within_cor"), .before = 1
    )
})

summary_tbl <- sweep |>
  summarise(
    n_tested = sum(!is.na(.data$p)),
    n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
    n_bh = sum(.data$bh < BH_ALPHA, na.rm = TRUE),
    min_bh = min(.data$bh, na.rm = TRUE),
    .by = c(label, contrast)
  ) |>
  left_join(
    labels_meta |> dplyr::select(label, internal, r2_alone, n_hi, n_lo),
    by = "label"
  ) |>
  arrange(desc(.data$r2_alone), .data$contrast)

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
  " across ", nrow(summary_tbl), " label-contrast cells"
)
