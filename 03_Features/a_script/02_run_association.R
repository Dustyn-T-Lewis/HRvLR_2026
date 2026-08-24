# Does a proteome change track an adaptation?
#
# Three feature levels by two windows by ten phenotypes. The training window
# asks whether the protein change over the block tracks how much a subject
# adapted; the acute window asks whether the response to a single bout in a
# trained muscle tracks the same thing. Both regress a per-subject change on a
# per-subject phenotype, so nobody is cut into a group and no cut point has to
# be defended.
#
# There is no baseline window. A baseline association compares levels between
# people, which answers a different question from whether a change tracks a
# change, and V1 already tested the baseline form across 54 cells without
# promoting anything.
#
# BH within each cell, never across the sweep: the phenotypes are correlated
# (three fibre-area measures share r > 0.9) and both windows draw on the same
# subjects. The cell count is reported instead, and 03_confirm_hits.R supplies
# what it takes to read that count - how many cells an all-null sweep of this
# shape would produce.

pacman::p_load(here, dplyr, purrr, tidyr, tibble, readr, openxlsx)

source(here("functions", "association.R"))

OUT_DIR <- here("03_Features", "c_data")
BH_ALPHA <- 0.05

pheno <- phenotype_table()
PHENOTYPES <- setdiff(names(pheno), c("subject", "group_arm"))

features <- feature_matrices()
annotation <- read_csv(
  here("02_Normalization", "c_data", "normalized.csv"),
  show_col_types = FALSE
) |>
  dplyr::select(feature = uniprot_id, gene, protein, description)

changes <- expand_grid(level = names(features), window = names(WINDOWS)) |>
  mutate(mat = map2(.data$level, .data$window, function(l, w) {
    subject_change(features[[l]], w)
  }))

cells <- expand_grid(
  level = names(features), window = names(WINDOWS), phenotype = PHENOTYPES
)

results <- pmap_dfr(cells, function(level, window, phenotype) {
  idx <- which(changes$level == level & changes$window == window)
  feat <- changes$mat[[idx]]
  res <- associate(feat, phenotype_vector(pheno, phenotype))
  res |>
    mutate(
      level = level, window = window, phenotype = phenotype,
      n = attr(res, "n"), .before = 1
    )
})

summary_tbl <- results |>
  summarise(
    n_subjects = dplyr::first(.data$n),
    n_features = dplyr::n(),
    n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
    n_bh = sum(.data$bh < BH_ALPHA, na.rm = TRUE),
    min_bh = min(.data$bh, na.rm = TRUE),
    .by = c(level, window, phenotype)
  ) |>
  arrange(.data$min_bh)

survivors <- results |>
  filter(.data$bh < BH_ALPHA) |>
  left_join(annotation, by = "feature") |>
  arrange(.data$bh)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(results, file.path(OUT_DIR, "02_association_full.csv"))
write.xlsx(
  list(association_summary = summary_tbl, survivors = survivors),
  file.path(OUT_DIR, "02_association.xlsx")
)

print(as.data.frame(head(summary_tbl, 12)), digits = 3)
message(
  "\n", nrow(survivors), " feature-cell survivors at BH < ", BH_ALPHA,
  " across ", nrow(summary_tbl), " cells (",
  dplyr::n_distinct(summary_tbl$level), " levels x ",
  dplyr::n_distinct(summary_tbl$window), " windows x ",
  dplyr::n_distinct(summary_tbl$phenotype), " phenotypes)"
)
