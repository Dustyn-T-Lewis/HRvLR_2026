# Three responder groups per phenotype, and what shape their difference takes.
#
# Each phenotype is cut into thirds, and every feature is fitted once against
# the ordered factor. That one fit yields three read-outs: a linear trend
# across LR < MR < HR, a quadratic deviation for a middle group that sits above
# or below both ends, and an omnibus F for any group difference at all.
#
# The quadratic term is the reason this design exists. A two-group split and a
# linear regression are both blind to a non-monotonic relationship, so it is
# the only question here that the continuous sweep could not already have
# answered. The linear term is reported alongside it as the comparison, and is
# expected to reproduce what the continuous fit found, with less power.
#
# Features are dropped unless every group has at least three observations. The
# protein matrix is 12% missing and a change window needs both timepoints, so
# without that filter a protein measured in one subject of one group still
# returns a group contrast. n_dropped_sparse records the cost per cell.
#
# BH within each cell and term, never pooled. 03_confirm_hits.R supplies the
# sweep-level correction, which is what decides.

pacman::p_load(here, dplyr, purrr, tidyr, tibble, readr, openxlsx)

source(here("functions", "tertile_groups.R"))

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

groups <- map(PHENOTYPES, function(v) {
  tertile_split(pheno[[v]], pheno$subject)
}) |>
  setNames(PHENOTYPES)

# How much do the ten group-sets repeat each other? Two phenotypes cutting the
# cohort identically are one test, not two.
partition_key <- function(g) {
  paste(as.character(g$group[order(g$subject)]), collapse = "")
}
partitions <- tibble(
  phenotype = PHENOTYPES,
  key = map_chr(groups, partition_key)
) |>
  mutate(partition = as.integer(factor(.data$key, levels = unique(.data$key))))

windows <- tibble(level = rep(names(features), each = length(WINDOWS))) |>
  mutate(window = rep(names(WINDOWS), times = length(features))) |>
  mutate(mat = map2(.data$level, .data$window, function(l, w) {
    subject_window(features[[l]], w)
  }))

cells <- expand_grid(
  level = names(features), window = names(WINDOWS), phenotype = PHENOTYPES
)

results <- pmap_dfr(cells, function(level, window, phenotype) {
  idx <- which(windows$level == level & windows$window == window)
  res <- fit_tertile(windows$mat[[idx]], groups[[phenotype]])
  res |>
    mutate(
      level = level, window = window, phenotype = phenotype,
      n = attr(res, "n"), sizes = attr(res, "sizes"),
      n_dropped = attr(res, "n_dropped"), .before = 1
    )
})

summary_tbl <- results |>
  summarise(
    n_subjects = dplyr::first(.data$n),
    sizes = dplyr::first(.data$sizes),
    n_features = dplyr::n(),
    n_dropped_sparse = dplyr::first(.data$n_dropped),
    n_nominal = sum(.data$p < 0.05, na.rm = TRUE),
    n_bh = sum(.data$bh < BH_ALPHA, na.rm = TRUE),
    min_bh = min(.data$bh, na.rm = TRUE),
    .by = c(level, window, phenotype, term)
  ) |>
  left_join(
    partitions |> dplyr::select(phenotype, partition),
    by = "phenotype"
  ) |>
  arrange(.data$min_bh)

survivors <- results |>
  filter(.data$bh < BH_ALPHA) |>
  left_join(annotation, by = "feature") |>
  arrange(.data$bh)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(results, file.path(OUT_DIR, "02_tertile_full.csv"))
write.xlsx(
  list(
    tertile_summary = summary_tbl, survivors = survivors,
    group_assignments = imap_dfr(groups, ~ mutate(.x, phenotype = .y)),
    partitions = partitions
  ),
  file.path(OUT_DIR, "02_tertile.xlsx")
)

print(as.data.frame(head(summary_tbl, 12)), digits = 3)
message(
  "\n", nrow(survivors), " feature-cell-term survivors at BH < ", BH_ALPHA,
  " across ", nrow(summary_tbl), " cells (",
  dplyr::n_distinct(summary_tbl$partition), " distinct partitions of ",
  length(PHENOTYPES), " phenotypes)"
)
message(
  "by term: ",
  paste(
    sprintf(
      "%s %d", names(table(survivors$term)), as.integer(table(survivors$term))
    ),
    collapse = " | "
  )
)
