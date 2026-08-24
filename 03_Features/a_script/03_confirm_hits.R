# What an all-null sweep of this shape produces.
#
# The association sweep returned one survivor across 60 cells. That number is
# unreadable on its own: BH controls the false discovery rate inside a cell
# given its own null, and a cell over 12 module eigengenes is a far weaker
# filter than the same alpha over 1900 proteins. So every cell is permuted, not
# only the one that hit, which gives both the per-cell p-value and the number
# of hit cells to expect when nothing is there.
#
# The phenotype is shuffled across subjects. That is the exact null here: under
# no association, which subject carries which adaptation is arbitrary, and the
# proteome change matrix is left exactly as measured.

pacman::p_load(here, dplyr, purrr, tidyr, tibble, readr, openxlsx)

source(here("functions", "association.R"))

OUT_DIR <- here("03_Features", "c_data")
BH_ALPHA <- 0.05
N_PERM <- 999

pheno <- phenotype_table()
PHENOTYPES <- setdiff(names(pheno), c("subject", "group_arm"))

features <- feature_matrices()

changes <- expand_grid(level = names(features), window = names(WINDOWS)) |>
  mutate(mat = map2(.data$level, .data$window, function(l, w) {
    subject_change(features[[l]], w)
  }))

cells <- expand_grid(
  level = names(features), window = names(WINDOWS), phenotype = PHENOTYPES
)

count_hits <- function(feat, y) {
  sum(associate(feat, y)$bh < BH_ALPHA, na.rm = TRUE)
}

set.seed(42)
confirmation <- pmap_dfr(cells, function(level, window, phenotype) {
  idx <- which(changes$level == level & changes$window == window)
  feat <- changes$mat[[idx]]
  y <- phenotype_vector(pheno, phenotype)
  observed <- count_hits(feat, y)

  null <- vapply(seq_len(N_PERM), function(perm) {
    count_hits(feat, setNames(sample(unname(y)), names(y)))
  }, numeric(1))

  tibble(
    level = level, window = window, phenotype = phenotype,
    n_bh = observed,
    null_any_hit = mean(null > 0),
    null_hits_mean = mean(null),
    p_empirical = (sum(null >= observed) + 1) / (N_PERM + 1)
  )
})

# Per-cell null rates differ by an order of magnitude across levels, because
# BH over 12 modules is a much weaker filter than BH over 1900 proteins. The
# expectation has to be summed cell by cell rather than taken from one rate.
sweep_calibration <- confirmation |>
  summarise(
    n_cells = dplyr::n(),
    observed_hit_cells = sum(.data$n_bh >= 1),
    expected_hit_cells = sum(.data$null_any_hit),
    .by = level
  ) |>
  bind_rows(
    tibble(
      level = "all",
      n_cells = nrow(confirmation),
      observed_hit_cells = sum(confirmation$n_bh >= 1),
      expected_hit_cells = sum(confirmation$null_any_hit)
    )
  )

confirmed <- confirmation |> filter(.data$n_bh >= 1, .data$p_empirical < 0.05)

write_csv(confirmation, file.path(OUT_DIR, "03_confirmation.csv"))
write.xlsx(
  list(confirmation = confirmation, sweep_calibration = sweep_calibration),
  file.path(OUT_DIR, "03_confirmation.xlsx")
)

print(as.data.frame(sweep_calibration), digits = 3)
print(as.data.frame(filter(confirmation, .data$n_bh >= 1)), digits = 3)
message(
  "\n", nrow(confirmed), " of ", sum(confirmation$n_bh >= 1),
  " hit cells survive a ", N_PERM, "-permutation phenotype null; ",
  "a sweep this size expects ",
  round(sweep_calibration$expected_hit_cells[
    sweep_calibration$level == "all"
  ], 1), " hit cells with nothing there"
)
