# The sweep run backwards: let the proteome draw the groups, then ask whether
# those groups differ on any phenotype.
#
# This is not circular. The split is fitted on protein data alone and the
# phenotypes are never shown to the clustering, so a difference afterwards
# would be a real association rather than a restatement of the input. It is the
# only direction in this project that could produce a proteome-defined
# responder group.
#
# 01_run_subtypes.R already showed the proteome carries no cluster structure
# that beats its own null, so a two-group split is being forced here. That is
# stated rather than hidden: the split exists, it just is not evidence of
# subtypes. Whether it tracks a phenotype is a separate question and worth
# asking, because a weak-but-real axis can fail a cluster-structure test and
# still correlate with something.
#
# Four views, because "the proteome" is not one thing: baseline and
# training-change, each on module eigengenes and on the top-variance proteins.
# Ten phenotypes each, BH within a view. The null is exact and cheap - under no
# association the cluster labels are exchangeable across subjects - so the hit
# count is calibrated by shuffling them.

pacman::p_load(
  here, dplyr, tidyr, purrr, tibble, readr, mclust, broom, openxlsx
)

source(here("functions", "feature_levels.R"))

OUT_DIR <- here("03_Features", "03_Subtypes", "c_data")

N_PC <- 2
N_PROTEINS <- 500
N_PERM <- 999
BH_ALPHA <- 0.05

pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
)
PHENOTYPES <- setdiff(names(pheno), c("subject", "group_arm"))

meta <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))$
  metadata |>
  as.data.frame() |>
  transmute(sample_id = Col_ID, subject = Subject_ID, timepoint = Timepoint)

eigengene_wide <- function() {
  read_csv(here("03_Features", "02_WGCNA", "c_data", "wgcna_eigengene.csv"),
    show_col_types = FALSE
  ) |>
    pivot_wider(names_from = "group_id", values_from = "ME")
}

# A subject-by-feature matrix for one timepoint, or the T2 minus T1 difference.
subject_matrix <- function(feature_wide, view) {
  d <- feature_wide |> left_join(meta, by = "sample_id")
  cols <- setdiff(names(feature_wide), "sample_id")
  if (view == "baseline") {
    d <- d |> filter(.data$timepoint == "T1")
    mat <- as.matrix(d[, cols])
    rownames(mat) <- d$subject
    return(mat)
  }
  paired <- d |>
    filter(.data$timepoint %in% c("T1", "T2")) |>
    pivot_longer(all_of(cols), names_to = "feature") |>
    pivot_wider(
      id_cols = c("subject", "feature"), names_from = "timepoint",
      values_from = "value"
    ) |>
    filter(!is.na(.data$T1), !is.na(.data$T2)) |>
    mutate(change = .data$T2 - .data$T1)
  wide <- paired |>
    pivot_wider(
      id_cols = "subject", names_from = "feature", values_from = "change"
    )
  mat <- as.matrix(wide[, setdiff(names(wide), "subject")])
  rownames(mat) <- wide$subject
  mat
}

protein_wide <- function() {
  mat <- protein_matrix()
  mat <- mat[stats::complete.cases(mat), , drop = FALSE]
  keep <- order(matrixStats::rowVars(mat), decreasing = TRUE)[
    seq_len(N_PROTEINS)
  ]
  as_tibble(t(mat[keep, , drop = FALSE]), rownames = "sample_id")
}

cluster_subjects <- function(mat) {
  scores <- stats::prcomp(mat, center = TRUE, scale. = TRUE)$x[, seq_len(N_PC),
    drop = FALSE
  ]
  fit <- mclust::Mclust(scores, G = 2, verbose = FALSE)
  tibble(subject = rownames(mat), cluster = paste0("c", fit$classification))
}

test_phenotypes <- function(assign_tbl) {
  d <- assign_tbl |> left_join(pheno, by = "subject")
  map_dfr(PHENOTYPES, function(v) {
    ok <- !is.na(d[[v]])
    if (dplyr::n_distinct(d$cluster[ok]) < 2 || sum(ok) < 6) {
      return(tibble())
    }
    broom::tidy(stats::t.test(d[[v]][ok] ~ d$cluster[ok])) |>
      transmute(
        phenotype = v, n = sum(ok),
        difference = estimate1 - estimate2,
        d = (estimate1 - estimate2) / stats::sd(d[[v]], na.rm = TRUE),
        p = p.value
      )
  }) |>
    mutate(bh = stats::p.adjust(.data$p, method = "BH"))
}

feature_sets <- list(eigengenes = eigengene_wide(), proteins = protein_wide())
views <- tidyr::expand_grid(
  space = names(feature_sets), view = c("baseline", "training_change")
)

set.seed(42)
results <- pmap(views, function(space, view) {
  mat <- subject_matrix(feature_sets[[space]], view)
  assign_tbl <- cluster_subjects(mat)
  observed <- test_phenotypes(assign_tbl) |>
    mutate(space = space, view = view, .before = 1)

  # Under no association the cluster labels are exchangeable across subjects,
  # so shuffling them is the exact null for this hit count.
  null_hits <- vapply(seq_len(N_PERM), function(i) {
    shuffled <- assign_tbl |> mutate(cluster = sample(.data$cluster))
    sum(test_phenotypes(shuffled)$bh < BH_ALPHA)
  }, numeric(1))

  list(
    tests = observed,
    calib = tibble(
      space = space, view = view,
      n_subjects = nrow(mat), n_features = ncol(mat),
      cluster_sizes = paste(sort(table(assign_tbl$cluster)), collapse = "/"),
      n_bh = sum(observed$bh < BH_ALPHA),
      null_any_hit = mean(null_hits > 0),
      p_empirical =
        (sum(null_hits >= sum(observed$bh < BH_ALPHA)) + 1) / (N_PERM + 1)
    ),
    assignment = assign_tbl |> mutate(space = space, view = view, .before = 1)
  )
})

tests <- map_dfr(results, "tests")
calibration <- map_dfr(results, "calib")
assignments <- map_dfr(results, "assignment")

write.xlsx(
  list(
    calibration = calibration,
    phenotype_tests = tests |> arrange(.data$bh),
    cluster_assignments = assignments
  ),
  file.path(OUT_DIR, "02_cluster_phenotype.xlsx")
)

print(as.data.frame(calibration), digits = 3)
message(
  "\n", sum(calibration$n_bh), " phenotype hits at BH < ", BH_ALPHA,
  " across ", nrow(tests), " tests in ", nrow(calibration), " proteome views"
)
