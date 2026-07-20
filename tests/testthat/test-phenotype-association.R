source(here::here("04_Figures", "functions", "f04_association.R"))
source(here::here("functions", "shared_singscore.R"))

make_meta <- function(subjects, timepoints = c("T1", "T2", "T3")) {
  meta <- expand.grid(
    Subject_ID = subjects, Timepoint = timepoints,
    stringsAsFactors = FALSE
  )
  meta$Col_ID <- paste(meta$Subject_ID, meta$Timepoint, sep = "_")
  meta
}

test_that("phase_subject_matrix returns baseline levels and phase deltas", {
  subjects <- paste0("S", 1:4)
  meta <- make_meta(subjects)
  mat <- matrix(seq_len(5 * nrow(meta)),
    nrow = 5,
    dimnames = list(paste0("f", 1:5), meta$Col_ID)
  )

  base <- phase_subject_matrix(mat, meta, "baseline")
  expect_setequal(colnames(base), subjects)
  expect_equal(base[, "S1"], mat[, "S1_T1"], ignore_attr = TRUE)

  train <- phase_subject_matrix(mat, meta, "training")
  expect_equal(train[, "S1"], mat[, "S1_T2"] - mat[, "S1_T1"], ignore_attr = TRUE)
})

test_that("training-phase delta drops subjects missing a timepoint", {
  meta <- make_meta(paste0("S", 1:4))
  meta <- meta[!(meta$Subject_ID == "S4" & meta$Timepoint == "T2"), ]
  mat <- matrix(rnorm(5 * nrow(meta)),
    nrow = 5,
    dimnames = list(paste0("f", 1:5), meta$Col_ID)
  )
  train <- phase_subject_matrix(mat, meta, "training")
  expect_false("S4" %in% colnames(train))
})

test_that("score_singscore drops sets below the size floor", {
  skip_if_not_installed("singscore")
  genes <- paste0("g", 1:40)
  mat <- matrix(rnorm(40 * 5), nrow = 40, dimnames = list(genes, paste0("s", 1:5)))
  sets <- list(big = genes[1:10], tiny = genes[1:2])
  scores <- score_singscore(mat, sets, min_size = 5)
  expect_true("big" %in% rownames(scores))
  expect_false("tiny" %in% rownames(scores))
})
