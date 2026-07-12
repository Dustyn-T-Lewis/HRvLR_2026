source(here::here("05_test", "phenotype_mapping", "a_script", "functions", "associate.R"))

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

test_that("associate_traits recovers a planted association", {
  set.seed(1)
  subjects <- paste0("S", 1:12)
  meta <- make_meta(subjects)
  trait <- rnorm(length(subjects))
  pheno <- data.frame(subject = subjects, resp = trait)

  mat <- matrix(rnorm(20 * nrow(meta)),
    nrow = 20,
    dimnames = list(paste0("f", 1:20), meta$Col_ID)
  )
  t1 <- meta$Timepoint == "T1"
  mat["f1", t1] <- 3 * trait + rnorm(length(subjects), sd = 0.1)

  res <- associate_traits(mat, meta, pheno, "resp", "baseline")
  expect_equal(res$feature[which.min(res$p)], "f1")
  expect_true(all(res$bh >= res$p - 1e-9))
  expect_true(res$bh[res$feature == "f1"] < 0.05)
})

test_that("singscore pathway score is independent of the sample set (no leakage)", {
  skip_if_not_installed("singscore")
  set.seed(2)
  genes <- paste0("g", 1:60)
  mat <- matrix(rnorm(60 * 6), nrow = 60, dimnames = list(genes, paste0("s", 1:6)))
  sets <- list(setA = genes[1:15])

  full <- score_pathways(mat, sets)
  solo <- score_pathways(mat[, 1, drop = FALSE], sets)
  expect_equal(unname(full["setA", "s1"]), unname(solo["setA", "s1"]), tolerance = 1e-9)
})

test_that("score_pathways drops sets below the size floor", {
  skip_if_not_installed("singscore")
  genes <- paste0("g", 1:40)
  mat <- matrix(rnorm(40 * 5), nrow = 40, dimnames = list(genes, paste0("s", 1:5)))
  sets <- list(big = genes[1:10], tiny = genes[1:2])
  scores <- score_pathways(mat, sets, min_size = 5)
  expect_true("big" %in% rownames(scores))
  expect_false("tiny" %in% rownames(scores))
})
