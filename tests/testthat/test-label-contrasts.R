test_that("label_cells crosses the label with timepoint and drops the unnamed", {
  source(here::here("functions", "label_contrasts.R"))

  meta <- tibble::tribble(
    ~sample_id, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2", "S01_T3", "S01", "T3",
    "S02_T1", "S02", "T1", "S02_T2", "S02", "T2", "S02_T3", "S02", "T3",
    "S03_T1", "S03", "T1"
  )
  label <- c(S01 = "hi", S02 = "lo")

  out <- label_cells(label, meta)

  expect_setequal(out$sample_id, setdiff(meta$sample_id, "S03_T1"))
  expect_equal(
    as.character(out$cell[out$sample_id == "S01_T2"]), "hi_T2"
  )
  expect_equal(
    as.character(out$cell[out$sample_id == "S02_T3"]), "lo_T3"
  )
  expect_equal(levels(out$cell), LABEL_CELLS)
})

test_that("label_design errors rather than fitting an unfillable cell", {
  source(here::here("functions", "label_contrasts.R"))

  meta <- tibble::tribble(
    ~sample_id, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2",
    "S02_T1", "S02", "T1", "S02_T2", "S02", "T2"
  )
  mat <- matrix(rnorm(8), nrow = 2, dimnames = list(
    c("p1", "p2"),
    meta$sample_id
  ))
  label <- c(S01 = "hi", S02 = "lo")

  expect_error(label_design(mat, label, meta = meta), "no samples in cell")
})

test_that("the two contrasts encode the differences they claim", {
  source(here::here("functions", "label_contrasts.R"))

  expect_equal(LABEL_CONTRAST_NAMES, c("Baseline", "Training_Interaction"))

  design <- diag(6)
  colnames(design) <- LABEL_CELLS
  cm <- limma::makeContrasts(contrasts = LABEL_CONTRASTS, levels = design)

  baseline <- cm[, 1]
  expect_equal(unname(baseline[c("hi_T1", "lo_T1")]), c(1, -1))
  expect_equal(sum(baseline), 0)

  interaction <- cm[, 2]
  expect_equal(
    unname(interaction[c("hi_T1", "hi_T2", "lo_T1", "lo_T2")]),
    c(-1, 1, 1, -1)
  )
  expect_equal(sum(interaction), 0)
})

test_that("a label that names no subject in a cell cannot silently pass", {
  source(here::here("functions", "label_contrasts.R"))

  meta <- tibble::tribble(
    ~sample_id, ~subject, ~timepoint,
    "S01_T1", "S01", "T1", "S01_T2", "S01", "T2", "S01_T3", "S01", "T3"
  )
  all_hi <- c(S01 = "hi")

  out <- label_cells(all_hi, meta)
  expect_setequal(as.character(out$cell), c("hi_T1", "hi_T2", "hi_T3"))
  expect_false(any(grepl("^lo_", out$cell)))
})
