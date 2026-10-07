# Checks on the saved outputs of a full run. They skip on a fresh clone before run_all.R.
roots <- here::here(c(
  "01_Processing", "02_HRvLR_groups", "03_Cluster_groups", "04_Continuous_no_group"
))
outputs <- list.files(roots, recursive = TRUE, full.names = TRUE)
outputs <- outputs[grepl("/(b_reports|c_data)/", outputs)]
workbooks <- outputs[grepl("[.]xlsx$", outputs)]
result_sheets <- \(f) setdiff(readxl::excel_sheets(f), c("read_me", "columns"))
# Reads only the named columns; the largest sheets hold hundreds of thousands of rows.
read_columns <- function(f, sheet, cols) {
  present <- names(readxl::read_excel(f, sheet, n_max = 0))
  types <- ifelse(present %in% cols, "guess", "skip")
  readxl::read_excel(f, sheet, col_types = types, guess_max = 1e5)
}
skip_without_outputs <- \() skip_if(length(workbooks) == 0, "no outputs; run run_all.R first")

test_that("every workbook opens with read_me and columns", {
  skip_without_outputs()
  for (f in workbooks) expect_equal(readxl::excel_sheets(f)[1:2], c("read_me", "columns"), info = f)
})

test_that("no output file keeps a step number in its name", {
  skip_without_outputs()
  expect_false(any(grepl("^[0-9]{2}_", basename(outputs))))
})

test_that("every estimate carries an interval that contains it", {
  skip_without_outputs()
  estimates <- c("logFC", "slope", "slope_difference", "interaction", "rho", "auc")
  for (f in workbooks) {
    for (sheet in result_sheets(f)) {
      cols <- names(readxl::read_excel(f, sheet, n_max = 0))
      # WGCNA's soft-threshold slope and log2 FC without a test are not estimates with intervals.
      found <- intersect(estimates, cols)
      if (sheet == "soft_threshold") found <- character()
      if (!"t" %in% cols) found <- setdiff(found, "logFC")
      if (!length(found)) next
      d <- read_columns(f, sheet, c(found, paste0(found, "_lower"), paste0(found, "_upper")))
      for (est in found) {
        lower <- d[[paste0(est, "_lower")]]
        upper <- d[[paste0(est, "_upper")]]
        info <- paste(basename(f), sheet, est)
        expect_false(is.null(lower) || is.null(upper), info = info)
        inside <- lower <= d[[est]] + 1e-9 & d[[est]] <= upper + 1e-9
        expect_true(all(inside, na.rm = TRUE), info = info)
      }
    }
  }
})

test_that("every design names its nine contrasts from its two groups", {
  skip_without_outputs()
  designs <- outputs[basename(outputs) == "design.rds"]
  for (f in designs) {
    d <- readRDS(f)
    a <- d$arms[1]
    b <- d$arms[2]
    expect_equal(colnames(d$dal$design$contrast_matrix), c(
      paste0("Training_", c(a, b)), paste0("Acute_", c(a, b)),
      paste0(c("Baseline_", "Trained_", "Acute_"), a, "v", b),
      "Training_Interaction", "Acute_Interaction"
    ), info = f)
    roles <- d$roles$role[d$roles$role != "descriptive"]
    expect_equal(sort(roles), c("floor", "primary", "secondary"), info = f)
    expect_true(all(d$dal$metadata$subject %in% names(d$group_of)), info = f)
  }
})

test_that("classification covers the five between-group tasks", {
  skip_without_outputs()
  for (f in workbooks[grepl("_classification[.]xlsx$", workbooks)]) {
    tasks <- unique(read_columns(f, "auc", "task")$task)
    expect_setequal(
      sub("_[^_]*$", "", tasks),
      c("Baseline", "Trained", "Acute", "Training_change", "Acute_change")
    )
  }
})

test_that("stacked proteome workbooks tag every comparison with its input and validity", {
  skip_without_outputs()
  stacked <- workbooks[grepl("_all_inputs[.]xlsx$", workbooks)]
  for (f in stacked) {
    for (sheet in setdiff(result_sheets(f), c("input_manifest", "package_versions"))) {
      d <- read_columns(f, sheet, c("input", "validity", "contrast", "task", "window"))
      info <- paste(basename(f), sheet)
      expect_false(anyNA(d$input), info = info)
      key <- intersect(c("contrast", "task", "window"), names(d))
      if (length(key)) expect_false(anyNA(d$validity[!is.na(d[[key[1]]])]), info = info)
    }
  }
})

test_that("removed methods leave no sheets behind", {
  skip_without_outputs()
  for (f in workbooks) {
    expect_false(any(c("treat", "ancova") %in% readxl::excel_sheets(f)), info = f)
  }
  for (f in workbooks[grepl("module_preservation[.]xlsx$", workbooks)]) {
    expect_false("calibration" %in% readxl::excel_sheets(f), info = f)
  }
})
