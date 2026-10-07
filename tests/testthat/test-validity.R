test_that("uses() reads the time points from contrast, task and window names", {
  names <- c(
    "Baseline_C1vC2", "Trained_C1vC2", "Acute_C1vC2", "Training_C1", "Training_Interaction",
    "Acute_C2", "Acute_change_C1vC2", "Training_change_C1vC2", "baseline_level", "post_level",
    "post_bout_level", "subject mean", "T2", "something_else"
  )
  expect_equal(uses(names), c(
    "T1", "T2", "T3", "training", "training", "acute", "acute", "training", "T1", "T2", "T3",
    "all", "T2", NA
  ))
})

test_that("validity() reproduces the table in the Proteome README", {
  used <- c(
    T1 = "Baseline_C1vC2", T2 = "Trained_C1vC2", T3 = "Acute_C1vC2", training = "Training_C1",
    acute = "Acute_C1", all = "subject mean"
  )
  expected <- rbind(
    T1 = c("circular", "clean", "clean", "caution", "clean", "caution"),
    T2 = c("clean", "circular", "clean", "caution", "caution", "caution"),
    T3 = c("clean", "clean", "circular", "clean", "caution", "caution"),
    training = c("caution", "caution", "clean", "circular", "caution", "caution"),
    acute = c("clean", "caution", "caution", "caution", "circular", "caution")
  )
  for (input in rownames(expected)) {
    expect_equal(unname(validity(input, used)), expected[input, ], info = input)
  }
  expect_true(is.na(validity("T1", "something_else")))
})

test_that("stack_inputs joins the five workbooks with input and validity on every row", {
  stage <- withr::local_tempdir()
  input_file <- withr::local_tempfile(lines = "x")
  for (input in inputs_used) {
    dir.create(file.path(stage, input, "c_data"), recursive = TRUE)
    write_workbook(
      list(
        results = tibble(contrast = c("Baseline_C1vC2", "Training_C1"), p = c(0.1, 0.2)),
        plain = tibble(p = 0.3),
        large = tibble(p = 0.5)
      ),
      c("One row per contrast.", "No comparison named.", "Kept per input."),
      file.path(stage, input, "c_data", "steps.xlsx"),
      c(table = "00_Input/table.csv"), input_file
    )
  }
  stack_inputs(stage, leave_out = "large")
  stacked_file <- file.path(stage, "c_data", "steps_all_inputs.xlsx")
  expect_false("large" %in% readxl::excel_sheets(stacked_file))
  expect_match(tail(readxl::read_excel(stacked_file, "read_me")$holds, 1), "Not stacked: large")
  expect_equal(names(readxl::read_excel(stacked_file, "plain")), c("input", "p"))
  stacked <- readxl::read_excel(stacked_file, "results")
  expect_equal(names(stacked)[1:2], c("input", "validity"))
  expect_equal(nrow(stacked), 10)
  expect_equal(
    stacked$validity[stacked$input == "T1"],
    unname(validity("T1", c("Baseline_C1vC2", "Training_C1")))
  )
  expect_equal(stacked$validity[stacked$input == "training"], c("caution", "circular"))
})
