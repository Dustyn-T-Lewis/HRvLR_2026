test_that("describe_columns stops on a column the dictionary does not hold", {
  expect_error(describe_columns(list(x = tibble(not_a_known_column = 1))), "not_a_known_column")
  described <- describe_columns(list(x = tibble(feature = "a", p = 0.1, fdr = 0.2)))
  expect_false(anyNA(described$meaning))
})

test_that("write_workbook writes read_me and columns first and the manifest and versions last", {
  input <- withr::local_tempfile(lines = "some input")
  out <- withr::local_tempfile(fileext = ".xlsx")
  write_workbook(
    list(results = tibble(feature = "a", p = 0.1)), "One row per feature.", out,
    c(table = "00_Input/table.csv"), input
  )
  expect_equal(
    readxl::excel_sheets(out),
    c("read_me", "columns", "results", "input_manifest", "package_versions")
  )
  manifest <- readxl::read_excel(out, "input_manifest")
  expect_equal(manifest$path, "00_Input/table.csv")
  expect_equal(manifest$md5, unname(tools::md5sum(input)))
})

test_that("paginate caps each group at 24 hits, 12 to a page", {
  hits <- tibble(task = rep(c("a", "b"), c(30, 5)))
  out <- paginate(hits, "task")
  expect_equal(count(out, task)$n, c(24, 5))
  expect_equal(max(out$position), 12)
  expect_equal(unique(out$pdf_page), 1:3)
})

test_that("write_workbook adds a note as the last read_me row", {
  input <- withr::local_tempfile(lines = "x")
  out <- withr::local_tempfile(fileext = ".xlsx")
  write_workbook(
    list(results = tibble(p = 0.1)), "One row.", out, c(table = "t.csv"), input,
    note = "Panel a reads no data."
  )
  read_me <- readxl::read_excel(out, "read_me")
  expect_equal(tail(read_me$sheet, 1), "note")
  expect_equal(tail(read_me$holds, 1), "Panel a reads no data.")
})

test_that("heat and heatmap_page build and draw an overview page", {
  suppressPackageStartupMessages(library(ComplexHeatmap))
  m <- matrix(c(-1, 0, 1, 2), 2, dimnames = list(c("a", "b"), c("x", "y")))
  effect <- heat(m, "log2 FC", c(-1, 0, 1), labels = c(a = "A", b = "B"), row_names = TRUE)
  status <- heat(
    matrix(c("p < 0.05", "BH < 0.05", NA, "not significant"), 2, dimnames = dimnames(m)),
    "status",
    labels = c(a = "A", b = "B")
  )
  expect_s4_class(effect, "Heatmap")
  page <- heatmap_page(effect + status, "Overview", "Caption.")
  expect_equal(page$title, "Overview")
  pdf(withr::local_tempfile(fileext = ".pdf"))
  expect_no_error(page$draw())
  dev.off()
})
