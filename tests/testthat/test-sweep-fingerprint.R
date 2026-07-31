pacman::p_load(here, testthat, withr, openxlsx)
source(here("functions", "sweep_grid.R"))

fake_bundle <- function(n_prot = 20, seed = 5) {
  set.seed(seed)
  m <- matrix(rnorm(n_prot * 6), nrow = n_prot)
  rownames(m) <- paste0("p", seq_len(n_prot))
  colnames(m) <- paste0("s", 1:6)
  list(feature_sets = list(proteins = m, singscore = m[1:5, ]))
}

B0 <- c(0L, 0L)

test_that("the fingerprint is stable for identical input", {
  expect_identical(
    sweep_fingerprint(fake_bundle(), B0), sweep_fingerprint(fake_bundle(), B0)
  )
})

test_that("the fingerprint moves when the feature set changes", {
  a <- sweep_fingerprint(fake_bundle(n_prot = 20), B0)
  b <- sweep_fingerprint(fake_bundle(n_prot = 19), B0)
  expect_false(identical(a, b))
})

test_that("the fingerprint moves when values change but shape does not", {
  a <- sweep_fingerprint(fake_bundle(seed = 1), B0)
  b <- sweep_fingerprint(fake_bundle(seed = 2), B0)
  expect_false(identical(a, b))
})

test_that("the fingerprint moves with the permutation grid", {
  b <- fake_bundle()
  expect_false(identical(
    sweep_fingerprint(b, b_grid = c(0L, 0L)),
    sweep_fingerprint(b, b_grid = c(0L, 200L))
  ))
  expect_identical(
    sweep_fingerprint(b, b_grid = c(0L, 200L)),
    sweep_fingerprint(b, b_grid = c(0L, 200L))
  )
})

seed_cell <- function(root, fingerprint, level = "proteins",
                      config = "total", phenotype = "d_mcsa", model = "rf") {
  write_sweep_cell(
    "X", level, config, phenotype, model,
    list(summary = data.frame(q2 = 1)),
    fingerprint = fingerprint, root_dir = root
  )
}

test_that("a B = 0 leaf never satisfies a B = 200 run", {
  root <- withr::local_tempdir()
  b <- fake_bundle()
  fast <- sweep_fingerprint(b, b_grid = c(0L, 0L))
  full <- sweep_fingerprint(b, b_grid = c(0L, 200L))
  seed_cell(root, fast)

  expect_true(leaf_done("X", "proteins", "total", "rf", fast, root_dir = root))
  expect_false(leaf_done("X", "proteins", "total", "rf", full, root_dir = root))
})

test_that("a leaf absent from the store is not done", {
  root <- withr::local_tempdir()
  expect_false(
    leaf_done("X", "proteins", "total", "rf", "abc", root_dir = root)
  )
})

test_that("a leaf is done only when its fingerprint matches", {
  root <- withr::local_tempdir()
  seed_cell(root, "abc")

  expect_true(
    leaf_done("X", "proteins", "total", "rf", "abc", root_dir = root)
  )
  expect_false(
    leaf_done("X", "proteins", "total", "rf", "zzz", root_dir = root)
  )
})

# A continuous leaf writes one cell per outcome, all from the same fit. If any
# of them predates the current input the whole leaf has to be refitted, so
# leaf_done() requires every phenotype to agree rather than the first it finds.
test_that("a leaf is not done when one of its phenotypes is stale", {
  root <- withr::local_tempdir()
  seed_cell(root, "abc", phenotype = "d_mcsa")
  seed_cell(root, "stale", phenotype = "d_fcsa_I")

  expect_false(
    leaf_done("X", "proteins", "total", "rf", "abc", root_dir = root)
  )
})

test_that("write_sweep_workbook keeps the sheets it was given", {
  root <- withr::local_tempdir()
  p <- file.path(root, "results.xlsx")
  write_sweep_workbook(
    p, list(metrics = data.frame(q2 = 1), folds = data.frame(i = 1:3)),
    fingerprint = "abc"
  )
  expect_true(all(c("metrics", "folds") %in% getSheetNames(p)))
  expect_equal(read.xlsx(p, "folds")$i, 1:3)
})

test_that("input_fingerprint is stable, order-sensitive and short", {
  a <- input_fingerprint(1:5, "x")
  expect_identical(a, input_fingerprint(1:5, "x"))
  expect_false(identical(a, input_fingerprint("x", 1:5)))
  expect_false(identical(a, input_fingerprint(1:5, "y")))
  expect_identical(nchar(a), 12L)
})

test_that("sweep_fingerprint is input_fingerprint over the bundle and grid", {
  b <- fake_bundle()
  expect_identical(
    sweep_fingerprint(b, c(0L, 200L)),
    input_fingerprint(b$feature_sets, c(0L, 200L))
  )
})

test_that("write_sweep_cell round-trips a cell through the store", {
  root <- withr::local_tempdir()
  sheets <- list(
    summary = data.frame(
      level = "proteins", config = "total", outcome = "d_mcsa", model = "rf",
      n = 16L, B = 200L, q2 = 0.41, perm_p_q2 = 0.015
    ),
    null = data.frame(outcome = "d_mcsa", model = "rf", q2 = c(0.1, 0.2))
  )
  write_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "rf", sheets,
    fingerprint = "abc123456789", root_dir = root
  )

  got <- read_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "rf", "summary",
    root_dir = root
  )
  expect_equal(nrow(got), 1L)
  expect_equal(got$q2, 0.41)
  expect_equal(
    nrow(read_sweep_cell(
      "X", "proteins", "total", "d_mcsa", "rf", "null",
      root_dir = root
    )),
    2L
  )
})

test_that("read_sweep_cell returns the sheet without the store's own columns", {
  root <- withr::local_tempdir()
  write_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "rf",
    list(null = data.frame(outcome = "d_mcsa", model = "rf", q2 = 0.1)),
    fingerprint = "abc", root_dir = root
  )
  got <- read_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "rf", "null",
    root_dir = root
  )
  expect_named(got, c("outcome", "model", "q2"))
})

test_that("rewriting a cell replaces its rows instead of stacking them", {
  root <- withr::local_tempdir()
  for (q in c(0.1, 0.9)) {
    write_sweep_cell(
      "X", "proteins", "total", "d_mcsa", "rf",
      list(summary = data.frame(q2 = q)),
      fingerprint = "abc", root_dir = root
    )
  }
  got <- read_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "rf", "summary",
    root_dir = root
  )
  expect_equal(nrow(got), 1L)
  expect_equal(got$q2, 0.9)
})

test_that("cells of one root do not read each other", {
  root <- withr::local_tempdir()
  write_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "rf",
    list(summary = data.frame(q2 = 0.4)),
    fingerprint = "abc", root_dir = root
  )
  write_sweep_cell(
    "X", "proteins", "total", "d_mcsa", "svm",
    list(summary = data.frame(q2 = 0.8)),
    fingerprint = "abc", root_dir = root
  )
  expect_equal(
    read_sweep_cell("X", "proteins", "total", "d_mcsa", "rf", "summary",
      root_dir = root
    )$q2,
    0.4
  )
  expect_equal(nrow(read_sweep_store("X", "summary", root_dir = root)), 2L)
})

# The composites rank cells and keep the top 12, so a tie is broken by table
# order. Sorting on write keeps that from depending on the order the sweep
# happened to run in.
test_that("the store is ordered by its key, not by write order", {
  root <- withr::local_tempdir()
  for (cfg in c("total", "T1", "acute")) {
    write_sweep_cell(
      "X", "proteins", cfg, "d_mcsa", "rf",
      list(summary = data.frame(q2 = 0.4)),
      fingerprint = "abc", root_dir = root
    )
  }
  expect_equal(
    read_sweep_store("X", "summary", root_dir = root)$config,
    sort(c("total", "T1", "acute"))
  )
})
