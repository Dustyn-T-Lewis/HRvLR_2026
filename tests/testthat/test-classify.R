# Eight subjects, four per arm, three timepoints. Feature `up` rises from T1 to
# T2 in every subject and sits higher in HR at baseline; `flat` is noise.
fixture <- function() {
  set.seed(42)
  subjects <- sprintf("S%02d", 1:8)
  meta <- tidyr::expand_grid(
    subject = subjects, timepoint = c("T1", "T2", "T3")
  )
  meta$arm <- ifelse(meta$subject %in% subjects[1:4], "HR", "LR")
  meta$sample_id <- paste(meta$subject, meta$timepoint, sep = "_")
  up <- 10 + (meta$timepoint != "T1") * 2 + (meta$arm == "HR") * 3 +
    rnorm(nrow(meta), sd = 0.1)
  mat <- rbind(up = up, flat = rnorm(nrow(meta)))
  colnames(mat) <- meta$sample_id
  list(meta = meta, mat = mat)
}

# Five subjects per arm at T1 only: `hi` separates the arms perfectly, `lo`
# the other way round, and `sparse` is observed in two HR subjects.
arm_fixture <- function() {
  meta <- tibble::tibble(
    subject = sprintf("S%02d", 1:10), timepoint = "T1",
    arm = rep(c("HR", "LR"), each = 5)
  )
  meta$sample_id <- paste(meta$subject, "T1", sep = "_")
  mat <- rbind(
    hi = c(6:10, 1:5), lo = c(1:5, 6:10),
    sparse = c(1, 2, rep(NA, 8))
  )
  colnames(mat) <- meta$sample_id
  list(meta = meta, mat = mat)
}

test_that("AUC reads above 0.5 as higher in the case group", {
  fx <- arm_fixture()
  expect_warning(
    res <- classify_features(
      fx$mat, TASKS[TASKS$task == "Baseline_HRvLR", ], fx$meta
    ),
    "less than 1"
  )
  get <- function(f, col) res[[col]][res$feature == f]
  expect_equal(get("hi", "auc"), 1)
  expect_equal(get("lo", "auc"), 0)
  expect_lt(get("hi", "p"), 0.05)
  expect_equal(
    get("hi", "p"),
    wilcox.test(6:10, 1:5, exact = FALSE)$p.value
  )
})

test_that("a feature with too few observations per group is not scored", {
  fx <- arm_fixture()
  expect_warning(
    res <- classify_features(
      fx$mat, TASKS[TASKS$task == "Baseline_HRvLR", ], fx$meta
    ),
    "less than 1"
  )
  expect_true(is.na(res$auc[res$feature == "sparse"]))
  expect_true(is.na(res$p[res$feature == "sparse"]))
  expect_equal(res$bh[1:2], p.adjust(res$p[1:2], "BH"))
})

test_that("paired tasks use complete pairs and the signed-rank p", {
  meta <- tidyr::expand_grid(
    subject = sprintf("S%02d", 1:5), timepoint = c("T1", "T2")
  )
  meta$arm <- "HR"
  meta$sample_id <- paste(meta$subject, meta$timepoint, sep = "_")
  t1 <- c(1, 2, NA, 4, 5)
  t2 <- c(2, 3, 9, 5, 6)
  mat <- matrix(c(rbind(t1, t2)),
    nrow = 1,
    dimnames = list("f", meta$sample_id)
  )
  res <- classify_features(mat, TASKS[TASKS$task == "Training_HR", ], meta)
  expect_equal(res$auc, 11 / 16)
  keep <- !is.na(t1)
  expect_equal(
    res$p,
    wilcox.test(t2[keep], t1[keep], paired = TRUE, exact = FALSE)$p.value
  )
})

test_that("a paired time task keeps one arm and matches subjects", {
  fx <- fixture()
  task <- TASKS[TASKS$task == "Training_HR", ]
  m <- task_matrices(fx$mat, task, fx$meta)
  expect_equal(colnames(m$control), sprintf("S%02d", 1:4))
  expect_identical(colnames(m$control), colnames(m$case))
})

test_that("an arm task splits subjects by arm on the right window", {
  fx <- fixture()
  m <- task_matrices(fx$mat, TASKS[TASKS$task == "Training_HRvLR", ], fx$meta)
  expect_equal(colnames(m$case), sprintf("S%02d", 1:4))
  expect_equal(colnames(m$control), sprintf("S%02d", 5:8))
  expect_equal(unname(m$case["up", ]), rep(2, 4), tolerance = 0.5)
})

test_that("classify_all recovers the planted effects and BHs within task", {
  fx <- fixture()
  res <- classify_all(fx$mat, fx$meta)
  expect_equal(nrow(res), nrow(TASKS) * 2)
  get <- function(t, f, col) res[[col]][res$task == t & res$feature == f]
  expect_equal(get("Training_HR", "up", "auc"), 1)
  expect_equal(get("Baseline_HRvLR", "up", "auc"), 1)
  one <- res[res$task == "Acute_LR", ]
  expect_equal(one$bh, p.adjust(one$p, "BH"))
})

test_that("associate_all crosses windows with phenotypes", {
  fx <- fixture()
  pheno <- tibble::tibble(
    subject = sprintf("S%02d", 1:8), a = rnorm(8), b = rnorm(8)
  )
  res <- associate_all(fx$mat, fx$meta, pheno,
    windows = c("T1", "training"), phenotypes = c("a", "b")
  )
  expect_equal(nrow(res), 2 * 2 * 2)
  expect_setequal(unique(res$window), c("T1", "training"))
  expect_true(all(res$n == 8))
})

test_that("chance_table counts nominal hits against alpha", {
  res <- tibble::tibble(
    task = rep(c("a", "b"), each = 4),
    p = c(0.01, 0.2, 0.5, NA, 0.9, 0.8, 0.7, 0.6),
    bh = c(0.04, 0.4, 0.5, NA, NA, NA, NA, NA)
  )
  out <- chance_table(res, "task")
  expect_equal(out$n_tested, c(3, 4))
  expect_equal(out$n_nominal, c(1, 0))
  expect_equal(out$n_expected, c(0.15, 0.2))
  expect_equal(out$n_bh, c(1, 0))
  expect_equal(out$min_bh, c(0.04, NA))
})

test_that("screen_level returns both screens and one chance row per cell", {
  fx <- fixture()
  pheno <- tibble::tibble(subject = sprintf("S%02d", 1:8))
  pheno[PHENOTYPES] <- replicate(length(PHENOTYPES), rnorm(8), simplify = FALSE)
  out <- screen_level(fx$mat, fx$meta, pheno)
  expect_named(
    out, c("classify", "associate", "chance_classify", "chance_associate")
  )
  expect_equal(nrow(out$chance_classify), nrow(TASKS))
  expect_equal(nrow(out$chance_associate), 3 * length(PHENOTYPES))
})
