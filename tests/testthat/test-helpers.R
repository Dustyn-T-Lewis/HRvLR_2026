suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(purrr)
  library(tibble)
})
source(test_path("..", "..", "R", "helpers.R"))

toy_meta <- function() {
  expand.grid(timepoint = c("T1", "T2", "T3"), subject = sprintf("S%02d", 1:10)) |>
    as_tibble() |>
    mutate(
      subject = as.character(subject), timepoint = as.character(timepoint),
      arm = if_else(subject %in% sprintf("S%02d", 1:5), "HR", "LR"),
      sample_id = paste(subject, timepoint, sep = "_"),
      group = factor(paste(arm, timepoint, sep = "_"))
    )
}

test_that("remove_blood recovers a noiseless slope and leaves the group means", {
  meta <- toy_meta()
  set.seed(1)
  blood <- rnorm(nrow(meta))
  blood <- blood - mean(blood)
  design <- model.matrix(~ 0 + group + blood, data.frame(group = meta$group, blood = blood))
  means <- as.numeric(meta$group)
  values <- rbind(
    exact = means + 0.8 * blood,
    noisy = means + rnorm(nrow(meta), sd = 0.3),
    sparse = c(1, 2, 3, rep(NA, nrow(meta) - 3))
  )
  colnames(values) <- meta$sample_id
  out <- suppressWarnings(remove_blood(values, design, meta$subject))
  expect_equal(unname(out$slope["exact"]), 0.8, tolerance = 1e-8)
  expect_equal(unname(out$values["exact", ]), means, tolerance = 1e-8)
  expect_equal(out$values["sparse", ], values["sparse", ])
})

test_that("ancova matches lm and needs three subjects per arm", {
  meta <- toy_meta()
  set.seed(2)
  values <- matrix(rnorm(2 * nrow(meta)), 2, dimnames = list(c("a", "b"), meta$sample_id))
  values["b", meta$arm == "HR" & meta$timepoint == "T1"][1:3] <- NA
  out <- ancova(values, meta)
  wide <- \(tp) values["a", meta$sample_id[meta$timepoint == tp]]
  arm <- factor(meta$arm[meta$timepoint == "T1"], c("LR", "HR"))
  reference <- summary(lm(wide("T2") ~ wide("T1") + arm))$coefficients[3, ]
  row <- filter(out, feature == "a", window == "training")
  expect_equal(row$effect, reference[["Estimate"]])
  expect_equal(row$p, reference[["Pr(>|t|)"]])
  expect_equal(c(row$n_HR, row$n_LR), c(5, 5))
  expect_false(any(out$feature == "b" & out$window == "training"))
  expect_true(any(out$feature == "b" & out$window == "acute"))
})

test_that("classify_row gives the Mann-Whitney AUC, DeLong limits and paired tests", {
  set.seed(3)
  hr <- matrix(rnorm(7, 1), 1, dimnames = list("a", NULL))
  lr <- matrix(rnorm(8), 1, dimnames = list("a", NULL))
  unpaired <- classify_row(list(positive = hr, negative = lr, paired = FALSE), "a")
  w <- wilcox.test(hr[1, ], lr[1, ])$statistic
  expect_equal(unname(unpaired["auc"]), unname(w) / 56)
  delong <- pROC::ci.auc(pROC::roc(
    controls = lr[1, ], cases = hr[1, ], direction = "<", quiet = TRUE
  ))
  expect_equal(unname(unpaired[c("auc_lower", "auc_upper")]), as.numeric(delong)[c(1, 3)])
  expect_equal(unname(unpaired["p"]), wilcox.test(hr[1, ], lr[1, ])$p.value)

  separated <- classify_row(list(positive = hr + 10, negative = lr, paired = FALSE), "a")
  expect_equal(unname(separated["auc"]), 1)
  expect_true(all(is.na(separated[c("auc_lower", "auc_upper")])))

  after <- hr[, 1:7, drop = FALSE] + 0.5
  paired <- classify_row(list(positive = after, negative = hr, paired = TRUE), "a")
  expect_equal(unname(paired["median_change"]), 0.5)
  expect_equal(unname(paired["p"]), wilcox.test(after[1, ], hr[1, ], paired = TRUE)$p.value)
  expect_true(is.na(paired["auc_lower"]))

  two_hr <- list(positive = hr[, 1:2, drop = FALSE], negative = lr, paired = FALSE)
  short <- classify_row(two_hr, "a")
  expect_true(is.na(short["p"]))
})

test_that("spearman_p reproduces cor.test with exact = FALSE", {
  set.seed(4)
  x <- rnorm(14)
  y <- x + rnorm(14)
  rho <- cor(x, y, method = "spearman")
  expect_equal(
    spearman_p(rho, 14),
    cor.test(x, y, method = "spearman", exact = FALSE)$p.value
  )
})

test_that("spearman_by_row skips rows seen in under two thirds of subjects", {
  set.seed(5)
  x <- matrix(rnorm(30), 2, dimnames = list(c("full", "sparse"), NULL))
  x["sparse", 1:6] <- NA
  out <- spearman_by_row(x, rnorm(15))
  expect_false(is.na(out$rho[out$feature == "full"]))
  expect_true(is.na(out$rho[out$feature == "sparse"]))
})

test_that("exact_spearman_p equals cor.test's exact p without ties", {
  set.seed(6)
  x <- rnorm(7)
  y <- x + rnorm(7)
  expect_equal(nrow(unique(orderings(4))), 24)
  expect_equal(
    exact_spearman_p(x, y),
    cor.test(x, y, method = "spearman", exact = TRUE)$p.value,
    tolerance = 1e-6
  )
})

test_that("subject_fit returns lm's coefficients and applies the floor", {
  set.seed(7)
  predictors <- data.frame(growth = rnorm(12), strength = rnorm(12))
  x <- rbind(a = 2 * predictors$growth + rnorm(12, sd = 0.1), b = c(rnorm(3), rep(NA, 9)))
  out <- subject_fit(x, predictors, ~ growth + strength, c("growth", "strength"))
  reference <- coef(lm(x["a", ] ~ growth + strength, predictors))
  expect_equal(
    out$effect[out$feature == "a"], unname(reference[c("growth", "strength")])
  )
  expect_false("b" %in% out$feature)
})

test_that("page and grid helpers place every row", {
  hits <- tibble(task = rep(c("a", "b"), c(13, 2)), feature = as.character(1:15))
  paged <- paginate(hits, "task")
  expect_equal(max(paged$pdf_page), 3)
  expect_equal(paged$position[13], 1)
  expect_length(slots(c("x", "y")), 12)
  grid <- as_grid(
    tibble(feature = c("f1", "f2"), task = c("a", "b"), value = 1:2), c("f1", "f2"),
    c("a", "b"), "task", \(d) d$value
  )
  expect_equal(grid["f2", "b"], 2)
  expect_true(is.na(grid["f1", "b"]))
  expect_equal(
    status_of(c(0.01, 0.01, 0.5, NA), c(0.01, 0.2, 0.9, NA)),
    c("BH < 0.05", "p < 0.05", "not significant", NA)
  )
})

test_that("display and workbook helpers produce what the scripts expect", {
  expect_equal(cap(c(-100, 1:99)), unname(quantile(c(100, 1:99), 0.99)))
  expect_equal(
    counts_line(
      list(n_nominal = 3, n_tested = 10, n_fdr_05 = 1), "proteins",
      tibble(page = 1, n_pages = 2)
    ),
    "3 of 10 proteins at p < 0.05, 1 at BH < 0.05. Page 1 of 2."
  )

  skip_if_not_installed("ComplexHeatmap")
  skip_if_not_installed("writexl")
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(grid)
  })
  m <- matrix(c(-1, 0, 1, 2), 2, dimnames = list(c("f1", "f2"), c("a", "b")))
  map_effect <- heat(m, "effect", c(-1, 0, 1), labels = c(f1 = "one", f2 = "two"))
  expect_s4_class(map_effect, "Heatmap")
  page <- heatmap_page(map_effect, "title", "caption")
  expect_equal(page$title, "title")
  pdf(NULL)
  page$draw()
  dev.off()

  inputs <- c(helpers = "R/helpers.R")
  paths <- test_path("..", "..", "R", "helpers.R")
  file <- tempfile(fileext = ".xlsx")
  write_workbook(list(results = tibble(x = 1)), "One row.", file, inputs, paths)
  expect_equal(
    readxl::excel_sheets(file), c("read_me", "results", "input_manifest", "package_versions")
  )
})

test_that("fit_contrasts matches a direct limma fit and contrast_rows lays it out", {
  meta <- toy_meta()
  set.seed(8)
  values <- matrix(rnorm(3 * nrow(meta)), 3, dimnames = list(c("a", "b", "c"), meta$sample_id))
  design <- model.matrix(~ 0 + group, meta)
  colnames(design) <- levels(meta$group)
  contrasts <- makeContrasts(HR_T2 - LR_T2, levels = design)
  fit <- fit_contrasts(values, design, contrasts, meta$subject)
  correlation <- duplicateCorrelation(values, design, block = meta$subject)$consensus.correlation
  direct <- lmFit(values, design, block = meta$subject, correlation = correlation) |>
    contrasts.fit(contrasts) |>
    eBayes()
  expect_equal(fit$p.value, direct$p.value)
  rows <- contrast_rows(fit, "x ")
  expect_equal(rows$test, rep("x HR_T2 - LR_T2", 3))
  expect_equal(rows$p, unname(direct$p.value[, 1]))
})

test_that("leave_one_out counts reruns that stay nominal and names the worst drop", {
  test <- function(keep) {
    tibble(
      feature = "f", test = "t",
      p = if ("s2" %in% keep) 0.01 else 0.2
    )
  }
  out <- leave_one_out(c("s1", "s2", "s3"), test)
  expect_equal(out$n_reruns, 3)
  expect_equal(out$n_nominal, 2)
  expect_equal(out$max_p, 0.2)
  expect_equal(out$worst, "s2")
})

test_that("fry_rows returns fry's p for each set and contrast", {
  meta <- toy_meta()
  set.seed(9)
  values <- matrix(rnorm(6 * nrow(meta)), 6, dimnames = list(letters[1:6], meta$sample_id))
  design <- model.matrix(~ 0 + group, meta)
  colnames(design) <- levels(meta$group)
  contrasts <- makeContrasts(up = HR_T2 - LR_T2, levels = design)
  index <- list(s1 = 1:3, s2 = 4:6)
  out <- fry_rows(values, index, design, contrasts, meta$subject)
  correlation <- duplicateCorrelation(values, design, block = meta$subject)$consensus.correlation
  direct <- fry(values, index, design, contrasts[, "up"],
    block = meta$subject, correlation = correlation, sort = "none"
  )
  expect_equal(out$test, c("fry up", "fry up"))
  expect_equal(out$p, direct$PValue)
})
