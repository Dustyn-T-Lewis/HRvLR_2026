pacman::p_load(here, testthat, dplyr)
source(here("functions", "blood_index_model.R"))

# 16 subjects, 8 per arm. S28 and S29 are partial, matching the roster, so the
# fit has to tolerate an unbalanced design rather than assume 3 per subject.
synthetic_blood <- function(arm_t3 = 0, seed = 21) {
  set.seed(seed)
  subjects <- sprintf("S%02d", 1:16)
  subject_effect <- setNames(rnorm(16, sd = 1.2), subjects)

  expand.grid(
    subject = subjects, Timepoint = c("T1", "T2", "T3"),
    stringsAsFactors = FALSE
  ) |>
    mutate(
      Group = factor(
        ifelse(match(.data$subject, subjects) <= 8, "HR", "LR"),
        levels = c("LR", "HR")
      ),
      Timepoint = factor(.data$Timepoint, levels = c("T1", "T2", "T3")),
      blood_index = 25 + subject_effect[.data$subject] +
        rnorm(dplyr::n(), sd = 0.8) +
        arm_t3 * (.data$Group == "HR" & .data$Timepoint == "T3")
    ) |>
    filter(!(.data$subject == "S15" & .data$Timepoint != "T1")) |>
    filter(!(.data$subject == "S16" & .data$Timepoint == "T1"))
}

test_that("blood_index_arm_test returns every fixed-effect term", {
  got <- blood_index_arm_test(synthetic_blood())

  expect_true(all(
    c("estimate", "se", "df", "t", "p") %in% names(got)
  ))
  expect_true("GroupHR:TimepointT3" %in% got$term)
  expect_true(all(got$p >= 0 & got$p <= 1))
})

test_that("LR is the reference arm and T1 the reference timepoint", {
  got <- blood_index_arm_test(synthetic_blood())
  expect_true("GroupHR" %in% got$term)
  expect_false(any(grepl("GroupLR", got$term)))
  expect_false(any(grepl("TimepointT1", got$term)))
})

test_that("a planted arm-by-T3 effect is recovered at the right magnitude", {
  got <- blood_index_arm_test(synthetic_blood(arm_t3 = 2.5))
  term <- got[got$term == "GroupHR:TimepointT3", ]

  expect_lt(term$p, 0.01)
  expect_equal(term$estimate, 2.5, tolerance = 0.6)
})

test_that("no planted interaction leaves the arm-by-T3 term null", {
  got <- blood_index_arm_test(synthetic_blood(arm_t3 = 0))
  expect_gt(got$p[got$term == "GroupHR:TimepointT3"], 0.05)
})

test_that("the fit tolerates the partial-roster subjects", {
  d <- synthetic_blood()
  per_subject <- table(d$subject)
  expect_true(any(per_subject == 1L))
  expect_true(any(per_subject == 2L))
  expect_no_error(blood_index_arm_test(d))
})

test_that("blood_index_data drops the samples the pipeline excluded", {
  skip_if_not(
    file.exists(here("01_Filtering", "c_data", "DAList_filtered.rds"))
  )
  d <- blood_index_data()
  dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))

  expect_identical(nrow(d), ncol(dal$data))
  expect_false(any(c("S29_T1", "S28_T2", "S28_T3") %in% d$Col_ID))
  expect_identical(levels(d$Group), c("LR", "HR"))
})
