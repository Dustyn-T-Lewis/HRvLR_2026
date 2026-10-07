test_that("relabel_samples keeps only grouped people and puts the new labels in arm and group", {
  meta <- toy_meta()
  out <- relabel_samples(meta, c(HR_S1 = "Lo", LR_S3 = "Hi", HR_S2 = "Hi"))
  expect_equal(out$sample_id, meta$sample_id[meta$subject != "LR_S4"])
  expect_equal(unique(out$arm), c("Lo", "Hi"))
  expect_equal(as.character(out$group[1:3]), c("Lo_T1", "Lo_T2", "Lo_T3"))
  expect_true(is.factor(out$group))
})

test_that("relabel_samples with the study's own labels changes nothing", {
  meta <- toy_meta()
  own <- with(distinct(meta, subject, arm), set_names(arm, subject))
  expect_identical(relabel_samples(meta, own), meta)
})

test_that("for_design cuts each kind of processing object to the grouped samples", {
  meta <- toy_meta()
  values <- matrix(seq_len(24), 2, 12, dimnames = list(c("P1", "P2"), meta$sample_id))
  design <- list(
    group_of = c(HR_S1 = "Hi", LR_S3 = "Lo"),
    dal = list(metadata = relabel_samples(meta, c(HR_S1 = "Hi", LR_S3 = "Lo")))
  )
  kept <- design$dal$metadata$sample_id

  dal <- for_design(list(metadata = meta, data = values, observed = values > 5), design)
  expect_equal(colnames(dal$data), kept)
  expect_equal(colnames(dal$observed), kept)
  expect_equal(dal$data, values[, kept])
  expect_equal(unique(dal$metadata$arm), c("Hi", "Lo"))

  modules <- for_design(list(meta = meta, eigengenes = values), design)
  expect_equal(colnames(modules$eigengenes), kept)

  scores <- for_design(list(scores = values), design)
  expect_equal(colnames(scores$scores), kept)
})

test_that("subject_values and subject_change give one column per person", {
  # HR_S2 has no T1 biopsy, so it drops out of the T1 values and the change.
  meta <- toy_meta()[-4, ]
  values <- matrix(seq_len(11), 1, 11, dimnames = list("P1", meta$sample_id))
  t1 <- subject_values(values, meta, "T1")
  expect_equal(colnames(t1), c("HR_S1", "LR_S3", "LR_S4"))
  expect_equal(unname(t1[1, ]), c(1, 6, 9))
  change <- subject_change(values, meta, "T1", "T2")
  expect_equal(colnames(change), c("HR_S1", "LR_S3", "LR_S4"))
  expect_equal(unname(change[1, ]), c(1, 1, 1))
})

test_that("permutation_p is (b + 1) / (m + 1) and never zero", {
  expect_equal(permutation_p(5, c(1, 2, 6, 7)), 3 / 5)
  expect_equal(permutation_p(10, 1:9), 1 / 10)
})
