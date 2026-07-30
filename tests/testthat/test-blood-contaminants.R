pacman::p_load(here, testthat, readr, dplyr, readxl)

# The curated blood list removes 95 proteins by accession and the muscle rescue
# cannot overturn it, so a typo in an accession silently keeps a contaminant and
# a stale one silently keeps nothing. Both are worth failing a build over.
blood_list <- read_csv(
  here("00_input", "blood_contaminants.csv"),
  show_col_types = FALSE
)

test_that("the curated blood list has the columns filter_config binds on", {
  expect_true(all(
    c("uniprot_id", "gene", "class", "reason") %in% names(blood_list)
  ))
  expect_false(any(is.na(blood_list$uniprot_id)))
  expect_false(any(is.na(blood_list$reason)))
})

test_that("every accession is unique and matches the measured proteome", {
  expect_false(any(duplicated(blood_list$uniprot_id)))

  raw <- read_excel(here("00_input", "HRvLR_raw.xlsx"))
  measured <- sub("-[0-9]+$", "", raw$uniprot_id)
  expect_true(all(sub("-[0-9]+$", "", blood_list$uniprot_id) %in% measured))
})

test_that("every member carries one of the five documented classes", {
  expect_setequal(
    unique(blood_list$class),
    c("plasma", "immunoglobulin", "erythrocyte", "complement", "leukocyte")
  )
})

test_that("the list does not collide with the cRAP curated entries", {
  source(here("01_Filtering", "a_script", "filter_config.R"))
  # filter_config binds the two, so a shared accession would duplicate rows and
  # inflate the removal counts in filter_log.
  expect_false(any(duplicated(contaminants$uniprot_id)))
})

test_that("no member is a muscle marker", {
  muscle <- c(
    "MB", "CKM", "ACTA1", "MYH1", "MYH2", "MYH7", "ALDOA", "CASQ1", "PYGM",
    "ATP2A1", "DES", "TNNT3", "TNNI2", "MYBPC1", "TNNC2", "MYL1", "PVALB",
    "ANKRD1", "ANKRD2", "CD36", "GPI", "ANXA2", "PPIA"
  )
  expect_length(intersect(blood_list$gene, muscle), 0L)
})

test_that("the sixteen returned proteins stayed off the list", {
  returned <- c(
    "RBMX", "SYNE2", "HNRNPA3", "DDI2", "GDI2", "STK38", "CPNE3", "CCS",
    "FBXO7", "LXN", "ARF1", "ANP32E", "SEPTIN6", "RAN", "DNAJB4", "FLOT1"
  )
  expect_length(intersect(blood_list$gene, returned), 0L)
})

test_that("the filtered matrix contains no curated blood protein", {
  skip_if_not(
    file.exists(here("01_Filtering", "c_data", "DAList_filtered.rds"))
  )
  dal <- readRDS(here("01_Filtering", "c_data", "DAList_filtered.rds"))
  expect_length(intersect(dal$annotation$gene, blood_list$gene), 0L)
})
