test_that("read_goslim_bp keeps the biological-process terms only", {
  skip_if_not_installed("GO.db")
  source(here::here("functions", "shared_pathway_utils.R"), local = TRUE)
  obo <- withr::local_tempfile(lines = c(
    "format-version: 1.2", "",
    "[Term]", "id: GO:0003012", "name: muscle system process",
    "namespace: biological_process", "",
    "[Term]", "id: GO:0005694", "name: chromosome",
    "namespace: cellular_component"
  ))
  out <- read_goslim_bp(obo)
  expect_equal(out$go_id, "GO:0003012")
  expect_equal(out$name, "muscle system process")
})

test_that("goslim_theme picks the narrowest containing slim term", {
  skip_if_not_installed("GO.db")
  source(here::here("functions", "shared_pathway_utils.R"), local = TRUE)
  slim <- tibble::tibble(
    go_id = c("GO:0008150", "GO:0003012"),
    name = c("biological_process", "muscle system process")
  )
  themes <- goslim_theme(c("GO:0006936", "GO:0003012", "GO:0006096"), slim)
  expect_equal(themes[1:2], rep("muscle system process", 2))
  expect_equal(themes[3], "biological_process")
})

test_that("goslim_sets names one set per slim term", {
  skip_if_not_installed("org.Hs.eg.db")
  source(here::here("functions", "shared_pathway_utils.R"), local = TRUE)
  slim <- tibble::tibble(go_id = "GO:0003012", name = "muscle system process")
  sets <- suppressMessages(goslim_sets(slim))
  expect_named(sets, "GOSLIM_MUSCLE_SYSTEM_PROCESS")
  expect_true("DMD" %in% sets[[1]])
})
