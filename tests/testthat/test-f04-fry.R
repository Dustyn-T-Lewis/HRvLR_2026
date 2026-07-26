pacman::p_load(here, testthat, limma, dplyr)
source(here("04_Figures", "shared", "WGCNA", "a_script", "enrichment.R"))

# A two-group design with one module planted to move and one planted flat. fry must find
# the first and not the second.
planted <- function(seed = 3, effect = 2.5) {
  set.seed(seed)
  n <- 20
  group <- rep(c("A", "B"), each = n / 2)
  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- c("A", "B")
  contrast_matrix <- limma::makeContrasts(BvA = B - A, levels = design)

  abund <- matrix(rnorm(60 * n), nrow = 60)
  rownames(abund) <- paste0("g", 1:60)
  colnames(abund) <- paste0("s", 1:n)

  # movers: shifted up in group B. flat: untouched.
  movers <- paste0("g", 1:20)
  abund[movers, group == "B"] <- abund[movers, group == "B"] + effect

  colors <- setNames(
    c(rep("movers", 20), rep("flat", 20), rep("grey", 20)),
    rownames(abund)
  )
  list(
    abund = abund, colors = colors, design = design,
    contrast_matrix = contrast_matrix, block = NULL
  )
}

test_that("module_fry finds a planted module and leaves a flat one alone", {
  d <- planted()
  res <- module_fry(
    d$colors, d$abund, d$design, d$contrast_matrix,
    block = NULL, correlation = NULL
  )

  expect_setequal(res$module, c("movers", "flat"))
  expect_identical(unique(res$contrast), "BvA")

  fdr <- setNames(res$fdr, res$module)
  expect_lt(fdr[["movers"]], 0.05)
  expect_gt(fdr[["flat"]], 0.05)

  expect_identical(res$direction[res$module == "movers"], "Up")
})

test_that("module_fry drops grey and modules smaller than three proteins", {
  d <- planted()
  d$colors[paste0("g", 41:59)] <- "tiny" # 19 grey -> tiny... leave g60 grey
  d$colors[paste0("g", 41:58)] <- "grey"
  d$colors[paste0("g", 59:60)] <- "pair" # only 2 proteins

  res <- module_fry(
    d$colors, d$abund, d$design, d$contrast_matrix,
    block = NULL, correlation = NULL
  )

  expect_false("grey" %in% res$module)
  expect_false("pair" %in% res$module)
})

test_that("module_fry reports the direction of a planted DOWN module", {
  d <- planted(effect = -2.5)
  res <- module_fry(
    d$colors, d$abund, d$design, d$contrast_matrix,
    block = NULL, correlation = NULL
  )
  expect_identical(res$direction[res$module == "movers"], "Down")
  expect_lt(res$fdr[res$module == "movers"], 0.05)
})
