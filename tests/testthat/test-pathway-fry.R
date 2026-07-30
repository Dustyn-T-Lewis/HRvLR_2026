pacman::p_load(here, testthat, limma, dplyr)
source(here("functions", "shared_pathway_utils.R"))

# One set planted to move, one planted flat, and a third whose members are
# mostly absent from the matrix. fry must find the first, clear the second, and
# never see the third.
planted <- function(seed = 7, effect = 2.5) {
  set.seed(seed)
  n <- 20
  group <- rep(c("A", "B"), each = n / 2)
  design <- model.matrix(~ 0 + factor(group))
  colnames(design) <- c("A", "B")
  contrasts <- limma::makeContrasts(BvA = B - A, levels = design)

  expr <- matrix(rnorm(60 * n), nrow = 60)
  rownames(expr) <- paste0("g", 1:60)
  colnames(expr) <- paste0("s", 1:n)
  movers <- paste0("g", 1:20)
  expr[movers, group == "B"] <- expr[movers, group == "B"] + effect

  gene_sets <- list(
    movers = movers,
    flat = paste0("g", 21:40),
    unmeasured = c(paste0("g", 41:45), paste0("absent", 1:20))
  )
  list(
    expr = expr, gene_sets = gene_sets, design = design,
    contrasts = contrasts
  )
}

test_that("pathway_fry finds a planted set and leaves a flat one alone", {
  d <- planted()
  res <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts)

  expect_setequal(res$pathway, c("movers", "flat"))
  expect_identical(unique(res$contrast), "BvA")

  fdr <- setNames(res$fdr, res$pathway)
  expect_lt(fdr[["movers"]], 0.05)
  expect_gt(fdr[["flat"]], 0.05)
  expect_identical(res$direction[res$pathway == "movers"], "Up")
})

test_that("pathway_fry reports the direction of a planted DOWN set", {
  d <- planted(effect = -2.5)
  res <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts)
  expect_identical(res$direction[res$pathway == "movers"], "Down")
  expect_lt(res$fdr[res$pathway == "movers"], 0.05)
})

test_that("pathway_fry counts detected members, not annotated ones", {
  d <- planted()
  res <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts)

  # 'unmeasured' has 25 annotated members but only 5 in the matrix.
  expect_false("unmeasured" %in% res$pathway)
  expect_identical(res$n[res$pathway == "movers"], 20L)
})

test_that("pathway_fry honours a custom detected-size floor", {
  d <- planted()
  res <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts,
    min_detected = 5L
  )
  expect_true("unmeasured" %in% res$pathway)
  expect_identical(res$n[res$pathway == "unmeasured"], 5L)
})

test_that("pathway_fry returns no rows when every set is under the floor", {
  d <- planted()
  res <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts,
    min_detected = 50L
  )
  expect_identical(nrow(res), 0L)
  expect_named(res, c("contrast", "pathway", "n", "direction", "p", "fdr"),
    ignore.order = TRUE
  )
})

test_that("pathway_fry applies BH within each contrast, never pooled", {
  d <- planted()
  cm <- limma::makeContrasts(BvA = B - A, AvB = A - B, levels = d$design)
  res <- pathway_fry(d$expr, d$gene_sets, d$design, cm)

  expect_setequal(unique(res$contrast), c("BvA", "AvB"))
  by_contrast <- split(res, res$contrast)
  for (r in by_contrast) {
    expect_equal(r$fdr, p.adjust(r$p, method = "BH"), tolerance = 1e-12)
  }
})

test_that("pathway_fry carries subject blocking through to fry", {
  d <- planted()
  block <- rep(seq_len(10), each = 2)
  cor_est <- limma::duplicateCorrelation(
    d$expr, d$design,
    block = block
  )$consensus
  res <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts,
    block = block, correlation = cor_est
  )
  unblocked <- pathway_fry(d$expr, d$gene_sets, d$design, d$contrasts)

  expect_setequal(res$pathway, unblocked$pathway)
  expect_false(isTRUE(all.equal(res$p, unblocked$p)))
})
