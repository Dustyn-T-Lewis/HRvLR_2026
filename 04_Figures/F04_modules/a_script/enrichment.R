# What the modules are (over-representation), and whether they move (limma::fry).
pacman::p_load(here, dplyr, tibble, purrr, limma)

# Over-representation of each module against the measured proteome. The universe is every
# protein in the network, not the whole genome - anything not quantified could never have
# been assigned to a module.
module_ora <- function(colors, top_n = 5) {
  # Namespace-qualified, never attached: org.Hs.eg.db brings AnnotationDbi::select onto the
  # search path, which masks dplyr::select for every panel downstream.
  stopifnot(
    requireNamespace("clusterProfiler", quietly = TRUE),
    requireNamespace("org.Hs.eg.db", quietly = TRUE)
  )
  universe <- names(colors)
  modules <- setdiff(unique(colors), "grey")

  purrr::map_dfr(modules, function(m) {
    genes <- names(colors)[colors == m]
    res <- clusterProfiler::enrichGO(
      gene = genes, universe = universe, keyType = "SYMBOL",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db, ont = "BP",
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2
    )
    if (is.null(res) || nrow(as.data.frame(res)) == 0) {
      return(tibble(module = m, term = NA_character_, padj = NA_real_, count = NA_integer_))
    }
    as.data.frame(res) |>
      dplyr::slice_min(p.adjust, n = top_n, with_ties = FALSE) |>
      dplyr::transmute(module = m, term = Description, padj = p.adjust, count = Count) |>
      tibble::as_tibble()
  })
}

# Does a module move in a contrast?
#
# fry is a rotation-based SELF-CONTAINED gene-set test: it asks whether the proteins in this
# predefined set move together in this contrast, against the same linear model that produced
# the DEP results, with the subject blocking and duplicateCorrelation carried through.
#
# The panel this replaces ran fgsea with the WGCNA modules as the gene sets, ranked by
# moderated t from the SAME matrix that defined those modules. That is circular - the sets
# were built to be co-expressed, so of course they rank together - and it returned
# padj = 7e-33 in a study whose minimum q elsewhere is 0.36. fry is the correct test and is
# already used correctly in this repo (04_Figures/functions/f00_concordance_panels.R).
module_fry <- function(colors, abund, design, contrast_matrix, block, correlation,
                       contrasts = colnames(contrast_matrix)) {
  shared <- intersect(names(colors), rownames(abund))
  colors <- colors[shared]
  mat <- abund[shared, , drop = FALSE]

  index <- split(seq_along(shared), colors[shared])
  index <- index[setdiff(names(index), "grey")]
  index <- index[vapply(index, length, integer(1)) >= 3L]

  purrr::map_dfr(contrasts, function(ct) {
    res <- limma::fry(
      mat,
      index = index, design = design, contrast = contrast_matrix[, ct],
      block = block, correlation = correlation
    )
    res |>
      tibble::rownames_to_column("module") |>
      dplyr::transmute(
        contrast = ct, module,
        n = NGenes, direction = Direction,
        p = PValue, fdr = FDR
      ) |>
      tibble::as_tibble()
  })
}
