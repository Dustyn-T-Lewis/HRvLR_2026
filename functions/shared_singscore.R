# Single-sample rank-based pathway scoring. Per-sample ranks make this
# leakage-free: a sample's score depends only on its own gene ranks.
score_singscore <- function(expr, gene_sets, min_size = 1L) {
  pacman::p_load(singscore)
  sizes <- vapply(gene_sets, function(g) sum(g %in% rownames(expr)), integer(1))
  gene_sets <- gene_sets[sizes >= min_size]
  ranks <- singscore::rankGenes(expr)
  scores <- do.call(rbind, lapply(
    gene_sets,
    function(g) singscore::simpleScore(ranks, upSet = g)$TotalScore
  ))
  rownames(scores) <- names(gene_sets)
  colnames(scores) <- colnames(expr)
  scores
}
