# The gene-set tree both pathway steps (02_HRvLR/02_Set_Contrasts, 03_Continuous/02_Trait_Pathways)
# use to group overlapping sets, and the rule that picks each cluster's representative. Expects
# dplyr attached.

# Candidate sets `d` (one block or window) grouped by the proteins they share: Jaccard distance
# (1 - shared / union; dist(method = "binary") on the set-by-gene incidence), average linkage, as
# EnrichmentMap and DAVID group them, cut at three heights. `sets` maps set_id to its genes. Adds
# cut_0.75, cut_0.85, cut_0.9 and the tree's row order.
set_tree <- function(d, sets) {
  if (nrow(d) == 1) {
    return(mutate(d, cut_0.75 = 1L, cut_0.85 = 1L, cut_0.9 = 1L, tree_order = 1L))
  }
  genes <- unique(unlist(sets[d$set_id]))
  incidence <- t(vapply(unname(sets[d$set_id]), \(s) genes %in% s, logical(length(genes))))
  tree <- hclust(dist(incidence, method = "binary"), method = "average")
  mutate(d,
    cut_0.75 = cutree(tree, h = 0.75), cut_0.85 = cutree(tree, h = 0.85),
    cut_0.9 = cutree(tree, h = 0.9), tree_order = order(tree$order)
  )
}

# Each cluster's representative: the smallest set whose p is within 10-fold of the cluster's best,
# so a narrow set nested in a broad one takes its place when the evidence is close.
most_specific <- function(cluster, p, size) {
  near_best <- p <= 10 * ave(p, cluster, FUN = min)
  pick <- order(cluster, !near_best, size, p)
  seq_along(p) %in% pick[!duplicated(cluster[pick])]
}
