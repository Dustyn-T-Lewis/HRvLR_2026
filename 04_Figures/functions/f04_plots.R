# Plotting helpers shared by the F04 association leaves: an effect-size
# volcano for the limma/HLM association tables, and a thin wrapper over
# singscore's own rank-density diagnostic.
pacman::p_load(ggplot2, dplyr, singscore, ggrepel)

# Effect-size volcano for one trait: slope vs -log10 BH, top features labelled.
assoc_volcano <- function(res, trait_label = unique(res$trait), bh_cut = 0.05,
                          top_n = 12) {
  res <- dplyr::mutate(res,
    sig = bh < bh_cut,
    nlog_bh = -log10(pmax(bh, .Machine$double.xmin))
  )
  lab <- dplyr::slice_min(dplyr::filter(res, sig), bh, n = top_n)
  ggplot(res, aes(beta, nlog_bh)) +
    geom_hline(yintercept = -log10(bh_cut), linetype = 2, colour = "grey60") +
    geom_point(aes(colour = sig), alpha = 0.7, show.legend = FALSE) +
    ggrepel::geom_text_repel(
      data = lab, aes(label = feature), size = 2.6,
      max.overlaps = 20, seed = 42
    ) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey75")) +
    labs(
      x = "association slope (per unit trait)", y = expression(-log[10] ~ BH),
      title = trait_label
    ) +
    theme_bw(base_size = 10)
}

# singscore's own rank-density plot for one sample against one gene set.
singscore_rank_density <- function(mat, gene_set, sample_id) {
  ranked <- singscore::rankGenes(mat)
  one <- ranked[, sample_id, drop = FALSE]
  singscore::plotRankDensity(one, upSet = gene_set)
}
