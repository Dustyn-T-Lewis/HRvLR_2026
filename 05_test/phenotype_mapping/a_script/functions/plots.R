# Figures for the phenotype-association unit. Two custom ggplots for the
# association result, plus thin wrappers over singscore's own diagnostic plots
# so the scoring step is shown, not just consumed.

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
      max.overlaps = 20
    ) +
    scale_colour_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey75")) +
    labs(
      x = "association slope (per unit trait)", y = expression(-log[10] ~ BH),
      title = trait_label
    ) +
    theme_bw(base_size = 10)
}

# Two pathway scores per sample as a 2D map, coloured by a sample attribute.
pathway_scatter <- function(scores, meta, set_x, set_y, colour_by = "Group",
                            id_col = "Col_ID") {
  df <- data.frame(
    id = colnames(scores), x = scores[set_x, ], y = scores[set_y, ]
  )
  df[[colour_by]] <- meta[[colour_by]][match(df$id, meta[[id_col]])]
  ggplot(df, aes(x, y, colour = .data[[colour_by]])) +
    geom_point(size = 2) +
    labs(x = set_x, y = set_y) +
    theme_bw(base_size = 10)
}

# singscore's own rank-density plot for one sample against one gene set.
singscore_rank_density <- function(mat, gene_set, sample_id) {
  ranked <- singscore::rankGenes(mat)
  one <- ranked[, sample_id, drop = FALSE]
  singscore::plotRankDensity(one, upSet = gene_set)
}

# singscore's score-vs-dispersion diagnostic across samples for one gene set.
singscore_dispersion <- function(mat, gene_set) {
  ranked <- singscore::rankGenes(mat)
  scored <- singscore::simpleScore(ranked, upSet = gene_set)
  singscore::plotDispersion(scored, annot = colnames(mat))
}
