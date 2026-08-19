# F03 supplement: hierarchical pathway trees, one per contrast with enriched
# terms Terms are clustered on Jaccard overlap of their member genes, so
# branches group sets that share proteins rather than sets that merely score
# alike. Tip colour is NES, tip size is -log10 padj.
#
# These are fgsea terms. METHOD_RANKING ranks fgsea below fry on null validity
# and keeps it for ranking and display only, and limma::fry returns zero over
# the same sets in every HR-vs-LR and interaction contrast. Each panel states
# which regime it is in.

pacman::p_load(
  here, dplyr, tidyr, readxl, ggplot2, ggtree, ape, patchwork, stringr
)

source(here("functions", "shared_style.R"))
source(here("functions", "shared_pathway_utils.R"))
source(here("functions", "shared_utils.R"))

TREE_TOP_N <- 25
TREE_K <- 5
GROUP_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

RPT_DIR <- here("04_Figures", "F03_pathway", "supp", "treeplots", "b_reports")
DAT_DIR <- here("04_Figures", "F03_pathway", "supp", "treeplots", "c_data")
clear_dir(RPT_DIR)
dir.create(DAT_DIR, recursive = TRUE, showWarnings = FALSE)

# fgsea_significant carries truncated display names; fgsea_all keeps the set ids
# the gene-set collection is keyed by, so the overlap join has something to
# match on.
fg <- read_excel(
  here("04_Figures", "F03_pathway", "c_data", "F03_pathway_source_data.xlsx"),
  sheet = "fgsea_all"
) |>
  filter(.data$padj < 0.05)
collection <- build_pathway_collection(
  min_size = SET_FLOOR, max_size = 500, include_goslim = TRUE,
  exclude_variants = TRUE
)

jaccard <- function(sets) {
  n <- length(sets)
  m <- matrix(0, n, n, dimnames = list(names(sets), names(sets)))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      m[i, j] <- length(intersect(sets[[i]], sets[[j]])) /
        length(union(sets[[i]], sets[[j]]))
    }
  }
  m
}

tidy_label <- function(x) {
  x |>
    str_remove("^(HALLMARK|GOBP|GOCC|GOMF|GOSLIM|REACTOME|KEGG_MEDICUS)_") |>
    str_replace_all("_", " ") |>
    str_to_sentence() |>
    str_trunc(46)
}

build_tree <- function(contrast) {
  d <- fg |>
    filter(
      .data$contrast == !!contrast, .data$pathway %in% names(collection)
    ) |>
    slice_min(.data$padj, n = TREE_TOP_N, with_ties = FALSE)
  if (nrow(d) < 4) {
    return(NULL)
  }

  sets <- collection[d$pathway]
  hc <- hclust(as.dist(1 - jaccard(sets)), method = "average")
  phy <- ape::as.phylo(hc)

  # %<+% matches the frame's first column against the tree's tip labels, so
  # label leads and nothing else is carried that ggtree already defines.
  tips <- d |>
    transmute(
      label = .data$pathway,
      display = tidy_label(.data$pathway),
      nes = .data$NES,
      weight = -log10(.data$padj)
    )

  gated <- contrast %in% GROUP_CONTRASTS
  ggtree(phy, size = 0.3) %<+% tips +
    geom_tippoint(aes(colour = nes, size = weight)) +
    geom_tiplab(
      aes(label = display),
      size = 1.7, offset = 0.02, colour = "grey20"
    ) +
    scale_colour_gradient2(
      low = DIR_COLORS[["Down"]], mid = "grey85", high = DIR_COLORS[["Up"]],
      midpoint = 0, name = "NES"
    ) +
    scale_size_continuous(range = c(0.7, 2.6), name = "-log10 padj") +
    labs(
      title = sprintf(
        "%s  |  %d terms shown of %d", contrast, nrow(d),
        sum(fg$contrast == contrast)
      ),
      subtitle = if (gated) {
        "limma::fry returns zero over these sets. fgsea is display-only here."
      } else {
        "Within-arm contrast. T3 samples carry roughly twice the blood of T1 and T2."
      }
    ) +
    xlim(0, 1.55) +
    theme_tree() +
    theme(
      plot.title = element_text(face = "bold", size = 6.5),
      plot.subtitle = element_text(
        size = 5.4, colour = if (gated) "#B2182B" else "grey40"
      ),
      legend.position = "none"
    )
}

contrasts_with_terms <- fg |>
  count(.data$contrast) |>
  filter(.data$n >= 4) |>
  arrange(desc(.data$n)) |>
  pull(.data$contrast)

trees <- Filter(Negate(is.null), stats::setNames(
  lapply(contrasts_with_terms, build_tree), contrasts_with_terms
))

sheet <- wrap_plots(trees, ncol = 2) +
  plot_annotation(
    title = "Pathway trees: terms clustered on shared member proteins",
    caption = paste(
      "Branch length is 1 - Jaccard overlap of member genes, average linkage.",
      "Tip colour is NES, tip size is -log10 padj.\nRed subtitles mark the contrasts",
      "that answer the responder question, where fry finds nothing over the same sets."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.caption = element_text(hjust = 0, size = 6, colour = "grey35")
    )
  )

save_panel(
  sheet, file.path(RPT_DIR, "F03_supp_treeplots"), 210,
  54 * ceiling(length(trees) / 2) + 22
)

write.csv(
  fg |>
    filter(.data$contrast %in% names(trees)) |>
    arrange(.data$contrast, .data$padj),
  file.path(DAT_DIR, "treeplot_terms.csv"),
  row.names = FALSE
)
cat(sprintf("F03 treeplots done: %d contrasts\n", length(trees)))
