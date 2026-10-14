# Volcanoes with fgsea rings
#
# Each page is one contrast's protein volcano with the surviving fgsea pathways ringed around it.
# The figure computes nothing. Two separate claims share a panel: point colour and the count badges
# read protein-level BH FDR, the ring reads set-level fgsea FDR. Only sets that survived
# `collapsePathways` are ringed. The baseline reference is not drawn; `03_Set_Contrasts` reports
# what it returns.
#
# fgsea ranks a gene set against all other proteins and returns a normalised enrichment score (NES).
# A set's leading edge is the proteins that drive its score. The pi score is Π = p^|logFC| (Xiao et
# al. 2014): small when p is small and the fold change large.
#
# Packages
#
# `write_workbook()` comes from `00_Setup/Shared_functions/analysis_reporting_helpers.R`; the
# volcano itself comes from `enrichVolcano`.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})
source(here("00_Setup", "Shared_functions", "analysis_reporting_helpers.R"))

# Input
#
# The input is `set_tests.rds` from `03_Set_Contrasts`, which holds both the per-protein results and
# the set tests. The chunk creates `b_reports/` and `c_data/` and stops if the file is missing.

stage <- here("05_Figures", "HRvLR", "S01_Volcano")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(set_tests = "03_Group_comparisons/HRvLR/c_data/03_set_tests.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 03_Set_Contrasts first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
fg <- readRDS(paths[["set_tests"]])

# Points and rings
#
# `plot_p` bounds only the logarithm for plotting: a p of zero becomes the smallest positive double.
# The ring candidates are the main fgsea results.

protein_results <- mutate(fg$protein_results, plot_p = pmax(p, .Machine$double.xmin))
fgsea_results <- filter(fg$set_tests, method == "fgsea", main)
stopifnot(nrow(fgsea_results) > 0, !is.null(protein_results$label))

# `volcano_ring()` matches leading-edge genes against the point labels, and a label carries its
# accession when a symbol sits on more than one protein. The chunk translates the edges into label
# space; otherwise the tick lines silently draw nothing.

gene_to_label <- fg$protein_results |>
  filter(selected, !is.na(gene)) |>
  distinct(gene, label) |>
  with(set_names(label, gene))
fgsea_results$leading_edge <- map(fgsea_results$leading_edge, \(genes) {
  unname(gene_to_label[genes[genes %in% names(gene_to_label)]])
})

# Theme and subtitles
#
# The default palette is red for up, blue for down, with a dark blue to dark red NES ramp. The two
# interactions are differences of differences, so they carry their algebra as a subtitle; the other
# contrast names read on their own.

plot_theme <- enrichVolcano::volcano_ring_theme(
  base_size = 12, base_family = "sans", palette = "default", ns = "#9a9a9a"
)
contrast_subtitle <- c(
  Training_Interaction = "(HR_T2 - HR_T1) - (LR_T2 - LR_T1)",
  Acute_Interaction = "(HR_T3 - HR_T2) - (LR_T3 - LR_T2)"
)

# One volcano
#
# `make_volcano()` takes a contrast name and a label rule and returns the plot, the ringed sets and
# the labelled proteins. The ring holds up to eight fgsea sets at FDR < 0.05, smallest FDR first.
# With `rank_by = "fdr"` it labels the five proteins with the smallest FDR among those at FDR <
# 0.05; with `rank_by = "pi"` it labels the five with the smallest pi score.

make_volcano <- function(contrast, rank_by = "fdr") {
  points <- filter(protein_results, .data$contrast == .env$contrast)
  ring <- fgsea_results |>
    filter(.data$contrast == .env$contrast, fdr < 0.05) |>
    slice_min(fdr, n = 8, with_ties = FALSE) |>
    # Overlapping collections mean two ringed sets can carry the same name, and
    # volcano_ring() strips the collection prefix before drawing. GOBP_MUSCLE_CONTRACTION and
    # REACTOME_MUSCLE_CONTRACTION then label two arcs identically.
    mutate(
      stem = sub("^[A-Z0-9]+_", "", pathway),
      pathway = if_else(
        duplicated(stem) | duplicated(stem, fromLast = TRUE),
        paste0(pathway, " ", collection), pathway
      )
    )
  # Indexing a named vector by a missing name returns NA, not NULL, so the interaction algebra
  # is dropped with na.omit() rather than defaulted with %||%.
  if (rank_by == "pi") {
    labels <- points |>
      arrange(pi_score, uniprot_id) |>
      slice_head(n = 5) |>
      pull(label)
    subtitle <- paste(
      c(
        na.omit(contrast_subtitle[contrast]),
        "Labels ranked by pi-score; colours and counts use BH FDR"
      ),
      collapse = "\n"
    )
  } else {
    labels <- points |>
      filter(fdr < 0.05) |>
      arrange(fdr, p, uniprot_id) |>
      slice_head(n = 5) |>
      pull(label)
    subtitle <- paste(
      c(na.omit(contrast_subtitle[contrast]), "Protein significance: BH FDR < 0.05"),
      collapse = "\n"
    )
  }
  volcano <- enrichVolcano::volcano_ring(
    volc_df = select(points, label, logFC, plot_p, fdr),
    enrich_df = ring,
    gene_col = "label", pval_col = "plot_p", padj_col = "fdr", nes_col = "nes",
    size_col = "n_proteins", volc_sig_col = "fdr", genes_col = "leading_edge",
    p_threshold = 0.05, logfc_threshold = 0,
    title = contrast, subtitle = subtitle,
    label_mode = if (length(labels)) "by_genes" else "none",
    label_genes = labels, label_n = 5,
    label_size = 3.1, axis_size = 3.1, count_size = 3.3,
    point_size = 1.2, theme = plot_theme
  ) +
    # volcano_ring() draws with clip = "off" against a square panel and puts the NES key on
    # the right, so a label anchored on that edge lands on top of the colourbar. Below the
    # plot there is nothing to collide with.
    guides(fill = guide_colorbar(direction = "horizontal", title.position = "top")) +
    labs(caption = paste(
      "Point colour: protein BH FDR < 0.05. Ring: up to eight fgsea sets at FDR < 0.05 that",
      "survived collapsePathways, ticks on their leading-edge proteins.",
      "Table: c_data/s01_volcano.xlsx."
    )) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "plain"),
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 9, hjust = 0.5),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(22, "mm"),
      plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)
    )
  list(
    figure = volcano,
    ring = transmute(ring, contrast, rank_by, set_id, collection, pathway, nes, fdr),
    labels = tibble(contrast, rank_by, label = labels)
  )
}

# All volcanoes
#
# The primary and secondary questions come first, then the four within-group responses. The last two
# pages repeat Training_HR and Training_LR with pi-ranked labels: same points, colours and rings. A
# pi label ranks and selects nothing; a protein named on these two panels is not a hit. The chunk
# writes all eight pages to `b_reports/s01_volcano.pdf`.

plot_order <- c(
  "Training_Interaction", "Acute_Interaction", "Training_HR", "Training_LR", "Acute_HR",
  "Acute_LR"
)
volcanoes <- c(
  map(plot_order, make_volcano),
  map(c("Training_HR", "Training_LR"), make_volcano, rank_by = "pi")
)
pdf(file.path(figure_dir, "s01_volcano.pdf"), width = 11, height = 8.5)
walk(volcanoes, \(v) print(v$figure))
invisible(dev.off())

# Save
#
# `c_data/s01_volcano.xlsx` holds the ringed sets and the labelled proteins for every page, with the
# input manifest.

sheets <- list(
  ringed_sets = list_rbind(map(volcanoes, "ring")),
  labelled_proteins = list_rbind(map(volcanoes, "labels"))
)
holds <- c(
  "The fgsea sets ringed on each volcano. Points come from set_tests.rds, protein_results.",
  "The proteins named on each volcano, by FDR or by pi score."
)
write_workbook(sheets, holds, file.path(out, "s01_volcano.xlsx"), inputs, paths)
