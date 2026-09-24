# Draw the protein volcanoes with the surviving fgsea pathways ringed. Computes nothing.
# Two separate claims share a panel: point colour and the count badges read protein-level BH
# FDR, the ring reads set-level fgsea FDR. Only sets that survived collapsePathways are ringed.
# The floor contrast is not drawn; 01_run_fgsea_and_fry reports what it returns.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

stage <- here("03_Pathway_Enrichment", "02_enrich_volcano_fgsea")
figure_dir <- file.path(stage, "b_reports")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
# Clear last run's figures so the bundle holds only this run's pages.
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

inputs <- c(set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_run_fgsea_and_fry first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
fg <- readRDS(paths[["set_tests"]])

# Keep exact p-values in the export; bound only the logarithm for plotting.
# A protein untested in a contrast has no p there and is not drawn.
protein_results <- fg$protein_results |>
  filter(!is.na(P.Value)) |>
  mutate(plot_p = pmax(P.Value, .Machine$double.xmin))
# the results table is long, one row per set per contrast per method; rings read fgsea
fgsea_results <- filter(fg$set_tests, method == "fgsea", main)
stopifnot(nrow(fgsea_results) > 0, !is.null(protein_results$label))

# volcano_ring() matches leading-edge genes against the point labels, and a label carries its
# accession when a symbol sits on more than one protein. Translate the edges into label space
# or the tick lines silently draw nothing.
# Built from every protein, tested or not, so a protein untested in one contrast still gets its
# tick in another.
gene_to_label <- fg$protein_results |>
  filter(!is.na(gene)) |>
  distinct(gene, label) |>
  with(set_names(label, gene))
fgsea_results$leadingEdge <- map(fgsea_results$leadingEdge, \(genes) {
  unname(gene_to_label[genes[genes %in% names(gene_to_label)]])
})

# The default palette is red for up, blue for down, with a dark blue to dark red NES ramp.
plot_theme <- enrichVolcano::volcano_ring_theme(
  base_size = 12, base_family = "sans", palette = "default", ns = "#9a9a9a"
)
# The two interactions are differences of differences, so they carry their algebra as a
# subtitle; the other contrast names read on their own.
contrast_subtitle <- c(
  Training_Interaction = "(HR_T2 - HR_T1) - (LR_T2 - LR_T1)",
  Acute_Interaction = "(HR_T3 - HR_T2) - (LR_T3 - LR_T2)"
)

# A named vector indexed by a missing name returns NA, not NULL, so the subtitle lookup is
# dropped with na.omit() rather than defaulted with %||%.
make_volcano <- function(contrast, rank_by = "fdr") {
  points <- filter(protein_results, .data$contrast == .env$contrast)
  ring <- fgsea_results |>
    filter(.data$contrast == .env$contrast, padj < 0.05) |>
    slice_min(padj, n = 8) |>
    # Overlapping databases mean two ringed sets can carry the same name, and
    # volcano_ring() strips the database prefix before drawing. GOBP_MUSCLE_CONTRACTION and
    # REACTOME_MUSCLE_CONTRACTION then label two arcs identically.
    mutate(
      stem = sub("^[A-Z0-9]+_", "", pathway),
      pathway = if_else(
        duplicated(stem) | duplicated(stem, fromLast = TRUE),
        paste0(pathway, " ", database), pathway
      )
    )
  if (rank_by == "pi") {
    labels <- points |>
      arrange(pi_score, protein) |>
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
      filter(adj.P.Val < 0.05) |>
      arrange(adj.P.Val, P.Value, protein) |>
      slice_head(n = 5) |>
      pull(label)
    subtitle <- paste(
      c(na.omit(contrast_subtitle[contrast]), "Protein significance: BH FDR < 0.05"),
      collapse = "\n"
    )
  }
  enrichVolcano::volcano_ring(
    volc_df = select(points, label, logFC, plot_p, adj.P.Val),
    enrich_df = ring,
    gene_col = "label", pval_col = "plot_p", padj_col = "padj",
    volc_sig_col = "adj.P.Val", genes_col = "leadingEdge",
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
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "plain"),
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 9, hjust = 0.5),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(22, "mm")
    )
}
save_volcano <- function(volcano, name) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)),
      plot = volcano, width = 7, height = 6.5, units = "in", dpi = 300, bg = "white"
    )
  })
}

# The primary and secondary questions first, then the four within-arm responses. The floor is
# fitted and reported upstream, not drawn here.
plot_order <- c(
  "Training_Interaction", "Acute_Interaction", "Training_HR", "Training_LR", "Acute_HR",
  "Acute_LR"
)
for (contrast in plot_order) {
  save_volcano(make_volcano(contrast), paste0("protein_volcano_fdr_", contrast))
  message("drew ", contrast)
}
# Same points, colours, rings and notation; only the labels move. A pi label ranks and selects
# nothing; a protein named on these two panels is not a hit.
for (contrast in c("Training_HR", "Training_LR")) {
  volcano <- make_volcano(contrast, rank_by = "pi")
  save_volcano(volcano, paste0("protein_volcano_pi_rank_", contrast))
  message("drew ", contrast, " ranked by pi")
}

drawn <- sort(basename(list.files(figure_dir, pattern = "[.]png$")))
print(tibble(file = drawn))
stopifnot(length(drawn) == length(plot_order) + 2)

combined <- file.path(figure_dir, "02_enrich_volcano_fgsea_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(length(pages), " figures bundled into ", qpdf::pdf_length(combined), " pages")
