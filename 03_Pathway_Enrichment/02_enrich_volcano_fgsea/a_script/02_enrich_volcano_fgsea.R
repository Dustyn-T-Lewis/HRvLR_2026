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
out <- file.path(stage, "c_data")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_run_fgsea_and_fry first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
fg <- readRDS(paths[["set_tests"]])

# Bound only the logarithm for plotting. A protein untested in a contrast has no p there and is
# not drawn.
protein_results <- fg$protein_results |>
  filter(!is.na(p)) |>
  mutate(plot_p = pmax(p, .Machine$double.xmin))
fgsea_results <- filter(fg$set_tests, method == "fgsea", main)
stopifnot(nrow(fgsea_results) > 0, !is.null(protein_results$label))

# volcano_ring() matches leading-edge genes against the point labels, and a label carries its
# accession when a symbol sits on more than one protein. Translate the edges into label space
# or the tick lines silently draw nothing.
# Built from every protein, tested or not, so a protein untested in one contrast still gets its
# tick in another.
gene_to_label <- fg$protein_results |>
  filter(selected, !is.na(gene)) |>
  distinct(gene, label) |>
  with(set_names(label, gene))
fgsea_results$leading_edge <- map(fgsea_results$leading_edge, \(genes) {
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
      "Table: c_data/02_enrich_volcano_fgsea.xlsx."
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
# The primary and secondary questions first, then the four within-arm responses. The last two
# repeat Training_HR and Training_LR with pi-ranked labels: same points, colours and rings. A pi
# label ranks and selects nothing; a protein named on these two panels is not a hit.
plot_order <- c(
  "Training_Interaction", "Acute_Interaction", "Training_HR", "Training_LR", "Acute_HR",
  "Acute_LR"
)
volcanoes <- c(
  map(plot_order, make_volcano),
  map(c("Training_HR", "Training_LR"), make_volcano, rank_by = "pi")
)
pdf(file.path(figure_dir, "02_enrich_volcano_fgsea_figures.pdf"), width = 11, height = 8.5)
walk(volcanoes, \(v) print(v$figure))
invisible(dev.off())

sheets <- list(
  ringed_sets = list_rbind(map(volcanoes, "ring")),
  labelled_proteins = list_rbind(map(volcanoes, "labels")),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "The fgsea sets ringed on each volcano. Points come from set_tests.rds, protein_results.",
  "The proteins named on each volcano, by FDR or by pi score.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
writexl::write_xlsx(
  c(list(read_me = read_me), sheets), file.path(out, "02_enrich_volcano_fgsea.xlsx")
)
message("drew ", length(volcanoes), " volcanoes")
