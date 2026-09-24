# Two set tests, reported side by side.
#
# fgsea is competitive: it ranks proteins by moderated t and asks whether a set piles up at one
# end, relative to every other protein. That assumes the ranked proteins are exchangeable, and
# co-regulated sets break the assumption.
#
# fry is self-contained: it asks whether a set moved at all, rotating the residuals of the fitted
# model rather than shuffling gene labels. It takes the design, the subject block and the
# within-subject correlation, so the repeated measures are built in. fry cannot take a missing
# value, so it reads the imputed matrix, with a correlation estimated on that same matrix.
#
# Both ship as rows of one table. Each method's BH runs over all 1,378 sets within a contrast,
# pooling the five collections, as BFR does.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(limma)
  library(ggplot2)
})

out <- here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  fit = "02_Differential_Expression/02_Differential/c_data/fit.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run the upstream stages first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
gs <- readRDS(paths[["gene_sets"]])
fit <- readRDS(paths[["fit"]])
d <- readRDS(paths[["design"]])
imputed <- readRDS(paths[["imputed"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

sets <- gs$sets
protein_map <- gs$protein_map
protein_ids <- protein_map$protein
eb <- fit$eBayes_fit
contrast_names <- colnames(fit$design$contrast_matrix)
design <- fit$design$design_matrix
subject <- fit$metadata$subject
# A stale fit still lines up by count, so these compare names rather than lengths.
stopifnot(
  identical(rownames(eb$coefficients), protein_ids),
  identical(colnames(eb$coefficients), contrast_names),
  identical(rownames(imputed$data), protein_ids),
  identical(colnames(imputed$data), rownames(design)),
  d$floor %in% contrast_names
)

# topTable rebuilds all nine contrasts from the saved fit, BH within each as upstream applied it.
# pi_score is the Xiao et al. (2014) score: it orders a contrast and selects nothing.
protein_results <- map(set_names(contrast_names), function(contrast) {
  topTable(eb, coef = contrast, number = Inf, adjust.method = "BH", sort.by = "none") |>
    rownames_to_column("protein") |>
    select(protein, logFC, AveExpr, t, P.Value, adj.P.Val)
}) |>
  list_rbind(names_to = "contrast") |>
  left_join(protein_map, by = "protein", relationship = "many-to-one") |>
  mutate(pi_score = P.Value^abs(logFC))
stopifnot(nrow(protein_results) == length(protein_ids) * length(contrast_names))

protein_summary <- protein_results |>
  summarise(
    tested = sum(!is.na(P.Value)),
    up = sum(adj.P.Val < 0.05 & logFC > 0, na.rm = TRUE),
    down = sum(adj.P.Val < 0.05 & logFC < 0, na.rm = TRUE),
    min_fdr = signif(min(adj.P.Val, na.rm = TRUE), 3), .by = contrast
  )
print(protein_summary)

# One vector per contrast, representative proteins only, keyed by gene symbol so fgsea's
# leadingEdge comes back in the namespace the volcanoes label with. A protein untested in a
# contrast has no t and drops out of that ranking.
ranked <- protein_results |>
  filter(selected, !is.na(gene), !is.na(t)) |>
  split(~contrast) |>
  map(\(result) set_names(result$t, result$gene))
ranked <- ranked[contrast_names]

# fry indexes rows of the matrix, not genes, so the sets are mapped back through the
# representative protein chosen in 00_build_gene_sets.
gene_map <- filter(protein_map, selected)
set_rows <- map(sets, \(genes) match(gene_map$protein[match(genes, gene_map$gene)], protein_ids))
stopifnot(!any(map_lgl(set_rows, anyNA)))

abundance <- as.matrix(imputed$data)
correlation <- duplicateCorrelation(abundance, design, block = subject)$consensus.correlation
message("within-subject correlation on the imputed matrix: ", round(correlation, 3))

set.seed(1)
fgsea_raw <- map(ranked, \(stats) fgsea::fgsea(sets, stats, minSize = 15, maxSize = 500))

# Both packages return their own BH column over the list they were given: fgsea's padj and
# fry's FDR. Neither is recomputed here.
set_tests <- bind_rows(
  map(contrast_names, \(contrast) {
    as_tibble(fgsea_raw[[contrast]]) |>
      transmute(
        set_id = pathway, contrast, method = "fgsea", n = size,
        direction = if_else(NES > 0, "Up", "Down"), NES, p = pval, padj,
        leadingEdge
      )
  }),
  map(contrast_names, \(contrast) {
    fry(abundance, set_rows, design, fit$design$contrast_matrix[, contrast],
      block = subject, correlation = correlation, sort = "none"
    ) |>
      rownames_to_column("set_id") |>
      transmute(
        set_id, contrast,
        method = "fry", n = NGenes,
        direction = Direction, NES = NA_real_, p = PValue, padj = FDR,
        leadingEdge = list(NULL)
      )
  })
)

# Five collections overlap, so glycolysis is tested in Hallmark, KEGG, Reactome and GO.
# collapsePathways re-runs each significant set conditioned on a more significant one's leading
# edge and keeps it only if it stays significant. It prunes after testing rather than
# re-adjusting, so the surviving list's FDR is conservative.
main_sets <- map(set_names(contrast_names), function(contrast) {
  significant <- fgsea_raw[[contrast]][padj < 0.05][order(pval)]
  if (nrow(significant) < 2) {
    return(significant$pathway)
  }
  fgsea::collapsePathways(significant, sets, ranked[[contrast]], pval.threshold = 0.05)$mainPathways
})

main_lookup <- imap(main_sets, \(ids, cn) tibble(contrast = cn, set_id = ids, kept = TRUE)) |>
  list_rbind()
set_tests <- set_tests |>
  left_join(main_lookup, by = c("contrast", "set_id")) |>
  mutate(main = if_else(method == "fgsea", coalesce(kept, FALSE), NA), kept = NULL) |>
  left_join(
    select(gs$set_catalog, set_id, database, pathway, source_size, description),
    by = "set_id"
  ) |>
  relocate(contrast, method, set_id, database, pathway)

# One row per contrast: how many sets each test called, and how many survived collapse. The floor
# is printed first.
set_summary <- set_tests |>
  summarise(
    sets = n_distinct(set_id),
    fgsea = sum(method == "fgsea" & padj < 0.05),
    collapsed = sum(method == "fgsea" & padj < 0.05 & main),
    fry = sum(method == "fry" & padj < 0.05),
    .by = contrast
  ) |>
  arrange(contrast != d$floor)
print(as.data.frame(set_summary))


# ---- figures -------------------------------------------------------------------------------

# One directory per collection plus all_db, one file per contrast. Each panel shows the ten
# strongest collapse survivors by adjusted p.
figure_root <- here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "b_reports")
# Clear last run's figures so the bundle holds only this run's pages.
unlink(list.files(figure_root, "[.](png|pdf)$", full.names = TRUE, recursive = TRUE))
save_figure <- function(figure, file, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(paste0(file, ".", extension), figure,
      width = width, height = height, dpi = 200, bg = "white"
    )
  })
}
shown <- c(
  "Training_Interaction", "Acute_Interaction", "Training_HR", "Training_LR", "Acute_HR",
  "Acute_LR"
)
survivors <- set_tests |>
  filter(method == "fgsea", padj < 0.05, main, contrast %in% shown)

draw_dotplot <- function(rows, colour_by, file) {
  top <- rows |>
    slice_min(padj, n = 10, with_ties = FALSE) |>
    mutate(label = enrichVolcano::ev_clean_label(pathway))
  figure <- ggplot(top, aes(NES, reorder(label, NES), size = n, colour = .data[[colour_by]])) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2, 6), name = "genes") +
    labs(
      x = "normalised enrichment score", y = NULL, title = unique(top$contrast),
      subtitle = sprintf(
        "%s, fgsea, BH within contrast",
        if (colour_by == "database") "all collections" else rows$database[1]
      ),
      caption = sprintf(
        paste(
          "%s of %d collapse survivor%s, ranked by adjusted p.",
          "Size is gene count, colour is -log10 FDR. Table: c_data/set_tests.csv."
        ),
        if (nrow(rows) > 10) "Ten strongest" else "All", nrow(rows),
        if (nrow(rows) == 1) "" else "s"
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.caption = element_text(size = 6.5, colour = "grey45", hjust = 0)
    )
  figure <- if (colour_by == "database") {
    figure + scale_colour_brewer(palette = "Dark2", name = NULL)
  } else {
    figure + scale_colour_viridis_c(
      option = "rocket", direction = -1, end = 0.9,
      name = expression(-log[10] ~ FDR)
    )
  }
  save_figure(figure, file, width = 7.5, height = max(3, 1.9 + 0.38 * nrow(top)))
}

collections <- unique(gs$set_catalog$database[gs$set_catalog$qualifies])
for (db in c(collections, "all_db")) {
  dir <- file.path(figure_root, db)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rows <- if (db == "all_db") survivors else filter(survivors, database == db)
  drawn <- intersect(shown, unique(rows$contrast))
  for (cn in drawn) {
    draw_dotplot(
      mutate(filter(rows, contrast == cn), `-log10 FDR` = -log10(padj)),
      if (db == "all_db") "database" else "-log10 FDR",
      file.path(dir, paste0("01_dotplot_", cn))
    )
  }
  message("drew ", db, ": ", length(drawn), " contrasts")
}

collapse_effect <- set_tests |>
  filter(method == "fgsea", padj < 0.05, contrast %in% shown) |>
  summarise(before = n(), after = sum(main), .by = c(contrast, database)) |>
  tidyr::pivot_longer(c(before, after), names_to = "stage", values_to = "sets") |>
  mutate(
    stage = factor(stage, c("before", "after")),
    contrast = factor(contrast, levels = shown)
  )
save_figure(
  ggplot(collapse_effect, aes(stage, sets, fill = database)) +
    geom_col(position = "dodge") +
    facet_wrap(~contrast, scales = "free_y", nrow = 1) +
    scale_fill_brewer(palette = "Dark2", name = NULL) +
    labs(
      x = NULL, y = "significant sets", title = "Significant sets before and after collapse",
      subtitle = "fgsea at FDR 0.05, then collapsePathways",
      caption = paste(
        "collapsePathways re-tests each significant set conditioned on a stronger set's leading",
        "edge and keeps it only if it stays significant. Table: c_data/set_tests.csv."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)),
  file.path(figure_root, "all_db", "02_collapse_before_after"),
  width = 12, height = 3
)

# Every set nominal under fry in at least one contrast, as a dot matrix. Fill is the fgsea NES,
# for direction; size is -log10 fry p; a black ring marks fry FDR < 0.05.
# One label per set, built before the join so two sets sharing a display name stay apart.
set_labels <- gs$set_catalog |>
  filter(qualifies) |>
  transmute(
    set_id,
    label = gsub("\n", " ", paste0("[", database, "] ", enrichVolcano::ev_clean_label(pathway)))
  ) |>
  mutate(label = if_else(duplicated(label) | duplicated(label, fromLast = TRUE), set_id, label))
hit_data <- set_tests |>
  filter(method == "fry") |>
  select(set_id, contrast, p, fdr = padj) |>
  left_join(
    set_tests |> filter(method == "fgsea") |> select(set_id, contrast, effect = NES),
    by = c("set_id", "contrast")
  ) |>
  left_join(set_labels, by = "set_id")
ranked_hits <- hit_data |>
  summarise(n_nominal = sum(p < 0.05), best = min(p), .by = label) |>
  filter(n_nominal > 0) |>
  arrange(desc(n_nominal), best)
page_of <- ceiling(seq_len(nrow(ranked_hits)) / 75)
hit_plots <- list()
for (page in seq_len(max(page_of, 0))) {
  rows <- ranked_hits$label[page_of == page]
  page_data <- hit_data |>
    filter(label %in% rows) |>
    mutate(
      label = factor(label, levels = rev(rows)),
      contrast = factor(contrast, levels = contrast_names), nominal = p < 0.05
    )
  hit_plots[[page]] <-
    ggplot(page_data, aes(contrast, label)) +
    geom_point(data = filter(page_data, !nominal), colour = "grey85", size = 0.6) +
    geom_point(
      data = filter(page_data, nominal), aes(size = -log10(p), fill = effect),
      shape = 21, colour = "grey30", stroke = 0.2
    ) +
    geom_point(
      data = filter(page_data, fdr < 0.05), aes(size = -log10(p)),
      shape = 21, colour = "black", stroke = 0.9, fill = NA
    ) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    scale_size_continuous(range = c(1, 4.5), limits = c(-log10(0.05), NA)) +
    labs(
      title = sprintf("Set contrast hits, fry (%d/%d)", page, max(page_of)),
      subtitle = sprintf(
        "fry p per contrast, filled by fgsea NES; %d sets nominal in at least one contrast",
        nrow(ranked_hits)
      ),
      x = NULL, y = NULL, fill = "NES", size = "-log10 p",
      caption = paste(
        "Rows: sets at fry p < 0.05 in at least one contrast, most contrasts first, then by",
        "best p. Filled dot: nominal, fill = fgsea NES, size = -log10 fry p. Black ring: fry",
        "FDR < 0.05. Grey speck: not nominal. Table: c_data/set_tests.csv."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = if (length(rows) > 30) 5.5 else 8),
      plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)
    )
}
# Paged, so one PDF and no PNG.
dir.create(file.path(figure_root, "hits"), showWarnings = FALSE)
pdf(file.path(figure_root, "hits", "03_set_hits.pdf"), width = 11, height = 8.5, bg = "white")
walk(hit_plots, print)
invisible(dev.off())
message("drew ", length(hit_plots), " hit-matrix pages")

packages <- c("here", "limma", "fgsea", "dplyr", "purrr", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
flat <- mutate(set_tests, leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

saveRDS(
  list(
    protein_results = protein_results, set_tests = set_tests, set_summary = set_summary,
    correlation = correlation,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "set_tests.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    set_summary = set_summary,
    significant = filter(flat, padj < 0.05),
    protein_summary = protein_summary,
    protein_results = protein_results,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "01_run_fgsea_and_fry.xlsx")
)
readr::write_csv(flat, file.path(out, "set_tests.csv"))
combined <- file.path(figure_root, "01_run_fgsea_and_fry_figures.pdf")
pages <- setdiff(list.files(figure_root, "[.]pdf$", recursive = TRUE, full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote set_tests.rds, 01_run_fgsea_and_fry.xlsx, set_tests.csv and ",
  length(pages), " figures bundled into ", qpdf::pdf_length(combined), " pages"
)
