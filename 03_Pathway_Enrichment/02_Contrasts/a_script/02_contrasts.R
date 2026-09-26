# fgsea and fry, side by side. fgsea is competitive on moderated t and
# assumes exchangeable proteins, which they are not. fry is self-contained and rotates residuals
# under the design, subject block and within-subject correlation. fry takes no missing value, so
# it reads the imputed matrix with a correlation estimated on that matrix. Each method's own BH
# runs within contrast over the sets it tested: all of them for fry, fewer for fgsea when a
# contrast leaves a set under minSize.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(limma)
  library(ggplot2)
})

out <- here("03_Pathway_Enrichment", "02_Contrasts", "c_data")
figure_dir <- here("03_Pathway_Enrichment", "02_Contrasts", "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  fit = "02_Differential_Expression/02_Contrasts/c_data/fit.rds",
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

sets <- gs$sets
protein_map <- gs$protein_map
protein_ids <- protein_map$uniprot_id
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
  identical(fit$metadata$sample_id, rownames(design)),
  d$floor %in% contrast_names
)

# topTable rebuilds all nine contrasts from the saved fit, BH within each as upstream applied it.
# pi_score is the Xiao et al. (2014) score: it orders a contrast and selects nothing.
protein_results <- map(set_names(contrast_names), function(contrast) {
  topTable(eb, coef = contrast, number = Inf, adjust.method = "BH", sort.by = "none") |>
    rownames_to_column("uniprot_id") |>
    select(uniprot_id, logFC, AveExpr, t, p = P.Value, fdr = adj.P.Val)
}) |>
  list_rbind(names_to = "contrast") |>
  left_join(protein_map, by = "uniprot_id", relationship = "many-to-one") |>
  mutate(pi_score = p^abs(logFC))
stopifnot(nrow(protein_results) == length(protein_ids) * length(contrast_names))

protein_summary <- protein_results |>
  summarise(
    n_tested = sum(!is.na(p)),
    n_up = sum(fdr < 0.05 & logFC > 0, na.rm = TRUE),
    n_down = sum(fdr < 0.05 & logFC < 0, na.rm = TRUE),
    min_fdr = signif(min(fdr, na.rm = TRUE), 3), .by = contrast
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
# representative protein chosen in 00_Gene_Sets.
gene_map <- filter(protein_map, selected)
set_rows <- map(sets, \(genes) match(gene_map$uniprot_id[match(genes, gene_map$gene)], protein_ids))
stopifnot(!any(map_lgl(set_rows, anyNA)))

abundance <- as.matrix(imputed$data)
correlation <- duplicateCorrelation(abundance, design, block = subject)$consensus.correlation
message("within-subject correlation on the imputed matrix: ", round(correlation, 3))

set.seed(1)
fgsea_raw <- map(ranked, \(stats) fgsea::fgsea(sets, stats, minSize = 15, maxSize = 500))

# Both packages return their own BH column over the list they were given, fgsea's padj and fry's
# FDR, kept here as fdr. Neither is recomputed.
set_tests <- bind_rows(
  map(contrast_names, \(contrast) {
    as_tibble(fgsea_raw[[contrast]]) |>
      transmute(
        set_id = pathway, contrast, method = "fgsea", n_proteins = size,
        direction = if_else(NES > 0, "Up", "Down"), nes = NES, p = pval, fdr = padj,
        leading_edge = leadingEdge
      )
  }),
  map(contrast_names, \(contrast) {
    fry(abundance, set_rows, design, fit$design$contrast_matrix[, contrast],
      block = subject, correlation = correlation, sort = "none"
    ) |>
      rownames_to_column("set_id") |>
      transmute(
        set_id, contrast,
        method = "fry", n_proteins = NGenes,
        direction = Direction, nes = NA_real_, p = PValue, fdr = FDR,
        leading_edge = list(NULL)
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
    select(gs$set_catalog, set_id, collection, pathway, source_size, description),
    by = "set_id"
  ) |>
  relocate(contrast, method, set_id, collection, pathway)

# One row per contrast: how many sets each test called, and how many survived collapse. The floor
# is printed first.
set_summary <- set_tests |>
  summarise(
    n_tested_fry = sum(method == "fry"),
    n_tested_fgsea = sum(method == "fgsea"),
    n_fdr_fgsea = sum(method == "fgsea" & fdr < 0.05),
    n_collapsed = sum(method == "fgsea" & fdr < 0.05 & main),
    n_fdr_fry = sum(method == "fry" & fdr < 0.05),
    .by = contrast
  ) |>
  arrange(contrast != d$floor)
print(as.data.frame(set_summary))


# ---- figures -------------------------------------------------------------------------------

# One dot plot per collection and contrast, plus one over all collections. Each shows the ten
# strongest collapse survivors by FDR.
shown <- c(
  "Training_Interaction", "Acute_Interaction", "Training_HR", "Training_LR", "Acute_HR",
  "Acute_LR"
)
survivors <- set_tests |>
  filter(method == "fgsea", fdr < 0.05, main, contrast %in% shown)

draw_dotplot <- function(rows, colour_by) {
  top <- rows |>
    slice_min(fdr, n = 10, with_ties = FALSE) |>
    mutate(label = enrichVolcano::ev_clean_label(pathway))
  figure <- top |>
    ggplot(aes(nes, reorder(label, nes), size = n_proteins, colour = .data[[colour_by]])) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2, 6), name = "proteins") +
    labs(
      x = "normalised enrichment score", y = NULL, title = unique(top$contrast),
      subtitle = sprintf(
        "%s, fgsea, BH within contrast",
        if (colour_by == "collection") "all collections" else rows$collection[1]
      ),
      caption = sprintf(
        "%s of %d collapse survivor%s, ranked by FDR. Size is protein count, colour is %s. %s",
        if (nrow(rows) > 10) "Ten strongest" else "All", nrow(rows),
        if (nrow(rows) == 1) "" else "s",
        if (colour_by == "collection") "collection" else "-log10 FDR",
        "Table: c_data/02_contrasts.xlsx, significant."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.major.y = element_blank(),
      plot.caption = element_text(size = 6.5, colour = "grey45", hjust = 0)
    )
  if (colour_by == "collection") {
    figure + scale_colour_brewer(palette = "Dark2", name = NULL)
  } else {
    figure + scale_colour_viridis_c(
      option = "rocket", direction = -1, end = 0.9,
      name = expression(-log[10] ~ FDR)
    )
  }
}

collections <- unique(gs$set_catalog$collection[gs$set_catalog$qualifies])
dotplots <- map(c("all", collections), \(which_collection) {
  rows <- if (which_collection == "all") {
    survivors
  } else {
    filter(survivors, collection == which_collection)
  }
  map(intersect(shown, unique(rows$contrast)), \(cn) {
    draw_dotplot(
      mutate(filter(rows, contrast == cn), `-log10 FDR` = -log10(fdr)),
      if (which_collection == "all") "collection" else "-log10 FDR"
    )
  })
})

collapse_effect <- set_tests |>
  filter(method == "fgsea", fdr < 0.05, contrast %in% shown) |>
  summarise(before = n(), after = sum(main), .by = c(contrast, collection)) |>
  tidyr::pivot_longer(c(before, after), names_to = "stage", values_to = "n_sets") |>
  mutate(
    stage = factor(stage, c("before", "after")),
    contrast = factor(contrast, levels = shown)
  )
collapse_figure <- ggplot(collapse_effect, aes(stage, n_sets, fill = collection)) +
  geom_col(position = "dodge") +
  facet_wrap(~contrast, scales = "free_y", nrow = 2) +
  scale_fill_brewer(palette = "Dark2", name = NULL) +
  labs(
    x = NULL, y = "significant sets", title = "Significant sets before and after collapse",
    subtitle = "fgsea at FDR 0.05, then collapsePathways",
    caption = paste(
      "collapsePathways re-tests each significant set conditioned on a stronger set's leading",
      "edge and keeps it only if it stays significant.",
      "Table: c_data/02_contrasts.xlsx, significant."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.caption = element_text(size = 7, colour = "grey45", hjust = 0))

# Every set nominal under fry in at least one contrast, as a dot matrix. Fill is the fgsea NES,
# for direction; size is -log10 fry p; a black ring marks fry FDR < 0.05.
# One label per set, built before the join so two sets sharing a display name stay apart.
set_labels <- gs$set_catalog |>
  filter(qualifies) |>
  transmute(
    set_id,
    label = gsub("\n", " ", paste0("[", collection, "] ", enrichVolcano::ev_clean_label(pathway)))
  ) |>
  mutate(label = if_else(duplicated(label) | duplicated(label, fromLast = TRUE), set_id, label))
hit_data <- set_tests |>
  filter(method == "fry") |>
  select(set_id, contrast, p, fdr) |>
  left_join(
    set_tests |> filter(method == "fgsea") |> select(set_id, contrast, effect = nes),
    by = c("set_id", "contrast")
  ) |>
  left_join(set_labels, by = "set_id")
ranked_hits <- hit_data |>
  summarise(n_nominal = sum(p < 0.05), best = min(p), .by = label) |>
  filter(n_nominal > 0) |>
  arrange(desc(n_nominal), best)
page_of <- ceiling(seq_len(nrow(ranked_hits)) / 75)
hit_plots <- map(seq_len(max(page_of, 0)), \(page) {
  rows <- ranked_hits$label[page_of == page]
  page_data <- hit_data |>
    filter(label %in% rows) |>
    mutate(
      label = factor(label, levels = rev(rows)),
      contrast = factor(contrast, levels = contrast_names), nominal = p < 0.05
    )
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
    scale_x_discrete(drop = FALSE) +
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
        "FDR < 0.05. Grey speck: tested, not nominal. Table: set_tests in c_data/set_tests.rds."
      )
    ) +
    theme_minimal(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      axis.text.y = element_text(size = if (length(rows) > 30) 5.5 else 8),
      plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)
    )
})

pages <- c(dotplots[[1]], list(collapse_figure), list_flatten(dotplots[-1]), hit_plots)
pdf(file.path(figure_dir, "02_contrasts_figures.pdf"), width = 11, height = 8.5)
walk(pages, print)
invisible(dev.off())

saveRDS(
  list(protein_results = protein_results, set_tests = set_tests),
  file.path(out, "set_tests.rds"),
  compress = "xz"
)
sheets <- list(
  set_summary = set_summary,
  significant = set_tests |>
    filter(fdr < 0.05) |>
    mutate(leading_edge = map_chr(leading_edge, paste, collapse = ";")),
  protein_summary = protein_summary,
  fry_correlation = tibble(matrix = "imputed", correlation = correlation),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Sets tested and called per contrast, by method, and fgsea survivors of collapsePathways.",
  "Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds.",
  "Proteins tested and at BH < 0.05 per contrast.",
  "Within-subject correlation fry used, estimated on the imputed matrix.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
writexl::write_xlsx(
  c(list(read_me = read_me), sheets), file.path(out, "02_contrasts.xlsx")
)
message("wrote set_tests.rds, 02_contrasts.xlsx and ", length(pages), " figure pages")
