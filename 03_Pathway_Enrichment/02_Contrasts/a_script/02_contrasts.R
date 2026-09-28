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


# Every fgsea collapse survivor at FDR < 0.05 and every set at fry p < 0.05, one dot plot per
# contrast and method, 40 sets to a page.
labelled <- set_tests |>
  mutate(
    contrast = factor(contrast, levels = contrast_names),
    label = gsub("\n", " ", enrichVolcano::ev_clean_label(pathway))
  )
dot_pages <- function(rows, x, x_label, method_label, caption) {
  rows |>
    arrange(contrast, .data[[x]]) |>
    mutate(page = (row_number() - 1) %/% 40 + 1, n_pages = max(page), .by = contrast) |>
    group_by(contrast, page) |>
    group_map(\(page_rows, key) {
      ggplot(page_rows, aes(.data[[x]], reorder(label, .data[[x]]), colour = collection)) +
        geom_vline(xintercept = 0, colour = "grey75") +
        geom_point(aes(size = n_proteins)) +
        scale_colour_brewer(palette = "Dark2", name = NULL) +
        scale_size_continuous(range = c(1.5, 5), name = "proteins") +
        labs(
          title = sprintf("%s, %s", key$contrast, method_label),
          subtitle = sprintf(
            "%d sets. Page %d of %d.",
            sum(rows$contrast == key$contrast), key$page, page_rows$n_pages[1]
          ),
          x = x_label, y = NULL, caption = caption
        ) +
        theme_minimal(base_size = 8) +
        theme(plot.caption = element_text(hjust = 0), plot.title = element_text(face = "bold"))
    })
}
fgsea_pages <- dot_pages(
  filter(labelled, method == "fgsea", fdr < 0.05, main), "nes",
  "normalised enrichment score", "fgsea collapse survivors at FDR < 0.05",
  paste(
    "fgsea on the moderated t ranking; BH within contrast over all collections; collapsePathways",
    "keeps the sets that stay significant given a stronger set. Table: c_data/02_contrasts.xlsx,",
    "significant."
  )
)
fry_pages <- dot_pages(
  labelled |>
    filter(method == "fry", p < 0.05) |>
    mutate(signed = if_else(direction == "Up", 1, -1) * -log10(p)),
  "signed", "signed -log10 p (positive: up)", "fry sets at p < 0.05",
  paste(
    "fry on the imputed matrix with the fitted design, subject block and within-subject",
    "correlation; BH within contrast. Table: set_tests in c_data/set_tests.rds."
  )
)
pages <- c(fgsea_pages, fry_pages)
stopifnot(
  sum(map_int(pages, \(g) nrow(g$data))) ==
    nrow(filter(set_tests, method == "fgsea", fdr < 0.05, main)) +
      nrow(filter(set_tests, method == "fry", p < 0.05))
)

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
