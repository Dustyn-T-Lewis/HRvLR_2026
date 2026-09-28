# fgsea and fry, side by side. fgsea is competitive on moderated t and
# assumes exchangeable proteins, which they are not. fry is self-contained and rotates residuals
# under the design, subject block and within-subject correlation. fry takes no missing value, so
# it reads the imputed matrix with a correlation estimated on that matrix. Each method's own BH
# runs within contrast over the sets it tested: all of them for fry, fewer for fgsea when a
# contrast leaves a set under minSize.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
  library(writexl)
})

stage <- here("03_Pathway_Enrichment", "02_Contrasts")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  fit = "02_Differential_Expression/02_Contrasts/c_data/fit.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  singscore = "03_Pathway_Enrichment/01_Scores/c_data/singscore.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
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

catalog <- filter(gs$set_catalog, qualifies)
sample_values <- readRDS(paths[["singscore"]])$scores[catalog$set_id, ]
column_meta <- imputed$metadata
feature_label <- with(catalog, set_names(enrichVolcano::ev_clean_label(pathway), set_id))
feature_noun <- "sets"
contrast_results <- set_tests |>
  filter(method == "fgsea") |>
  transmute(feature = set_id, contrast, effect = nes, p, fdr)
# Display order within each collection: average-linkage clustering on 1 - Pearson r of singscore.
row_groups <- split(catalog$set_id, factor(catalog$collection, unique(catalog$collection))) |>
  map(\(ids) ids[hclust(as.dist(1 - cor(t(sample_values[ids, ]))), "average")$order])
# Heatmap and workbook helpers, the same in every contrasts, classify and associate script.
# Display only: row order, colour limits and centring touch no test.
ht_opt$message <- FALSE
diverging <- c("#5E3C99", "white", "#E66101")
status_colours <- c("not significant" = "white", "p < 0.05" = "#35978F", "BH < 0.05" = "black")
status_legend <- Legend(
  labels = c(names(status_colours), "not tested"), title = "status",
  legend_gp = gpar(fill = c(status_colours, "grey80"))
)
status_of <- \(p, fdr) {
  case_when(fdr < 0.05 ~ "BH < 0.05", p < 0.05 ~ "p < 0.05", !is.na(p) ~ "not significant")
}
# Colours stop at the 99th percentile of |x| so a few extremes do not wash out the rest.
cap <- \(x) unname(quantile(abs(x), 0.99, na.rm = TRUE))
# Long results to a feature-by-column matrix; untested cells stay NA and draw grey.
as_grid <- function(d, rows, cols, col, value) {
  m <- matrix(NA, length(rows), length(cols), dimnames = list(rows, cols))
  d <- filter(d, feature %in% rows, .data[[col]] %in% cols)
  m[cbind(d$feature, d[[col]])] <- value(d)
  m
}
# An effect map on a diverging scale when breaks are given, otherwise a status map.
heat <- function(m, title, breaks = NULL, legend = title, row_names = FALSE,
                 column_title = title, column_names_gp = gpar(fontsize = 7), ...) {
  rows <- rownames(m)
  Heatmap(m,
    name = title, na_col = "grey80",
    col = if (is.null(breaks)) status_colours else circlize::colorRamp2(breaks, diverging),
    show_heatmap_legend = !is.null(breaks) && !is.null(legend),
    heatmap_legend_param = list(title = legend, at = if (!is.null(breaks)) signif(breaks, 2)),
    cluster_rows = FALSE, cluster_columns = FALSE, column_title = column_title,
    column_title_gp = gpar(fontsize = 8), column_names_gp = column_names_gp,
    show_row_names = row_names && length(rows) <= 80,
    row_labels = gsub("\n", " ", feature_label[rows]),
    row_names_gp = gpar(fontsize = if (length(rows) <= 40) 7 else 5),
    use_raster = length(rows) > 200, ...
  )
}
heatmap_page <- function(maps, title, caption, legends = list(status_legend)) {
  list(title = title, draw = function() {
    draw(maps,
      column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
      merge_legend = TRUE, heatmap_legend_list = legends, padding = unit(c(14, 4, 4, 4), "mm")
    )
    grid.text(stringr::str_wrap(caption, 190),
      x = unit(4, "mm"), y = unit(3, "mm"), just = c("left", "bottom"), gp = gpar(fontsize = 7)
    )
  })
}
write_workbook <- function(sheets, holds, file) {
  sheets <- c(sheets, list(
    input_manifest = tibble(
      input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
    ),
    package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
      as_tibble() |>
      select(package, version = loadedversion, source)
  ))
  holds <- c(holds, "Files read, with md5.", "Packages loaded at run time.")
  write_xlsx(c(list(read_me = tibble(sheet = names(sheets), holds = holds)), sheets), file)
}

# Overviews: every feature across the 45 biopsies, then across the nine contrasts. Biopsy traits
# sit on the T1 and T2 columns; T3 has none. comp_hypertrophy is one value per subject.
arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
trait_table <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
biopsy_traits <- trait_table |>
  select(subject, matches("_(pre|post)_")) |>
  pivot_longer(-subject, names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*") |>
  pivot_wider(names_from = trait) |>
  transmute(
    subject,
    timepoint = if_else(when == "pre", "T1", "T2"),
    type1_share = 100 * fibres_type1 / fibres_mixed, fcsa_mixed, mcsa, leg_press_1rm
  )
columns <- as_tibble(column_meta) |>
  transmute(sample_id, subject, arm = as.character(arm), timepoint = as.character(timepoint)) |>
  mutate(group = factor(paste(arm, timepoint, sep = "_"), paste(
    rep(c("HR", "LR"), each = 3), c("T1", "T2", "T3"),
    sep = "_"
  ))) |>
  left_join(biopsy_traits, by = c("subject", "timepoint")) |>
  left_join(select(trait_table, subject, comp_hypertrophy), by = "subject") |>
  arrange(group, subject)
stopifnot(!anyDuplicated(columns$sample_id), setequal(columns$sample_id, colnames(sample_values)))
grey_ramp <- \(x) circlize::colorRamp2(range(x, na.rm = TRUE), c("#F7F7F7", "#252525"))
tracks <- list(
  "type I share (%)" = columns$type1_share, "fCSA mixed (um2)" = columns$fcsa_mixed,
  "mCSA (cm2)" = columns$mcsa, "leg press 1RM (kg)" = columns$leg_press_1rm,
  "comp_hypertrophy, per subject" = columns$comp_hypertrophy
)
sample_annotation <- HeatmapAnnotation(
  df = data.frame(arm = columns$arm, timepoint = columns$timepoint, tracks, check.names = FALSE),
  col = c(
    list(arm = arm_colours, timepoint = c(T1 = "grey85", T2 = "grey55", T3 = "grey25")),
    map(tracks, grey_ramp)
  ),
  na_col = "grey80", annotation_name_gp = gpar(fontsize = 7), annotation_name_side = "left",
  simple_anno_size = unit(3, "mm")
)
row_order <- paste(
  "Rows are ordered within collection by average-linkage clustering (1 - Pearson r)",
  "of singscore, for display only."
)
effect_limit <- cap(contrast_results$effect)
overviews <- imap(row_groups, \(rows, group) {
  m <- sample_values[rows, columns$sample_id, drop = FALSE]
  m <- sweep(m, 1, rowMeans(m))
  limit <- cap(m)
  biopsy_map <- heat(m, "biopsies", c(-limit, 0, limit), "singscore, centred",
    row_names = TRUE, column_title = character(0), top_annotation = sample_annotation,
    column_split = columns$group, cluster_column_slices = FALSE,
    column_labels = columns$subject, column_names_gp = gpar(fontsize = 5)
  )
  contrast_grid <- \(value) as_grid(contrast_results, rows, contrast_names, "contrast", value)
  contrast_maps <- heat(contrast_grid(\(d) d$effect), "fgsea NES",
    c(-effect_limit, 0, effect_limit), "NES"
  ) +
    heat(contrast_grid(\(d) status_of(d$p, d$fdr)), "fgsea p and BH q", row_names = TRUE)
  n_rows <- format(length(rows), big.mark = ",")
  list(
    heatmap_page(biopsy_map, sprintf("Biopsy overview: %s", group), paste(
      sprintf("All %s %s by %d biopsies,", n_rows, feature_noun, ncol(m)),
      "grouped by arm and time point; column labels are subjects.",
      "Fill: singscore minus the set's mean over biopsies.",
      sprintf("Colours stop at the 99th percentile of |value|, %.2g.", limit),
      row_order, "Grey: missing."
    ), legends = list()),
    heatmap_page(contrast_maps, sprintf("Contrast overview: %s", group), paste(
      sprintf("All %s %s by the nine contrasts.", n_rows, feature_noun),
      "Fill: fgsea normalised enrichment score; status: p and BH q within contrast,",
      "collections pooled.",
      sprintf("Colours stop at the 99th percentile of |value|, %.2g.", effect_limit),
      row_order, "Grey: not tested. Table: c_data/set_tests.rds, set_tests."
    ))
  )
}) |>
  list_flatten()
stopifnot(length(overviews) == 2 * length(row_groups))

pdf(file.path(figures, "02_contrasts_figures.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(pages, print)
invisible(dev.off())

saveRDS(
  list(protein_results = protein_results, set_tests = set_tests),
  file.path(out, "set_tests.rds"),
  compress = "xz"
)
write_workbook(
  list(
    set_summary = set_summary,
    significant = set_tests |>
      filter(fdr < 0.05) |>
      mutate(leading_edge = map_chr(leading_edge, paste, collapse = ";")),
    protein_summary = protein_summary,
    fry_correlation = tibble(matrix = "imputed", correlation = correlation)
  ),
  c(
    "Sets tested and called per contrast, by method, and fgsea survivors of collapsePathways.",
    "Every set at FDR < 0.05, either method. The full table is set_tests in set_tests.rds.",
    "Proteins tested and at BH < 0.05 per contrast.",
    "Within-subject correlation fry used, estimated on the imputed matrix."
  ),
  file.path(out, "02_contrasts.xlsx")
)
