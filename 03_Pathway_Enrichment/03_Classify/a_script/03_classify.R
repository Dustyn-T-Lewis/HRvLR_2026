# How well each gene set's singscore separates the study's groups, on the protein level's eight
# tasks. Each collection is its own screen, so BH runs within task and collection.
# Within-arm tasks pair each subject's two biopsies: p is the Wilcoxon signed-rank test and the
# effect is the median paired change. Between-arm tasks compare HR with LR subjects on a level or
# on their own change: p is the rank-sum test and AUC is W / (n1 * n2), above 0.5 when higher in
# HR. The auc column is filled for every task; for paired tasks it ignores the pairing and is not
# drawn. Between-arm AUCs carry a DeLong 95% interval.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(grid)
  library(ComplexHeatmap)
  library(writexl)
})
source(here("R", "helpers.R"))

stage <- here("03_Pathway_Enrichment", "03_Classify")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  singscore = "03_Pathway_Enrichment/01_Scores/c_data/singscore.rds",
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  proteins = "01_Preprocess/03_Quantification/c_data/DAList_corrected.rds"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
catalog <- readRDS(paths[["gene_sets"]])$set_catalog |> filter(qualifies)
values <- readRDS(paths[["singscore"]])$scores[catalog$set_id, ]
meta <- readRDS(paths[["proteins"]])$metadata
stopifnot(identical(colnames(values), meta$sample_id))
feature_label <- with(catalog, set_names(enrichVolcano::ev_clean_label(pathway), set_id))
collection_of <- with(catalog, set_names(collection, set_id))
# ev_clean_label() breaks names every few words; panel headers re-wrap them at 50 characters.
panel_label <- \(x) stringr::str_wrap(gsub("\n", " ", x), 50)
feature_noun <- "sets"
value_label <- "singscore"
# Display order within each collection: average-linkage clustering on 1 - Pearson r of singscore.
row_groups <- split(catalog$set_id, factor(catalog$collection, unique(catalog$collection))) |>
  map(\(ids) ids[hclust(as.dist(1 - cor(t(values[ids, ]))), "average")$order])

by_subject <- function(timepoint, arm = c("HR", "LR")) {
  rows <- filter(meta, timepoint == .env$timepoint, arm %in% .env$arm)
  set_names(as.data.frame(values[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
}
paired <- function(from, to, arm = c("HR", "LR")) {
  a <- by_subject(from, arm)
  b <- by_subject(to, arm)
  both <- intersect(colnames(a), colnames(b))
  list(before = a[, both, drop = FALSE], after = b[, both, drop = FALSE])
}
arm_of <- set_names(as.character(meta$arm), meta$subject)
within_arm <- function(from, to, arm, label) {
  p <- paired(from, to, arm)
  list(
    positive = p$after, negative = p$before, paired = TRUE, label = label,
    levels = c(from, to)
  )
}
between_arm <- function(level, label) {
  hr <- arm_of[colnames(level)] == "HR"
  list(
    positive = level[, hr, drop = FALSE], negative = level[, !hr, drop = FALSE],
    paired = FALSE, label = label
  )
}
change <- \(from, to) with(paired(from, to), after - before)
tasks <- list(
  Training_HR = within_arm("T1", "T2", "HR", "Training, HR (T1 to T2)"),
  Training_LR = within_arm("T1", "T2", "LR", "Training, LR (T1 to T2)"),
  Acute_HR = within_arm("T2", "T3", "HR", "Acute bout, HR (T2 to T3)"),
  Acute_LR = within_arm("T2", "T3", "LR", "Acute bout, LR (T2 to T3)"),
  Baseline_HRvLR = between_arm(by_subject("T1"), "HR vs LR at T1"),
  Trained_HRvLR = between_arm(by_subject("T2"), "HR vs LR at T2"),
  Training_change_HRvLR = between_arm(change("T1", "T2"), "HR vs LR, training change"),
  Acute_change_HRvLR = between_arm(change("T2", "T3"), "HR vs LR, acute change")
)

auc <- imap(tasks, \(spec, name) {
  fits <- vapply(rownames(values), \(id) classify_row(spec, id), numeric(7))
  tibble(
    task = name, task_label = spec$label, paired = spec$paired, feature = rownames(values),
    !!!as_tibble(t(fits))
  )
}) |>
  list_rbind() |>
  mutate(collection = collection_of[feature], label = feature_label[feature], .after = feature) |>
  mutate(fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = c(task, collection))

# The smallest p a task can reach: signed-rank on n pairs gives 2 / 2^n, rank-sum on n1 and n2
# gives 2 / choose(n1 + n2, n1). HR has 6 training and 7 acute pairs, LR 8 of each.
counts <- auc |>
  filter(!is.na(p)) |>
  summarise(
    paired = first(paired), n_positive = max(n_positive), n_negative = max(n_negative),
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(task, task_label, collection)
  ) |>
  mutate(
    min_p = if_else(paired, 2 / 2^n_positive, 2 / choose(n_positive + n_negative, n_positive)),
    .after = n_negative
  )



# Overview: every feature by task, median paired change and AUC each on its own scale.
paired_tasks <- names(tasks)[map_lgl(tasks, "paired")]
task_labels <- map_chr(tasks, "label")
change_limit <- cap(auc$median_change)
overviews <- imap(row_groups, \(rows, group) {
  task_grid <- \(cols, value) as_grid(auc, rows, cols, "task", value)
  between <- setdiff(names(tasks), paired_tasks)
  maps <- heat(task_grid(paired_tasks, \(d) d$median_change), "Median paired change",
    c(-change_limit, 0, change_limit), paste("median change,", value_label),
    column_labels = task_labels[paired_tasks]
  ) +
    heat(task_grid(between, \(d) d$auc), "AUC, HR against LR", c(0, 0.5, 1), "AUC",
      column_labels = task_labels[between]
    ) +
    heat(task_grid(names(tasks), \(d) status_of(d$p, d$fdr)), "Wilcoxon p and BH q",
      row_names = TRUE, column_labels = task_labels
    )
  heatmap_page(maps, sprintf("Classification overview: %s", group), paste(
    sprintf("All %s %s, one row each.", format(length(rows), big.mark = ","), feature_noun),
    "Rows are ordered within collection by average-linkage clustering (1 - Pearson r)",
    "of singscore, for display only.",
    sprintf("Median change colours stop at the 99th percentile of |change|, %.2g.", change_limit),
    "Grey: not tested. Table: c_data/03_classify.xlsx, auc."
  ))
})

# Grids: every p < 0.05 hit, strongest first (p, then the larger effect, then feature ID).
# Within-arm tasks join each subject's two biopsies; HR-vs-LR tasks draw an ROC curve.
task_counts <- counts |>
  summarise(
    across(c(task_label, paired, n_positive, n_negative, min_p), first),
    across(c(n_tested, n_nominal, n_fdr_05), sum),
    .by = task
  )
hits <- auc |>
  filter(p < 0.05) |>
  mutate(
    task = factor(task, names(tasks)),
    name = paste0(collection, ": ", panel_label(label)),
    panel = if_else(paired,
      sprintf(
        "%s\nmedian change %+.2f  p %.2g  q %.2g\nn pairs %d", name, median_change, p, fdr,
        n_positive
      ),
      sprintf(
        "%s\nAUC %.2f %s  p %.2g  q %.2g\nHR %d, LR %d", name, auc,
        if_else(is.na(auc_lower), "[no CI]", sprintf("[%.2f, %.2f]", auc_lower, auc_upper)),
        p, fdr, n_positive, n_negative
      )
    )
  ) |>
  arrange(task, p, desc(abs(if_else(paired, median_change, auc - 0.5))), feature) |>
  paginate("task")
task_page <- function(hits) {
  spec <- tasks[[as.character(hits$task[1])]]
  info <- filter(task_counts, task == as.character(hits$task[1]))
  panels <- slots(hits$panel)
  if (spec$paired) {
    points <- map2(hits$feature, hits$panel, \(id, panel) {
      tibble(
        panel,
        subject = colnames(spec$positive),
        before = spec$negative[id, ], after = spec$positive[id, ]
      )
    }) |>
      list_rbind() |>
      filter(!is.na(before), !is.na(after)) |>
      pivot_longer(c(before, after), names_to = "biopsy") |>
      mutate(
        biopsy = factor(spec$levels[match(biopsy, c("before", "after"))], spec$levels),
        panel = factor(panel, panels)
      )
    plot <- ggplot(points, aes(biopsy, value, group = subject)) +
      geom_line(colour = "grey55") +
      geom_point(size = 1.4) +
      facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free_y", drop = FALSE) +
      labs(x = NULL, y = value_label, caption = sprintf(paste(
        "Wilcoxon signed-rank test on up to %d pairs; smallest attainable p %.2g. Header: median",
        "paired change, p, BH q within task and collection,",
        "and the feature's complete pairs.",
        "Table: c_data/03_classify.xlsx, auc."
      ), info$n_positive, info$min_p))
  } else {
    curves <- map(hits$feature, \(id) {
      s <- sides(spec, id)
      pROC::roc(controls = s$neg, cases = s$pos, direction = "<", quiet = TRUE)
    }) |>
      set_names(hits$panel)
    plot <- pROC::ggroc(curves, aes = "group", legacy.axes = TRUE, colour = "black") +
      geom_abline(
        aes(intercept = 0, slope = 1),
        data = tibble(name = hits$panel), linetype = 2, colour = "grey60"
      ) +
      facet_wrap(~ factor(name, levels = panels), nrow = 3, ncol = 4, drop = FALSE) +
      coord_equal() +
      labs(x = "1 - specificity", y = "sensitivity", caption = sprintf(paste(
        "HR (up to %d) against LR (up to %d); AUC above 0.5 means higher in HR. Wilcoxon rank-sum",
        "test; smallest attainable p %.2g. Header: AUC with its DeLong 95%% interval, p,",
        "none at AUC 0 or 1, where DeLong's variance is zero; at this n it runs narrow.",
        "BH q within task and collection,",
        "and the feature's n.",
        "Table: c_data/03_classify.xlsx, auc."
      ), info$n_positive, info$n_negative, info$min_p))
  }
  plot +
    labs(title = info$task_label, subtitle = counts_line(info, feature_noun, hits)) +
    page_theme
}
pages <- map(split(hits, hits$pdf_page), task_page)

plot_index <- bind_rows(
  tibble(task = "overview", label = map_chr(overviews, "title"), pdf_page = seq_along(overviews)),
  transmute(hits,
    task = as.character(task), collection, feature, label, pdf_page = pdf_page + length(overviews),
    position
  )
)
stopifnot(
  nrow(hits) == sum(auc$p < 0.05, na.rm = TRUE),
  !anyDuplicated(select(filter(plot_index, !is.na(feature)), -label, -pdf_page, -position)),
  max(plot_index$pdf_page) == length(overviews) + length(pages)
)

pdf(file.path(figures, "03_classify_figures.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(pages, print)
invisible(dev.off())

write_workbook(
  list(counts = counts, auc = arrange(auc, task, collection, p), plot_index = plot_index),
  c(
    "Per task: pairing, n per side, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts.",
    "Per set and task: AUC, median paired change (paired tasks), Wilcoxon p, BH and Storey q.",
    "PDF page of each overview, and page and panel position (1 to 12) of every p < 0.05 row."
  ),
  file.path(out, "03_classify.xlsx"),
  inputs, paths
)
