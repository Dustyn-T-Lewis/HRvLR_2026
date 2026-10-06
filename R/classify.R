# Single-feature discrimination between two groups, shared by the protein, set and module levels.
# Nothing is fitted: each feature's AUC is the Mann-Whitney statistic W / (n1 * n2), above 0.5 when
# the feature is higher in the first group, with the Wilcoxon rank-sum p, BH and Storey q within
# task (and collection when collection_of is given). Five tasks compare the groups: the T1, T2 and
# T3 levels, and the training (T2 - T1) and acute (T3 - T2) changes. Every AUC carries DeLong's
# 95% interval, which runs narrow at this n; the features drawn in the report also carry pROC's
# stratified bootstrap interval. The share of features at p < 0.05 is a multiple-testing summary,
# not a measure of how well the groups separate.
classify_step <- function(stage, inputs, paths, values, meta, feature_label, feature_noun,
                          feature_unit, value_label, row_groups, row_order,
                          collection_of = NULL, panel_name = \(label, collection) label,
                          observed_share = NULL, groups = c("HR", "LR"), boot_n = 2000) {
  out <- file.path(stage, "c_data")
  figures <- file.path(stage, "b_reports")
  step <- tolower(basename(stage))
  workbook <- paste0(step, ".xlsx")
  for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  by_collection <- if (is.null(collection_of)) character(0) else "collection"
  labelled_heat <- \(...) heat(..., labels = feature_label)
  bh_scope <- if (length(by_collection)) "task and collection," else "task,"
  labelled <- function(d) {
    if (length(by_collection)) d <- mutate(d, collection = collection_of[feature], .after = feature)
    d <- mutate(d, label = feature_label[feature], .after = all_of(c("feature", by_collection)))
    if (is.null(observed_share)) {
      return(d)
    }
    mutate(d, observed_share = observed_share[feature], .after = label)
  }

  arm_of <- set_names(as.character(meta$arm), meta$subject)
  by_subject <- function(timepoint) {
    rows <- filter(meta, timepoint == .env$timepoint, arm %in% groups)
    set_names(as.data.frame(values[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
  }
  change <- function(from, to) {
    before <- by_subject(from)
    after <- by_subject(to)
    both <- intersect(colnames(before), colnames(after))
    after[, both, drop = FALSE] - before[, both, drop = FALSE]
  }
  split_groups <- function(level, label) {
    first <- arm_of[colnames(level)] == groups[1]
    list(
      positive = level[, first, drop = FALSE], negative = level[, !first, drop = FALSE],
      label = label
    )
  }
  versus <- paste(groups, collapse = "v")
  label_versus <- paste(groups, collapse = " vs ")
  tasks <- set_names(list(
    split_groups(by_subject("T1"), paste(label_versus, "at T1")),
    split_groups(by_subject("T2"), paste(label_versus, "at T2")),
    split_groups(by_subject("T3"), paste(label_versus, "at T3")),
    split_groups(change("T1", "T2"), paste0(label_versus, ", training change")),
    split_groups(change("T2", "T3"), paste0(label_versus, ", acute change"))
  ), paste0(c("Baseline_", "Trained_", "Acute_", "Training_change_", "Acute_change_"), versus))

  auc <- imap(tasks, \(spec, name) {
    fits <- vapply(rownames(values), \(id) classify_row(spec, id), numeric(6))
    tibble(task = name, task_label = spec$label, feature = rownames(values), !!!as_tibble(t(fits)))
  }) |>
    list_rbind() |>
    labelled() |>
    mutate(fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = c(task, all_of(by_collection)))

  # The smallest p the rank-sum test can reach with n1 and n2 subjects: 2 / choose(n1 + n2, n1).
  counts <- auc |>
    filter(!is.na(p)) |>
    summarise(
      n_positive = max(n_positive), n_negative = max(n_negative),
      n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
      .by = c(task, task_label, all_of(by_collection))
    ) |>
    mutate(min_p = 2 / choose(n_positive + n_negative, n_positive), .after = n_negative)

  task_labels <- map_chr(tasks, "label")
  overviews <- imap(row_groups, \(rows, group) {
    task_grid <- \(value) as_grid(auc, rows, names(tasks), "task", value)
    maps <- labelled_heat(
      task_grid(\(d) d$auc), sprintf("AUC, %s", label_versus), c(0, 0.5, 1), "AUC",
      column_labels = task_labels
    ) +
      labelled_heat(task_grid(\(d) status_of(d$p, d$fdr)), "Wilcoxon p and BH q",
        row_names = TRUE, column_labels = task_labels
      )
    heatmap_page(maps, sprintf("Discrimination overview: %s", group), paste(
      sprintf("All %s %s, one row each.", format(length(rows), big.mark = ","), feature_noun),
      row_order, paste0("Grey: not tested. Table: c_data/", workbook, ", auc.")
    ))
  })

  # The report's hits: the strongest p < 0.05 features per task, each with a bootstrap interval.
  hits <- auc |>
    filter(p < 0.05) |>
    mutate(task = factor(task, names(tasks))) |>
    arrange(task, p, desc(abs(auc - 0.5)), feature) |>
    paginate("task")
  set.seed(1)
  boot <- pmap(list(as.character(hits$task), hits$feature), \(task, id) {
    s <- sides(tasks[[task]], id)
    pROC::roc(controls = s$neg, cases = s$pos, direction = "<", quiet = TRUE) |>
      pROC::ci.auc(
        method = "bootstrap", boot.n = boot_n, boot.stratified = TRUE, progress = "none"
      ) |>
      as.numeric()
  })
  hits <- hits |>
    mutate(
      boot_lower = map_dbl(boot, 1), boot_upper = map_dbl(boot, 3),
      name = panel_name(label, pick(everything())[["collection"]]),
      panel = sprintf(
        "%s\nAUC %.2f [%.2f, %.2f]  p %.2g  q %.2g\n%s %d, %s %d", name, auc, boot_lower,
        boot_upper, p, fdr, groups[1], n_positive, groups[2], n_negative
      )
    )
  auc <- left_join(
    auc, select(hits, task, feature, boot_lower, boot_upper) |> mutate(task = as.character(task)),
    by = c("task", "feature")
  )
  task_page <- function(page_hits) {
    spec <- tasks[[as.character(page_hits$task[1])]]
    info <- filter(counts, task == as.character(page_hits$task[1])) |>
      summarise(across(c(task_label, n_positive, n_negative, min_p), first),
        across(c(n_tested, n_nominal, n_fdr_05), sum),
        .by = task
      )
    curves <- map(page_hits$feature, \(id) {
      s <- sides(spec, id)
      pROC::roc(controls = s$neg, cases = s$pos, direction = "<", quiet = TRUE)
    }) |>
      set_names(page_hits$panel)
    pROC::ggroc(curves, aes = "group", legacy.axes = TRUE, colour = "black") +
      geom_abline(
        aes(intercept = 0, slope = 1),
        data = tibble(name = page_hits$panel), linetype = 2, colour = "grey60"
      ) +
      facet_wrap(~ factor(name, levels = slots(page_hits$panel)),
        nrow = 3, ncol = 4, drop = FALSE
      ) +
      coord_equal() +
      labs(
        title = info$task_label, subtitle = counts_line(info, feature_noun, page_hits),
        x = "1 - specificity", y = "sensitivity", caption = sprintf(paste(
          "%s (up to %d) against %s (up to %d); AUC above 0.5 means higher in %s. Header: AUC with",
          "its stratified bootstrap 95%% interval (%d resamples), Wilcoxon rank-sum p (smallest",
          "attainable %.2g), BH q within", bh_scope, "and the feature's n.",
          paste0("Table: c_data/", workbook, ", auc.")
        ), groups[1], info$n_positive, groups[2], info$n_negative, groups[1], boot_n, info$min_p)
      ) +
      page_theme
  }
  pages <- map(split(hits, hits$pdf_page), task_page)

  plot_index <- bind_rows(
    tibble(task = "overview", label = map_chr(overviews, "title"), pdf_page = seq_along(overviews)),
    transmute(hits,
      task = as.character(task), pick(all_of(by_collection)), feature, label,
      pdf_page = pdf_page + length(overviews), position
    )
  )
  stopifnot(
    nrow(hits) == capped_count(filter(auc, p < 0.05), "task"),
    !anyDuplicated(select(filter(plot_index, !is.na(feature)), -label, -pdf_page, -position)),
    max(plot_index$pdf_page) == length(overviews) + length(pages)
  )

  pdf(file.path(figures, paste0(step, "_figures.pdf")), width = 11, height = 8.5)
  walk(overviews, \(page) page$draw())
  walk(pages, print)
  invisible(dev.off())

  write_workbook(
    list(
      counts = counts, auc = arrange(auc, task, pick(all_of(by_collection)), p),
      plot_index = plot_index
    ),
    c(
      "Per task: n per group, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts.",
      sprintf(paste(
        "Per %s and task: AUC with DeLong 95%% interval, Wilcoxon p, BH and Storey q; bootstrap",
        "95%% interval for the features drawn in the report."
      ), feature_unit),
      "PDF page of each overview, and page and panel position (1 to 12) of every drawn row."
    ),
    file.path(out, workbook),
    inputs, paths
  )
}
