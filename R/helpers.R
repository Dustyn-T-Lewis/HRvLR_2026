# Helpers every stage script sources after its library() calls: the statistics written for this
# study, and the heatmap, page and workbook plumbing the contrasts, classify and associate steps
# share. tests/testthat/test-helpers.R checks the statistics against hand-computed answers.


# Removes each feature's blood slope. The design carries the group means beside a centred `blood`
# column, so the slope comes from spread within cells and no contrast leaks into it. A feature with
# no estimable slope stays as it was.
remove_blood <- function(values, design, subject) {
  correlation <- duplicateCorrelation(values, design, block = subject)$consensus.correlation
  fit <- lmFit(values, design, block = subject, correlation = correlation) |>
    eBayes()
  slope <- fit$coefficients[, "blood"]
  list(
    values = values - outer(coalesce(slope, 0), design[, "blood"]),
    slope = slope, p = fit$p.value[, "blood"]
  )
}

# ANCOVA per feature and window: the later value on the earlier one plus arm, over subjects with
# both biopsies and at least three per arm. The arm term is HR minus LR at equal starting value.
ancova <- function(values, meta) {
  by_subject <- function(timepoint) {
    rows <- filter(meta, timepoint == .env$timepoint)
    set_names(as.data.frame(values[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
  }
  windows <- list(training = c("T1", "T2"), acute = c("T2", "T3"))
  imap(windows, \(timepoints, window) {
    before <- by_subject(timepoints[1])
    after <- by_subject(timepoints[2])
    both <- intersect(colnames(before), colnames(after))
    arm <- factor(meta$arm[match(both, meta$subject)], c("LR", "HR"))
    map(rownames(values), \(feature) {
      x <- before[feature, both]
      y <- after[feature, both]
      ok <- !is.na(x) & !is.na(y)
      n <- table(arm[ok])
      if (min(n) < 3) {
        return(NULL)
      }
      coefs <- summary(lm(y[ok] ~ x[ok] + arm[ok]))$coefficients
      tibble(
        feature, window,
        n_HR = n[["HR"]], n_LR = n[["LR"]], effect = coefs[3, 1], t = coefs[3, 3], p = coefs[3, 4]
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}

# A feature needs three observations per side. Ties push wilcox.test to its normal approximation
# and it warns on every call.
sides <- function(spec, id) {
  pos <- spec$positive[id, ]
  neg <- spec$negative[id, ]
  if (spec$paired) {
    kept <- !is.na(pos) & !is.na(neg)
    return(list(pos = pos[kept], neg = neg[kept]))
  }
  list(pos = pos[!is.na(pos)], neg = neg[!is.na(neg)])
}
classify_row <- function(spec, id) {
  s <- sides(spec, id)
  n <- c(n_positive = length(s$pos), n_negative = length(s$neg))
  if (min(n) < 3) {
    return(c(auc = NA, auc_lower = NA, auc_upper = NA, median_change = NA, p = NA, n))
  }
  w <- suppressWarnings(wilcox.test(s$pos, s$neg)$statistic)
  # direction = "<" pins HR as the higher side; pROC otherwise flips any AUC below 0.5. At perfect
  # separation DeLong's variance is zero and the interval collapses to a point, so none is given.
  interval <- if (spec$paired || w %in% c(0, prod(n))) {
    c(NA, NA)
  } else {
    pROC::roc(controls = s$neg, cases = s$pos, direction = "<", quiet = TRUE) |>
      pROC::ci.auc(method = "delong") |>
      as.numeric() |>
      _[c(1, 3)]
  }
  c(
    auc = unname(w) / prod(n), auc_lower = interval[1], auc_upper = interval[2],
    median_change = if (spec$paired) median(s$pos - s$neg) else NA,
    p = suppressWarnings(wilcox.test(s$pos, s$neg, paired = spec$paired)$p.value),
    n
  )
}

# cor.test(exact = FALSE)'s p, step for step: S statistic, rho back from S, then the t tail.
spearman_p <- function(rho, n) {
  s <- (n^3 - n) * (1 - rho) / 6
  r <- 1 - s / ((n * (n^2 - 1)) / 6)
  t <- r / sqrt((1 - r^2) / (n - 2))
  pmin(2 * if_else(s > (n^3 - n) / 6, pt(t, n - 2), pt(t, n - 2, lower.tail = FALSE)), 1)
}
# A feature is tested only when two thirds of the window's subjects saw it.
spearman_by_row <- function(x, outcome) {
  usable <- !is.na(outcome)
  x <- x[, usable, drop = FALSE]
  n <- rowSums(!is.na(x))
  rho <- suppressWarnings(
    cor(t(x), outcome[usable], method = "spearman", use = "pairwise.complete.obs")[, 1]
  )
  rho[n < ceiling(2 * sum(usable) / 3)] <- NA
  tibble(feature = rownames(x), n_subjects = unname(n), rho = unname(rho), p = spearman_p(rho, n))
}
# Exact p from every ordering of one variable against the other's ranks, ties kept; without ties
# it equals cor.test(exact = TRUE). Only used within arms, at 9 subjects or fewer.
orderings <- function(n) {
  if (n == 1) {
    return(matrix(1L))
  }
  shorter <- orderings(n - 1)
  do.call(rbind, map(seq_len(n), \(first) cbind(first, shorter + (shorter >= first))))
}
null_rho <- new.env()
exact_spearman_p <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  rx <- rank(x[ok])
  ry <- rank(y[ok])
  stopifnot(length(rx) <= 9)
  key <- paste(paste(sort(rx), collapse = ","), paste(sort(ry), collapse = ","))
  if (is.null(null_rho[[key]])) {
    shuffled <- matrix(rx[orderings(length(rx))], ncol = length(rx))
    centred <- ry - mean(ry)
    rho <- (shuffled - mean(rx)) %*% centred / sqrt(sum((rx - mean(rx))^2) * sum(centred^2))
    null_rho[[key]] <- sort(abs(rho[, 1]))
  }
  null <- null_rho[[key]]
  (length(null) - findInterval(abs(cor(rx, ry)) - 1e-9, null)) / length(null)
}

# limma on one value per subject, with the two-thirds floor: the chosen terms' effect, t and p.
subject_fit <- function(x, predictors, formula, terms) {
  ok <- complete.cases(predictors)
  x <- x[, ok, drop = FALSE]
  x <- x[rowSums(!is.na(x)) >= ceiling(2 * ncol(x) / 3), , drop = FALSE]
  fit <- lmFit(x, model.matrix(formula, predictors[ok, , drop = FALSE])) |>
    eBayes(robust = TRUE)
  map(terms, \(term) {
    tibble(
      feature = rownames(x), coefficient = term, n_subjects = unname(rowSums(!is.na(x))),
      effect = unname(fit$coefficients[, term]), t = unname(fit$t[, term]),
      p = unname(fit$p.value[, term])
    )
  }) |>
    list_rbind() |>
    filter(!is.na(p))
}


# Display only: row order, colour limits and centring touch no test.
ComplexHeatmap::ht_opt(message = FALSE)
diverging <- c("#5E3C99", "white", "#E66101")
status_colours <- c("not significant" = "white", "p < 0.05" = "#35978F", "BH < 0.05" = "black")
status_legend <- ComplexHeatmap::Legend(
  labels = c(names(status_colours), "not tested"), title = "status",
  legend_gp = grid::gpar(fill = c(status_colours, "grey80"))
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
# An effect map on a diverging scale when breaks are given, otherwise a status map. Row labels
# come from the calling script's feature_label.
heat <- function(m, title, breaks = NULL, legend = title, row_names = FALSE,
                 column_title = title, column_names_gp = gpar(fontsize = 7),
                 labels = feature_label, ...) {
  rows <- rownames(m)
  Heatmap(m,
    name = title, na_col = "grey80",
    col = if (is.null(breaks)) status_colours else circlize::colorRamp2(breaks, diverging),
    show_heatmap_legend = !is.null(breaks) && !is.null(legend),
    heatmap_legend_param = list(title = legend, at = if (!is.null(breaks)) signif(breaks, 2)),
    cluster_rows = FALSE, cluster_columns = FALSE, column_title = column_title,
    column_title_gp = gpar(fontsize = 8), column_names_gp = column_names_gp,
    show_row_names = row_names && length(rows) <= 80,
    row_labels = gsub("\n", " ", labels[rows]),
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


# Twelve panels to a page; blank slots keep a short last page at the same panel size.
slots <- \(panels) c(panels, strrep(" ", seq_len(12 - length(panels))))
paginate <- function(hits, by) {
  hits |>
    mutate(
      page = (row_number() - 1) %/% 12 + 1, n_pages = max(page),
      position = (row_number() - 1) %% 12 + 1, .by = all_of(by)
    ) |>
    mutate(pdf_page = consecutive_id(pick(all_of(by)), page))
}
page_theme <- ggplot2::theme_minimal(base_size = 9) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(size = 7),
    plot.title = ggplot2::element_text(face = "bold"),
    plot.caption = ggplot2::element_text(hjust = 0), legend.position = "top"
  )
counts_line <- \(info, noun, hits) {
  sprintf(
    "%d of %d %s at p < 0.05, %d at BH < 0.05. Page %d of %d.", info$n_nominal, info$n_tested,
    noun, info$n_fdr_05, hits$page[1], hits$n_pages[1]
  )
}


# Every workbook opens on read_me and ends with the files read and the packages loaded.
write_workbook <- function(sheets, holds, file, inputs, paths) {
  sheets <- c(sheets, list(
    input_manifest = tibble(
      input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
    ),
    package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
      as_tibble() |>
      select(package, version = loadedversion, source)
  ))
  holds <- c(holds, "Files read, with md5.", "Packages loaded at run time.")
  writexl::write_xlsx(c(list(read_me = tibble(sheet = names(sheets), holds = holds)), sheets), file)
}
