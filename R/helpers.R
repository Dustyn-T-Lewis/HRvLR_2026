# Helpers every stage script sources after its library() calls: the statistics written for this
# study, and the heatmap, page and workbook plumbing the contrasts, classify and associate steps
# share.

# Blood-corrected values for the methods that cannot carry a covariate: each feature's blood slope,
# fitted beside the given design, is removed (limma::removeBatchEffect). The design protects what
# the analysis compares; pass one without arm when the labels must not inform the correction.
blood_corrected <- function(values, blood, design) {
  limma::removeBatchEffect(values, covariates = blood, design = design)
}

# limpa's fit on its row-wise quantification: precision weights from the standard errors, subject
# blocked, then the contrasts with robust moderation. Every protein-level contrast uses it.
fit_limpa <- function(proteins, design, contrasts, subject) {
  limpa::dpcDE(proteins, design, block = subject, plot = FALSE) |>
    contrasts.fit(contrasts) |>
    eBayes(robust = TRUE)
}

# Storey's q-value within one screen, beside BH as an exploratory correction: BH scaled by the
# estimated share of true nulls (pi0, bootstrap estimate, which copes with discrete p-values). A
# screen of fewer than 100 tests cannot estimate pi0, so there pi0 is 1 and q equals BH.
storey_q <- function(p) {
  q <- rep(NA_real_, length(p))
  tested <- !is.na(p)
  q[tested] <- if (sum(tested) < 100) {
    p.adjust(p[tested], "BH")
  } else {
    qvalue::qvalue(p[tested], pi0.method = "bootstrap")$qvalues
  }
  q
}

# The four contrasts that compare HR with LR in the candidate tables.
hr_lr_contrasts <- c("Training_Interaction", "Trained_HRvLR", "Acute_Interaction", "Acute_HRvLR")

# limma on a complete matrix (module eigengenes) with subject blocked and robust moderation.
fit_contrasts <- function(values, design, contrasts, subject) {
  correlation <- duplicateCorrelation(values, design, block = subject)$consensus.correlation
  lmFit(values, design, block = subject, correlation = correlation) |>
    contrasts.fit(contrasts) |>
    eBayes(robust = TRUE)
}
# One row per feature and contrast from a contrasts fit: logFC, t and p.
contrast_rows <- function(fit, prefix = "") {
  tibble(
    feature = rep(rownames(fit$coefficients), ncol(fit$coefficients)),
    test = rep(paste0(prefix, colnames(fit$coefficients)), each = nrow(fit$coefficients)),
    effect = c(fit$coefficients), t = c(fit$t), p = c(fit$p.value)
  )
}

# fry per contrast on a complete matrix, sets given as row indices: one row per set and contrast.
fry_rows <- function(values, index, design, contrasts, subject) {
  correlation <- duplicateCorrelation(values, design, block = subject)$consensus.correlation
  map(colnames(contrasts), \(contrast) {
    fry(values, index, design, contrasts[, contrast],
      block = subject, correlation = correlation, sort = "none"
    ) |>
      rownames_to_column("feature") |>
      transmute(feature, test = paste("fry", contrast), p = PValue)
  }) |>
    list_rbind()
}

# Reruns test() once per subject left out; test(keep) returns feature, test and p for the kept
# subjects. Per feature and test: reruns made, reruns still at p < 0.05, the largest p, and the
# subject whose removal gave it.
leave_one_out <- function(subjects, test) {
  map(set_names(subjects), \(dropped) test(setdiff(subjects, dropped))) |>
    list_rbind(names_to = "dropped") |>
    filter(!is.na(p)) |>
    summarise(
      n_reruns = n(), n_nominal = sum(p < 0.05), max_p = max(p), worst = dropped[which.max(p)],
      .by = c(feature, test)
    )
}

# One column per subject at a time point, and the change between two time points over the
# subjects sampled at both.
subject_values <- function(values, meta, timepoint) {
  rows <- filter(meta, timepoint == .env$timepoint)
  set_names(as.data.frame(values[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
}
subject_change <- function(values, meta, from, to) {
  before <- subject_values(values, meta, from)
  after <- subject_values(values, meta, to)
  both <- intersect(colnames(before), colnames(after))
  after[, both, drop = FALSE] - before[, both, drop = FALSE]
}

# Which tests feed which candidate list (docs/decisions.md, 2026-09-29); Acute_Interaction feeds
# both. The baseline list applies the same rules at T1, before training.
candidate_lists <- list(
  training = list(
    contrasts = c("Training_Interaction", "Trained_HRvLR", "Acute_Interaction"),
    ancova = "training", tasks = c("Training_change_HRvLR", "Trained_HRvLR"), window = "training"
  ),
  acute = list(
    contrasts = c("Acute_Interaction", "Acute_HRvLR"),
    ancova = "acute", tasks = "Acute_change_HRvLR", window = "acute"
  ),
  baseline = list(
    contrasts = "Baseline_HRvLR", ancova = character(0), tasks = "Baseline_HRvLR",
    window = "baseline"
  )
)
# The rows of a long evidence table (columns contrast, task, window) that fall on each list.
on_lists <- function(evidence) {
  imap(candidate_lists, \(rule, name) {
    evidence |>
      filter(
        contrast %in% rule$contrasts | window %in% c(paste0("ancova_", rule$ancova), rule$window) |
          task %in% rule$tasks
      ) |>
      mutate(list = name)
  }) |>
    list_rbind()
}

# ANCOVA per feature and window: the later value on the earlier one plus arm, over subjects with
# both biopsies and at least three per arm. The arm term is HR minus LR at equal starting value.
# Residual variances are moderated across features within each window (limma::squeezeVar), as
# limma moderates its own fits, so a feature that looks quiet by chance does not earn a large t.
ancova <- function(values, meta) {
  windows <- list(training = c("T1", "T2"), acute = c("T2", "T3"))
  imap(windows, \(timepoints, window) {
    before <- subject_values(values, meta, timepoints[1])
    after <- subject_values(values, meta, timepoints[2])
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
      fit <- summary(lm(y[ok] ~ x[ok] + arm[ok]))
      tibble(
        feature, window,
        n_HR = n[["HR"]], n_LR = n[["LR"]], effect = fit$coefficients[3, 1],
        unscaled = fit$coefficients[3, 2] / fit$sigma, s2 = fit$sigma^2, df = fit$df[2]
      )
    }) |>
      list_rbind()
  }) |>
    list_rbind() |>
    split(~window) |>
    map(\(rows) {
      squeezed <- squeezeVar(rows$s2, rows$df)
      # As eBayes does, total df never exceed the pooled residual df.
      df_total <- pmin(rows$df + squeezed$df.prior, sum(rows$df))
      mutate(rows,
        t = effect / (unscaled * sqrt(squeezed$var.post)), p = 2 * pt(-abs(t), df_total)
      )
    }) |>
    list_rbind() |>
    select(feature, window, n_HR, n_LR, effect, t, p)
}

# Exploratory elastic-net panel: repeated stratified k-fold CV around nestedcv's glmnet, alpha 0.5
# fixed and lambda.1se from an inner leave-one-out. Stratified outer folds keep the HR:LR balance
# of every training set; pooled leave-one-out does not, and an empty model then scores AUC 0
# (Airola et al. 2011). x has one row per subject. Returns each repeat's out-of-fold AUC, how
# often each feature was kept across the repeat-by-fold models, and the share that kept none.
panel_cv <- function(x, y, repeats = 10, k = 5, seed = 1) {
  inner <- length(y) - ceiling(length(y) / k)
  runs <- map(seq_len(repeats), \(r) {
    set.seed(seed + r)
    fit <- suppressWarnings(nestedcv::nestcv.glmnet(y, x,
      family = "binomial", outer_folds = caret::createFolds(y, k = k), n_inner_folds = inner,
      alphaSet = 0.5, min_1se = 1, finalCV = FALSE
    ))
    list(
      auc = pROC::roc(fit$output$testy, fit$output$predyp,
        levels = levels(y), direction = "<", quiet = TRUE
      ) |>
        pROC::auc() |>
        as.numeric(),
      kept = map(fit$outer_result, \(fold) setdiff(names(fold$coef), "(Intercept)"))
    )
  })
  kept <- list_flatten(map(runs, "kept"))
  list(
    auc = map_dbl(runs, "auc"),
    frequency = table(factor(unlist(kept), colnames(x))) / length(kept),
    empty = mean(lengths(kept) == 0)
  )
}
# Permutation p that is never zero: (b + 1) / (m + 1) (Phipson & Smyth 2010).
permutation_p <- function(observed, null) (sum(null >= observed) + 1) / (length(null) + 1)

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
# A feature is tested only when two thirds of the window's subjects saw it. exact = TRUE swaps in
# the exact p, for within-arm tests at nine subjects or fewer.
spearman_by_row <- function(x, outcome, exact = FALSE) {
  usable <- !is.na(outcome)
  x <- x[, usable, drop = FALSE]
  n <- rowSums(!is.na(x))
  rho <- suppressWarnings(
    cor(t(x), outcome[usable], method = "spearman", use = "pairwise.complete.obs")[, 1]
  )
  rho[n < ceiling(2 * sum(usable) / 3)] <- NA
  p <- spearman_p(rho, n)
  if (exact) {
    tested <- which(!is.na(rho))
    p[tested] <- map_dbl(tested, \(i) exact_spearman_p(x[i, ], outcome[usable]))
  }
  tibble(feature = rownames(x), n_subjects = unname(n), rho = unname(rho), p = unname(p))
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
