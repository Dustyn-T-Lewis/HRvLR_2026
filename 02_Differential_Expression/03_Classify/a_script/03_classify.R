# How well each protein separates the study's groups, on eight tasks that mirror the contrasts.
# Within-arm tasks pair each subject's two biopsies: p is the Wilcoxon signed-rank test and the
# effect is the median paired change. Between-arm tasks compare HR with LR subjects on a level or
# on their own change: p is the rank-sum test and AUC is W / (n1 * n2), above 0.5 when higher in
# HR. The auc column is filled for every task; for paired tasks it ignores the pairing and is not
# drawn. BH within task.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
})

stage <- here("02_Differential_Expression", "03_Classify")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(design = "02_Differential_Expression/01_Design/c_data/design.rds")
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
dal <- readRDS(paths[["design"]])$dal
values <- as.matrix(dal$data)
meta <- dal$metadata
feature_label <- dal$annotation |>
  mutate(label = if_else(
    is.na(gene) | duplicated(gene) | duplicated(gene, fromLast = TRUE),
    paste0(coalesce(gene, "?"), " (", uniprot_id, ")"), gene
  )) |>
  with(set_names(label, uniprot_id))
feature_noun <- "proteins"
value_label <- "log2 abundance"

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
  Baseline_HRvLR = between_arm(by_subject("T1"), "HR vs LR at T1 (floor)"),
  Trained_HRvLR = between_arm(by_subject("T2"), "HR vs LR at T2"),
  Training_change_HRvLR = between_arm(change("T1", "T2"), "HR vs LR, training change"),
  Acute_change_HRvLR = between_arm(change("T2", "T3"), "HR vs LR, acute change")
)

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
    return(c(auc = NA, median_change = NA, p = NA, n))
  }
  w <- suppressWarnings(wilcox.test(s$pos, s$neg)$statistic)
  c(
    auc = unname(w) / prod(n),
    median_change = if (spec$paired) median(s$pos - s$neg) else NA,
    p = suppressWarnings(wilcox.test(s$pos, s$neg, paired = spec$paired)$p.value),
    n
  )
}
auc <- imap(tasks, \(spec, name) {
  fits <- vapply(rownames(values), \(id) classify_row(spec, id), numeric(5))
  tibble(
    task = name, task_label = spec$label, paired = spec$paired, feature = rownames(values),
    !!!as_tibble(t(fits))
  )
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = task) |>
  mutate(label = feature_label[feature], .after = feature)

# The smallest p a task can reach: signed-rank on n pairs gives 2 / 2^n, rank-sum on n1 and n2
# gives 2 / choose(n1 + n2, n1). HR has 6 training and 7 acute pairs, LR 8 of each.
counts <- auc |>
  filter(!is.na(p)) |>
  summarise(
    paired = first(paired), n_positive = max(n_positive), n_negative = max(n_negative),
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(task, task_label)
  ) |>
  mutate(
    min_p = if_else(paired, 2 / 2^n_positive, 2 / choose(n_positive + n_negative, n_positive)),
    .after = n_negative
  )

page_theme <- theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 7), plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0)
  )
page_labels <- function(info, page, n_pages, caption) {
  labs(
    title = info$task_label,
    subtitle = sprintf(
      "%d of %d %s at p < 0.05, %d at BH < 0.05. Page %d of %d.",
      info$n_nominal, info$n_tested, feature_noun, info$n_fdr_05, page, n_pages
    ),
    caption = caption
  )
}
# Paired tasks: each subject's two biopsies joined by a line.
paired_page <- function(hits, spec, info, page, n_pages) {
  points <- map2(hits$feature, hits$panel, \(id, panel) {
    tibble(
      panel = panel, subject = colnames(spec$positive),
      before = spec$negative[id, ], after = spec$positive[id, ]
    )
  }) |>
    list_rbind() |>
    filter(!is.na(before), !is.na(after)) |>
    pivot_longer(c(before, after), names_to = "biopsy") |>
    mutate(
      biopsy = factor(if_else(biopsy == "before", spec$levels[1], spec$levels[2]), spec$levels),
      panel = factor(panel, hits$panel)
    )
  ggplot(points, aes(biopsy, value, group = subject)) +
    geom_line(colour = "grey55") +
    geom_point(size = 1.4) +
    facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free_y") +
    labs(x = NULL, y = value_label) +
    page_labels(info, page, n_pages, sprintf(
      paste(
        "Wilcoxon signed-rank test on %d pairs; smallest attainable p %.2g. Header: median",
        "paired change, p, and BH q within task. Table: c_data/03_classify.xlsx, auc."
      ),
      info$n_positive, info$min_p
    )) +
    page_theme
}
# Between-arm tasks: one ROC curve per feature, HR as cases.
roc_page <- function(hits, spec, info, page, n_pages) {
  curves <- map(hits$feature, \(id) {
    s <- sides(spec, id)
    pROC::roc(controls = s$neg, cases = s$pos, direction = "<", quiet = TRUE)
  }) |>
    set_names(hits$panel)
  pROC::ggroc(curves, aes = "group", legacy.axes = TRUE, colour = "#B2182B") +
    geom_abline(linetype = 2, colour = "grey60") +
    facet_wrap(~ factor(name, levels = hits$panel), nrow = 3, ncol = 4) +
    coord_equal() +
    labs(x = "1 - specificity", y = "sensitivity") +
    page_labels(info, page, n_pages, sprintf(
      paste(
        "HR (%d) against LR (%d); AUC above 0.5 means higher in HR. Wilcoxon rank-sum test;",
        "smallest attainable p %.2g. Header: AUC, p, and BH q within task.",
        "Table: c_data/03_classify.xlsx, auc."
      ),
      info$n_positive, info$n_negative, info$min_p
    )) +
    page_theme
}

hits <- auc |>
  filter(p < 0.05, task != "Baseline_HRvLR") |>
  mutate(
    task = factor(task, names(tasks)),
    effect = if_else(paired, sprintf("median change %+.2f", median_change),
      sprintf("AUC %.2f", auc)
    ),
    panel = sprintf("%s\n%s  p %.2g  q %.2g", label, effect, p, fdr)
  ) |>
  arrange(task, p) |>
  mutate(page = (row_number() - 1) %/% 12 + 1, n_pages = max(page), .by = task)
pages <- hits |>
  group_by(task, page) |>
  group_map(\(page_hits, key) {
    spec <- tasks[[as.character(key$task)]]
    info <- filter(counts, task == as.character(key$task))
    draw <- if (spec$paired) paired_page else roc_page
    draw(page_hits, spec, info, key$page, page_hits$n_pages[1])
  })
stopifnot(length(pages) == nrow(distinct(hits, task, page)))

pdf(file.path(figures, "03_classify_figures.pdf"), width = 11, height = 8.5)
walk(pages, print)
invisible(dev.off())

sheets <- list(
  counts = counts,
  auc = arrange(auc, task, p),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Per task: pairing, n per side, smallest attainable p, tested, p < 0.05 and BH < 0.05 counts.",
  "Per protein and task: AUC, median paired change (paired tasks), Wilcoxon p and BH FDR.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "03_classify.xlsx"))
