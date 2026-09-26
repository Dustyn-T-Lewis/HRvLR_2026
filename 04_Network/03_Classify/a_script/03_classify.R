# How well each module eigengene separates the study's groups, on the protein level's eight tasks.
# Within-arm tasks pair each subject with itself; between-arm tasks compare HR with LR subjects on
# a level or on their own change. AUC is the rank-sum statistic over n1 * n2, so above 0.5 means
# higher at the later timepoint or in HR. p is Wilcoxon, signed-rank when paired; BH within task.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
})

stage <- here("04_Network", "03_Classify")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(modules = "04_Network/01_Build_Modules/c_data/modules.rds")
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
modules <- readRDS(paths[["modules"]])
values <- modules$eigengenes
meta <- modules$meta
stopifnot(identical(colnames(values), meta$sample_id))
feature_label <- set_names(rownames(values))

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
  list(positive = p$after, negative = p$before, paired = TRUE, label = label, favours = to)
}
between_arm <- function(level, label) {
  hr <- arm_of[colnames(level)] == "HR"
  list(
    positive = level[, hr, drop = FALSE], negative = level[, !hr, drop = FALSE],
    paired = FALSE, label = label, favours = "HR"
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

# Ties push wilcox.test to its normal approximation and it warns on every call.
classify_row <- function(pos, neg, paired) {
  kept <- if (paired) !is.na(pos) & !is.na(neg) else NULL
  pos <- if (paired) pos[kept] else pos[!is.na(pos)]
  neg <- if (paired) neg[kept] else neg[!is.na(neg)]
  if (min(length(pos), length(neg)) < 3) {
    return(c(auc = NA, p = NA, n_positive = length(pos), n_negative = length(neg)))
  }
  w <- suppressWarnings(wilcox.test(pos, neg)$statistic)
  p <- suppressWarnings(wilcox.test(pos, neg, paired = paired)$p.value)
  c(auc = unname(w) / (length(pos) * length(neg)), p = p, n_positive = length(pos),
    n_negative = length(neg))
}
auc <- imap(tasks, \(spec, name) {
  fits <- vapply(rownames(values), \(id) {
    classify_row(spec$positive[id, ], spec$negative[id, ], spec$paired)
  }, numeric(4))
  tibble(task = name, task_label = spec$label, feature = rownames(values), !!!as_tibble(t(fits)))
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = task) |>
  mutate(label = feature_label[feature], .after = feature)

# The smallest p a task can reach: signed-rank on n pairs gives 2 / 2^n, rank-sum on n1 and n2
# gives 2 / choose(n1 + n2, n1). HR has 6 training and 7 acute pairs, LR 8 of each.
chance_expectation <- auc |>
  filter(!is.na(p)) |>
  summarise(
    n_positive = max(n_positive), n_negative = max(n_negative), n_tested = n(),
    n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05), .by = c(task, task_label)
  ) |>
  mutate(
    paired = map_lgl(task, \(t) tasks[[t]]$paired),
    min_p = if_else(paired, 2 / 2^n_positive, 2 / choose(n_positive + n_negative, n_positive)),
    n_expected = 0.05 * n_tested, ratio = round(n_nominal / n_expected, 2), .after = n_nominal
  )

# One ROC curve per nominal feature, twelve to a page, by p. The floor task is not drawn.
roc_of <- function(spec, id) {
  pos <- spec$positive[id, ]
  neg <- spec$negative[id, ]
  if (spec$paired) {
    kept <- !is.na(pos) & !is.na(neg)
    pos <- pos[kept]
    neg <- neg[kept]
  }
  pROC::roc(controls = na.omit(neg), cases = na.omit(pos), direction = "<", quiet = TRUE)
}
roc_pages <- auc |>
  filter(p < 0.05, task != "Baseline_HRvLR") |>
  mutate(task = factor(task, names(tasks))) |>
  arrange(task, p) |>
  mutate(
    panel = sprintf("%s\nAUC %.2f  p %.2g  q %.2g", label, auc, p, fdr),
    page = (row_number() - 1) %/% 12, .by = task
  ) |>
  group_by(task, page) |>
  group_map(\(hits, key) {
    spec <- tasks[[as.character(key$task)]]
    info <- filter(chance_expectation, task == as.character(key$task))
    curves <- set_names(map(hits$feature, \(id) roc_of(spec, id)), hits$panel)
    pROC::ggroc(curves, aes = "group", legacy.axes = TRUE, colour = "#B2182B") +
      geom_abline(linetype = "22", colour = "grey60") +
      facet_wrap(~ factor(name, levels = hits$panel), nrow = 3, ncol = 4) +
      coord_equal() +
      labs(
        title = spec$label, x = "1 - specificity", y = "sensitivity",
        subtitle = sprintf(
          "eigengene; AUC above 0.5 favours %s; %d vs %d %s; page %d of %d",
          spec$favours, info$n_positive, info$n_negative,
          if (spec$paired) "paired samples" else "subjects", key$page + 1,
          ceiling(info$n_nominal / 12)
        ),
        caption = sprintf(
          paste(
            "%d of %d modules reach p < 0.05; chance predicts %.1f. Smallest attainable p: %.2g.",
            "p from the Wilcoxon %s test, q from BH within task.",
            "Table: c_data/03_classify.xlsx, auc."
          ),
          info$n_nominal, info$n_tested, info$n_expected, info$min_p,
          if (spec$paired) "signed-rank" else "rank-sum"
        )
      ) +
      theme_minimal(base_size = 9) +
      theme(strip.text = element_text(size = 7), plot.caption = element_text(hjust = 0))
  })
n_drawn <- sum(map_int(roc_pages, \(g) n_distinct(g$data$name)))
stopifnot(n_drawn == sum(auc$p < 0.05 & auc$task != "Baseline_HRvLR", na.rm = TRUE))

chance_figure <- chance_expectation |>
  mutate(task_label = factor(task_label, rev(map_chr(tasks, "label")))) |>
  ggplot(aes(ratio, task_label)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_col(fill = "grey60", width = 0.6) +
  geom_text(aes(label = sprintf("%d of %d, min p %.2g", n_nominal, n_tested, min_p)),
    hjust = -0.05, size = 3
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.4))) +
  labs(
    title = "Nominal classifiers against chance", x = "nominal hits / chance expectation",
    y = NULL, caption = paste(
      "Chance is 5% of the 12 modules tested. With 6 to 8 pairs the Wilcoxon p is coarse, so fewer",
      "than 5% of null tests reach 0.05 and a ratio below 1 is not worse than chance."
    )
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.caption = element_text(hjust = 0))

pdf(file.path(figures, "03_classify_figures.pdf"), width = 11, height = 8.5)
walk(c(list(chance_figure), roc_pages), print)
invisible(dev.off())

sheets <- list(
  chance_expectation = chance_expectation,
  auc = arrange(auc, task, p),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Per task: group sizes, smallest attainable p, tested, nominal, chance-expected and BH counts.",
  "AUC, Wilcoxon p and BH FDR per module and task.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "03_classify.xlsx"))
