# The module eigengenes put through the classification tasks and phenotype association the
# proteins and sets answered: eight tasks mirroring the contrasts, and Spearman against the ten
# phenotypes in three windows. Twelve modules make BH a far weaker filter than 1,900 proteins;
# each table reports the count chance would give beside it. The contrasts are in 04_test_modules.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

stage <- here("04_Network", "05_classify_and_associate_modules")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
# Paging changes how many files a figure writes, so the previous run's figures go first.
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

inputs <- c(
  modules = "04_Network/01_build_modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_build_modules first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
modules <- readRDS(paths[["modules"]])
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)
me <- modules$eigengenes


# ---- the eight classification tasks --------------------------------------------------------

meta <- modules$meta
arm_of <- deframe(distinct(meta, subject, arm))
by_subject <- function(mat, timepoint) {
  rows <- filter(meta, timepoint == .env$timepoint)
  set_names(as.data.frame(mat[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
}
change <- function(mat, from, to) {
  a <- by_subject(mat, from)
  b <- by_subject(mat, to)
  both <- intersect(colnames(a), colnames(b))
  b[, both, drop = FALSE] - a[, both, drop = FALSE]
}
windows <- list(
  training = change(me, "T1", "T2"), baseline = by_subject(me, "T1"), acute = change(me, "T2", "T3")
)
within_arm <- function(from, to, arm, label) {
  a <- by_subject(me, from)
  b <- by_subject(me, to)
  both <- intersect(colnames(a), colnames(b))
  both <- both[arm_of[both] == arm]
  list(
    positive = b[, both, drop = FALSE], negative = a[, both, drop = FALSE], paired = TRUE,
    label = label
  )
}
between_arm <- function(values, label) {
  list(
    positive = values[, arm_of[colnames(values)] == "HR", drop = FALSE],
    negative = values[, arm_of[colnames(values)] == "LR", drop = FALSE], paired = FALSE,
    label = label
  )
}
tasks <- list(
  Training_HR = within_arm("T1", "T2", "HR", "Training, HR"),
  Training_LR = within_arm("T1", "T2", "LR", "Training, LR"),
  Acute_HR = within_arm("T2", "T3", "HR", "Acute bout, HR"),
  Acute_LR = within_arm("T2", "T3", "LR", "Acute bout, LR"),
  Baseline_HRvLR = between_arm(windows$baseline, "HR vs LR at T1"),
  Trained_HRvLR = between_arm(by_subject(me, "T2"), "HR vs LR at T2"),
  Training_change_HRvLR = between_arm(windows$training, "HR vs LR, training change"),
  Acute_change_HRvLR = between_arm(windows$acute, "HR vs LR, acute change")
)
fit_roc <- function(positive, negative) {
  pROC::roc(controls = negative, cases = positive, direction = "<", quiet = TRUE)
}
# direction = "<" pins pROC's orientation, so AUC above 0.5 means higher in the later timepoint
# or in HR.
module_auc <- imap(tasks, \(spec, name) {
  tibble(module = rownames(me), task = name, paired = spec$paired) |>
    mutate(
      auc = map_dbl(module, \(m) {
        as.numeric(pROC::auc(fit_roc(spec$positive[m, ], spec$negative[m, ])))
      }),
      p = map_dbl(module, \(m) {
        wilcox.test(spec$positive[m, ], spec$negative[m, ], paired = spec$paired)$p.value
      })
    )
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = task)


# ---- phenotype association -----------------------------------------------------------------

outcomes <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_fcsa_mixed", "d_nfibre_mixed",
  "d_nfibre_I", "d_mcsa", "d_1rm_legpress", "d_1rm_ext", "volume_load"
)
outcome_of <- function(name, subjects) deframe(select(phenotype, subject, all_of(name)))[subjects]
# The t approximation: above nine observations cor.test's "exact" Spearman p is an Edgeworth
# series that returns 0 in the far tail.
spearman_by_row <- function(values, outcome) {
  usable <- !is.na(outcome)
  fits <- apply(values[, usable, drop = FALSE], 1, \(row) {
    test <- cor.test(row, outcome[usable], method = "spearman", exact = FALSE)
    c(r = unname(test$estimate), p = test$p.value)
  })
  tibble(module = rownames(values), n = sum(usable), r = fits["r", ], p = fits["p", ])
}
module_association <- imap(windows, \(values, window) {
  map(set_names(outcomes), \(name) spearman_by_row(values, outcome_of(name, colnames(values)))) |>
    list_rbind(names_to = "outcome") |>
    mutate(window = window, .before = 1)
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(window, outcome))

by_arm <- map(set_names(c("HR", "LR")), \(arm) {
  imap(windows, \(values, window) {
    kept <- values[, arm_of[colnames(values)] == arm, drop = FALSE]
    map(set_names(outcomes), \(name) spearman_by_row(kept, outcome_of(name, colnames(kept)))) |>
      list_rbind(names_to = "outcome") |>
      mutate(window = window, .before = 1)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "arm")

chance_expectation <- bind_rows(
  summarise(module_auc,
    tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05),
    .by = task
  ) |>
    transmute(analysis = "classification", comparison = task, tested, nominal, fdr_sig),
  summarise(module_association,
    tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05), .by = c(window, outcome)
  ) |>
    transmute(
      analysis = paste0("association: ", window), comparison = outcome, tested,
      nominal, fdr_sig
    )
) |>
  mutate(expected = 0.05 * tested, ratio = round(nominal / expected, 2), .after = tested)
print(as.data.frame(
  chance_expectation |>
    summarise(across(c(tested, expected, nominal, fdr_sig), sum), .by = analysis) |>
    mutate(ratio = round(nominal / expected, 2))
))


# ---- figures -------------------------------------------------------------------------------

figure_theme <- theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 7.6, lineheight = 1.3, margin = margin(3, 3, 5, 3)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(5, "mm"),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40"),
    plot.tag = element_text(size = 9, colour = "grey30"),
    plot.tag.position = "topright",
    legend.position = "top"
  )
# One PDF holding every page; a single-page figure also gets a PNG.
save_pages <- function(pages, name, width, height) {
  if (length(pages) == 1) {
    ggsave(file.path(figure_dir, paste0(name, ".png")), pages[[1]],
      width = width, height = height, dpi = 220, bg = "white"
    )
  }
  pdf(file.path(figure_dir, paste0(name, ".pdf")), width = width, height = height, bg = "white")
  walk(pages, print)
  invisible(dev.off())
}
save_paged <- function(data, draw, name, scales = "fixed", panel = 2.6, header = 2.1,
                       strip = 0) {
  n_panels <- nlevels(data$panel)
  columns <- min(4, n_panels)
  rows <- min(4, ceiling(n_panels / columns))
  page_of <- (as.integer(data$panel) - 1) %/% (columns * rows) + 1
  n_pages <- max(page_of)
  pages <- map(seq_len(n_pages), \(page) {
    draw(droplevels(data[page_of == page, ])) +
      facet_wrap(~panel, ncol = columns, nrow = rows, scales = scales) +
      labs(tag = if (n_pages > 1) sprintf("page %d of %d", page, n_pages)) +
      figure_theme
  })
  save_pages(pages, name, width = 1.2 + columns * panel, height = header + rows * (panel + strip))
  message("drew ", name, ": ", n_panels, " panels on ", n_pages, " pages")
  n_panels
}
module_order <- modules$module_summary$module

# One tile per module and column, every module drawn; the label marks nominal p and a black
# border BH < 0.05.
module_tiles <- function(data, columns, fill_label, midpoint, limits, title, subtitle, caption) {
  data |>
    mutate(
      module = factor(module, levels = rev(module_order)),
      column = factor(column, levels = columns)
    ) |>
    ggplot(aes(column, module, fill = effect)) +
    geom_tile(colour = "white") +
    geom_tile(data = \(x) filter(x, fdr < 0.05), colour = "black", linewidth = 0.8, fill = NA) +
    geom_text(aes(label = if_else(p < 0.05, sprintf("%.3f", p), "")), size = 2.6) +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = midpoint, limits = limits
    ) +
    labs(
      x = NULL, y = NULL, fill = fill_label, title = title, subtitle = subtitle,
      caption = paste(
        caption, "Label: nominal p where below 0.05. Black border: BH < 0.05.",
        "Table: c_data/05_classify_and_associate_modules.xlsx."
      )
    ) +
    figure_theme +
    theme(axis.text.x = element_text(angle = 35, hjust = 1))
}
save_pages(list(module_tiles(
  rename(module_auc, column = task, effect = auc), names(tasks), "AUC", 0.5, c(0, 1),
  "Module classification, every module and task",
  "pROC AUC, Wilcoxon p (signed-rank within arm), BH within task",
  "Fill: AUC, above 0.5 higher at the later timepoint or in HR."
)), "01_classification_tiles", width = 9, height = 5.5)
window_titles <- c(
  training = "training change (T2 - T1)", baseline = "level at T1",
  acute = "acute change (T3 - T2)"
)
association_tiles <- imap(window_titles, \(title, window) {
  module_tiles(
    module_association |>
      filter(window == .env$window) |>
      rename(column = outcome, effect = r),
    outcomes, "rho", 0, c(-1, 1), paste("Module association,", title),
    "Spearman, subjects pooled, BH within outcome",
    "Fill: Spearman rho between eigengene and phenotype."
  )
})
save_pages(unname(association_tiles), "02_association_tiles", width = 9, height = 5.5)

# An ROC panel for every module reaching nominal p on a task, the floor excluded.
roc_hits <- module_auc |>
  filter(p < 0.05, task != "Baseline_HRvLR") |>
  mutate(task = factor(task, levels = names(tasks))) |>
  arrange(task, p)
if (nrow(roc_hits)) {
  curves <- pmap(roc_hits, function(module, task, auc, p, fdr, ...) {
    spec <- tasks[[as.character(task)]]
    coordinates <- pROC::coords(fit_roc(spec$positive[module, ], spec$negative[module, ]), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = paste0(
        module, ": ", spec$label,
        "\nAUC ", sprintf("%.2f", auc), "   p ", signif(p, 2), "   q ", signif(fdr, 2)
      )
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel))) |>
    arrange(panel, fpr, tpr)
  invisible(save_paged(curves, \(page) {
    ggplot(page, aes(fpr, tpr, colour = direction, fill = direction)) +
      geom_abline(linetype = "22", linewidth = 0.35, colour = "grey60") +
      geom_ribbon(aes(ymin = 0, ymax = tpr), alpha = 0.16, colour = NA) +
      geom_step(linewidth = 0.8, direction = "hv") +
      scale_colour_manual(values = c(higher = "#B2182B", lower = "#2166AC"), guide = "none") +
      scale_fill_manual(values = c(higher = "#B2182B", lower = "#2166AC"), guide = "none") +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      labs(
        x = "1 - specificity", y = "sensitivity", title = "Module classification, nominal hits",
        subtitle = "eigengene per module, AUC with direction fixed",
        caption = stringr::str_wrap(width = 150, paste(
          "Every module reaching nominal p on a task, by task then p; the floor is excluded.",
          "Red separates higher at the later timepoint or in HR, blue lower. p from the Wilcoxon",
          "test (signed-rank within arm), q from BH over the twelve modules within the task."
        ))
      )
  }, "03_roc_nominal", strip = 0.6))
}

# A scatter for every module-outcome pair reaching nominal p, with the within-arm correlations.
association_hits <- module_association |>
  filter(p < 0.05) |>
  mutate(window = factor(window, levels = names(window_titles))) |>
  arrange(window, p)
if (nrow(association_hits)) {
  points <- pmap(association_hits, function(window, outcome, module, n, r, p, fdr, ...) {
    window <- as.character(window)
    values <- windows[[window]]
    arms <- by_arm |>
      filter(module == !!module, outcome == !!outcome, window == !!window) |>
      mutate(text = sprintf("%-3s n=%d  r=%+.2f  p=%.3f", arm, n, r, p))
    tibble(
      panel = paste0(
        module, " vs ", outcome, ", ", window,
        "\nr = ", sprintf("%+.2f", r), "   p ", signif(p, 2), "   q ", signif(fdr, 2)
      ),
      score = values[module, ], outcome_value = outcome_of(outcome, colnames(values)),
      arm = arm_of[colnames(values)], caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    filter(!is.na(outcome_value)) |>
    mutate(panel = factor(panel, levels = unique(panel)))
  invisible(save_paged(points, \(page) {
    ggplot(page, aes(score, outcome_value, colour = arm, fill = arm)) +
      geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
      geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.12, linewidth = 0.6) +
      geom_point(size = 1.9, alpha = 0.9) +
      geom_text(
        data = distinct(page, panel, caption), inherit.aes = FALSE,
        aes(x = -Inf, y = Inf, label = caption), family = "mono",
        hjust = -0.05, vjust = 1.25, size = 2.3, lineheight = 1.25, colour = "grey25"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.3))) +
      scale_colour_manual(values = c(HR = "#2166AC", LR = "#B2182B"), name = NULL) +
      scale_fill_manual(values = c(HR = "#2166AC", LR = "#B2182B"), name = NULL) +
      labs(
        x = "eigengene, one point per subject", y = "phenotype",
        title = "Module eigengene against phenotype, nominal hits",
        subtitle = "Spearman, subjects pooled across arms",
        caption = stringr::str_wrap(width = 150, paste(
          "Every module-outcome pair reaching nominal p, by window then p. Lines are fitted within",
          "each arm; the header r is the pooled correlation, q from BH over the twelve modules",
          "within window and outcome."
        ))
      )
  }, "04_association_nominal", scales = "free", panel = 3.1, header = 2.4))
}

save_pages(list(
  chance_expectation |>
    summarise(across(c(tested, expected, nominal), sum), .by = analysis) |>
    mutate(ratio = nominal / expected, analysis = factor(analysis, rev(unique(analysis)))) |>
    ggplot(aes(ratio, analysis, fill = ratio > 1)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%d of %d", nominal, tested)),
      hjust = -0.12, size = 2.8, colour = "grey25"
    ) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(
      x = "observed nominal hits / chance expectation", y = NULL,
      title = "Module nominal hits relative to chance",
      subtitle = "twelve modules, uncorrected p",
      caption = paste(
        "Bar length is nominal hits over 5% of tests. Red clears 1, grey does not.",
        "Table: c_data/05_classify_and_associate_modules.xlsx, chance_expectation."
      )
    ) +
    figure_theme +
    theme(panel.grid.major.y = element_blank())
), "05_chance", width = 7, height = 3.5)


packages <- c("here", "pROC", "dplyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    auc = module_auc, association = module_association, by_arm = by_arm,
    chance_expectation = chance_expectation,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "module_results.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    chance_expectation = chance_expectation, classification = arrange(module_auc, p),
    association = arrange(module_association, p), by_arm = by_arm,
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "05_classify_and_associate_modules.xlsx")
)
combined <- file.path(figure_dir, "05_classify_and_associate_modules_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote module_results.rds, 05_classify_and_associate_modules.xlsx and ", length(pages),
  " figures bundled into ", qpdf::pdf_length(combined), " pages"
)
