# For every tested set: how well its score separates the study's groups, and whether it tracks
# the phenotype. The unit is the set, with no collapse and no grouping,
# and results are read per database so each collection carries its own chance expectation.
#
# Eight tasks mirror the contrasts. The four within-arm tasks compare a subject's later biopsy
# with its earlier one, so they are paired. The four between-arm tasks compare HR subjects with LR
# subjects on a level or on their own change, so they are not. AUC comes from pROC as the
# effect-size descriptor; the p-value from the Wilcoxon test, signed-rank when paired.
#
# Nominal p is read against chance_expectation, the count a table that size returns under the
# null, per database. BH within database and task is reported beside it. Baseline_HRvLR is the
# floor.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

stage <- here("03_Pathway_Enrichment", "05_classify_and_associate_sets")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
# Paging changes how many files a figure writes, so the previous run's figures go first.
unlink(list.files(figure_dir, "[.](png|pdf)$", full.names = TRUE))

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds",
  singscore = "03_Pathway_Enrichment/04_run_singscore/c_data/singscore.rds",
  proteins = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 00 through 04 first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
set_catalog <- readRDS(paths[["gene_sets"]])$set_catalog
fgsea_results <- readRDS(paths[["set_tests"]])$set_tests |> filter(method == "fgsea")
set_score <- readRDS(paths[["singscore"]])$scores
metadata <- readRDS(paths[["proteins"]])$metadata
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

catalog <- set_catalog |>
  filter(qualifies) |>
  select(set_id, database, pathway, measured_size)
stopifnot(setequal(catalog$set_id, rownames(set_score)))
set_score <- set_score[catalog$set_id, ]
collection_sizes <- count(catalog, database)
message(nrow(catalog), " sets across ", nrow(collection_sizes), " collections")


# ---- subject-level matrices ----------------------------------------------------------------

meta <- as_tibble(metadata) |>
  select(sample_id = Col_ID, subject = Subject_ID, arm = Group, timepoint = Timepoint)
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
  training = change(set_score, "T1", "T2"),
  baseline = by_subject(set_score, "T1"),
  acute = change(set_score, "T2", "T3")
)


# ---- the eight classification tasks --------------------------------------------------------

# Each task carries its two matrices, columns aligned by subject when paired. `favours` is the
# group an AUC above 0.5 points to.
within_arm <- function(from, to, arm, label) {
  a <- by_subject(set_score, from)
  b <- by_subject(set_score, to)
  both <- intersect(colnames(a), colnames(b))
  both <- both[arm_of[both] == arm]
  list(
    positive = b[, both, drop = FALSE], negative = a[, both, drop = FALSE], paired = TRUE,
    label = label, favours = to, unit = "paired subjects"
  )
}
between_arm <- function(values, label) {
  list(
    positive = values[, arm_of[colnames(values)] == "HR", drop = FALSE],
    negative = values[, arm_of[colnames(values)] == "LR", drop = FALSE],
    paired = FALSE, label = label, favours = "HR", unit = "subjects"
  )
}
tasks <- list(
  Training_HR = within_arm("T1", "T2", "HR", "Training, HR (T1 to T2)"),
  Training_LR = within_arm("T1", "T2", "LR", "Training, LR (T1 to T2)"),
  Acute_HR = within_arm("T2", "T3", "HR", "Acute bout, HR (T2 to T3)"),
  Acute_LR = within_arm("T2", "T3", "LR", "Acute bout, LR (T2 to T3)"),
  Baseline_HRvLR = between_arm(windows$baseline, "HR vs LR at T1 (floor)"),
  Trained_HRvLR = between_arm(by_subject(set_score, "T2"), "HR vs LR at T2"),
  Training_change_HRvLR = between_arm(windows$training, "HR vs LR, training change"),
  Acute_change_HRvLR = between_arm(windows$acute, "HR vs LR, acute change")
)

# pROC::roc() auto-orients by default, which flips every below-chance set. direction = "<" pins
# it, so an AUC above 0.5 always means higher in the favoured group.
fit_roc <- function(positive, negative) {
  pROC::roc(controls = negative, cases = positive, direction = "<", quiet = TRUE)
}
classify <- function(spec) {
  tibble(
    set_id = rownames(spec$positive), n_positive = ncol(spec$positive),
    n_negative = ncol(spec$negative)
  ) |>
    mutate(
      auc = map_dbl(seq_along(set_id), \(i) {
        as.numeric(pROC::auc(fit_roc(spec$positive[i, ], spec$negative[i, ])))
      }),
      p = map_dbl(seq_along(set_id), \(i) {
        suppressWarnings(stats::wilcox.test(
          spec$positive[i, ], spec$negative[i, ],
          paired = spec$paired
        )$p.value)
      })
    )
}
# The argument is `spec`, not `task`: mutate() creates a `task` column first, and a later
# `task$label` in the same call would read that column instead of the argument.
set_auc <- imap(tasks, function(spec, name) {
  classify(spec) |>
    mutate(task = name, task_label = spec$label, paired = spec$paired, .before = 1)
}) |>
  list_rbind() |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(task, database))


# ---- phenotype association -----------------------------------------------------------------

outcomes <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_fcsa_mixed", "d_nfibre_mixed",
  "d_nfibre_I", "d_mcsa", "d_1rm_legpress", "d_1rm_ext", "volume_load"
)
outcome_of <- function(name, subjects) deframe(select(phenotype, subject, all_of(name)))[subjects]

# The t approximation, not cor.test's default: above nine observations its "exact" Spearman p is
# an Edgeworth series that returns exactly 0 in the far tail. Reading the two fields off the htest
# rather than tidying it is ten times faster.
spearman_by_row <- function(values, outcome) {
  usable <- !is.na(outcome)
  values <- values[, usable, drop = FALSE]
  outcome <- outcome[usable]
  fits <- apply(values, 1, \(row) {
    test <- stats::cor.test(row, outcome, method = "spearman", exact = FALSE)
    c(r = unname(test$estimate), p = test$p.value)
  })
  tibble(
    set_id = rownames(values), n = length(outcome),
    r = unname(fits["r", ]), p = unname(fits["p", ])
  )
}
set_association <- imap(windows, \(values, window) {
  map(set_names(outcomes), \(name) spearman_by_row(values, outcome_of(name, colnames(values)))) |>
    list_rbind(names_to = "outcome") |>
    mutate(window = window, .before = 1)
}) |>
  list_rbind() |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(window, outcome, database))

by_arm <- map(set_names(c("HR", "LR")), function(arm) {
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
  set_auc |>
    summarise(tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05),
      .by = c(task_label, database)
    ) |>
    transmute(analysis = "classification", comparison = task_label, database, tested, nominal,
      fdr_sig
    ),
  set_association |>
    summarise(tested = n(), nominal = sum(p < 0.05), fdr_sig = sum(fdr < 0.05),
      .by = c(window, outcome, database)
    ) |>
    transmute(analysis = paste0("association: ", window), comparison = outcome, database,
      tested, nominal, fdr_sig
    )
) |>
  mutate(expected = 0.05 * tested, ratio = round(nominal / expected, 2), .after = tested)
print(as.data.frame(
  chance_expectation |>
    filter(analysis == "classification") |>
    select(comparison, database, tested, nominal, ratio, fdr_sig)
))


# ---- figures: every nominal result, 16 panels to a page -------------------------------------

# Every set reaching nominal p gets a panel, ordered by collection then p.
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

# One PDF holding every page. A single-page figure also gets a PNG; paged ones do not, since a
# PNG per page adds tens of megabytes to the repository on every re-run.
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

# `draw` builds the plot from one page's rows; a paginating facet would build every panel for
# every page.
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
chance_line <- paste0(
  "Nominal p, uncorrected. Each collection is read against its own expectation: ",
  paste(sprintf(
    "%s %.1f of %d", collection_sizes$database,
    collection_sizes$n * 0.05, collection_sizes$n
  ), collapse = ", "), "."
)

draw_roc_figure <- function(task_name) {
  spec <- tasks[[task_name]]
  hits <- set_auc |>
    filter(task == task_name, p < 0.05) |>
    arrange(database, p)
  if (!nrow(hits)) {
    message("no set reaches nominal p for ", task_name)
    return(0L)
  }
  curves <- pmap(hits, function(set_id, database, pathway, auc, p, fdr, ...) {
    coordinates <- pROC::coords(fit_roc(spec$positive[set_id, ], spec$negative[set_id, ]), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      # An AUC below 0.5 separates the other way; it is coloured, not hidden.
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = paste0(
        database, ": ", enrichVolcano::ev_clean_label(pathway),
        "\nAUC ", sprintf("%.2f", auc), "   p ", signif(p, 2), "   q ", signif(fdr, 2)
      )
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel))) |>
    arrange(panel, fpr, tpr)

  save_paged(curves, \(page) {
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
        x = "1 - specificity", y = "sensitivity", title = spec$label,
        subtitle = sprintf(
          "singscore per set, %d %s, AUC with direction fixed", ncol(spec$positive), spec$unit
        ),
        caption = stringr::str_wrap(width = 190, sprintf(
          paste(
            "All %d sets reaching nominal p, by collection then p. Red separates higher in",
            "%s, blue lower. Shading is area under the curve, dashed line chance.",
            "p from the Wilcoxon %s test, q from BH within collection and task. %s"
          ),
          nrow(hits), spec$favours, if (spec$paired) "signed-rank" else "rank-sum", chance_line
        ))
      )
  }, sprintf("01_roc_%02d_%s", match(task_name, names(tasks)), task_name), strip = 0.6)
}
# The floor gets no figure: chance_expectation reports what the method returns when nothing is
# there.
drawn_tasks <- setdiff(names(tasks), "Baseline_HRvLR")
roc_drawn <- map_int(set_names(drawn_tasks), draw_roc_figure)

window_titles <- c(
  training = "training change (T2 - T1)", baseline = "level at T1",
  acute = "acute change (T3 - T2)"
)
draw_association_figure <- function(which_window) {
  hits <- set_association |>
    filter(window == which_window, p < 0.05) |>
    arrange(database, p)
  if (!nrow(hits)) {
    message("no set reaches nominal p for ", which_window)
    return(0L)
  }
  values <- windows[[which_window]]
  points <- pmap(hits, function(set_id, database, pathway, outcome, n, r, p, fdr, ...) {
    arms <- by_arm |>
      filter(set_id == !!set_id, outcome == !!outcome, window == which_window) |>
      mutate(text = sprintf("%-3s n=%d  r=%+.2f  p=%.3f", arm, n, r, p))
    tibble(
      panel = paste0(
        database, ": ", enrichVolcano::ev_clean_label(pathway), "   vs   ", outcome,
        "\nr = ", sprintf("%+.2f", r), "   p ", signif(p, 2), "   q ", signif(fdr, 2)
      ),
      score = values[set_id, ], outcome_value = outcome_of(outcome, colnames(values)),
      arm = arm_of[colnames(values)], caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    filter(!is.na(outcome_value)) |>
    mutate(panel = factor(panel, levels = unique(panel)))

  save_paged(points, \(page) {
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
        x = paste0("set score, ", window_titles[[which_window]], ", one point per subject"),
        y = "phenotype",
        title = paste("Set score against phenotype,", window_titles[[which_window]]),
        subtitle = "Spearman, subjects pooled across arms",
        caption = stringr::str_wrap(width = 220, sprintf(
          paste(
            "All %d set-outcome pairs reaching nominal p, by collection then p. Lines are fitted",
            "within each arm; the header r is the pooled correlation, q from BH within collection,",
            "window and outcome. %s"
          ),
          nrow(hits), chance_line
        ))
      )
  }, sprintf("02_association_%02d_%s", match(which_window, names(window_titles)), which_window),
  scales = "free", panel = 3.1, header = 2.4
  )
}
association_drawn <- map_int(set_names(names(window_titles)), draw_association_figure)

# Every set reaching nominal p on a drawn task or window has a panel.
stopifnot(
  roc_drawn == map_int(drawn_tasks, \(task_name) sum(set_auc$task == task_name & set_auc$p < 0.05)),
  association_drawn == map_int(names(window_titles), \(w) {
    sum(set_association$window == w & set_association$p < 0.05)
  })
)

# Whether a collection clears chance gets its own panel, not just a sheet.
chance_figure <- chance_expectation |>
  filter(analysis == "classification") |>
  mutate(comparison = factor(comparison, levels = map_chr(tasks, "label")))
save_pages(
  list(ggplot(chance_figure, aes(ratio, database, fill = ratio > 1)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%d of %d", nominal, tested)),
      hjust = -0.12, size = 2.5, colour = "grey25"
    ) +
    facet_wrap(~comparison, ncol = 2) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(
      x = "observed nominal hits / chance expectation", y = NULL,
      title = "Nominal hits relative to chance, by collection",
      subtitle = sprintf("Wilcoxon per set, %d collections, uncorrected p", nrow(collection_sizes)),
      caption = paste(
        "Bar length is observed nominal hits divided by the count that collection returns under",
        "the null. Red clears 1, grey does not. Labels give observed of tested.",
        "Table: c_data/05_classify_and_associate_sets.xlsx, chance_expectation sheet."
      )
    ) +
    figure_theme +
    theme(panel.grid.major.y = element_blank())),
  "03_chance_expectation",
  width = 10, height = 9
)


# ---- one workbook --------------------------------------------------------------------------

packages <- c("here", "pROC", "dplyr", "purrr", "ggplot2", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
set_detail <- catalog |>
  left_join(
    fgsea_results |>
      filter(contrast == "Training_Interaction") |>
      select(set_id, interaction_nes = NES, interaction_p = p),
    by = "set_id"
  ) |>
  arrange(database, interaction_p)

sheets <- list(
  chance_expectation = chance_expectation,
  set_auc = arrange(set_auc, p),
  set_association = arrange(set_association, p),
  set_by_arm = filter(by_arm, p < 0.05),
  set_catalog = set_detail,
  set_scores = rownames_to_column(as.data.frame(set_score), "set_id"),
  input_manifest = manifest,
  package_versions = versions
)
descriptions <- c(
  chance_expectation = "Nominal hits against chance, per collection. Read this first.",
  set_auc = "How well each set separates each task. AUC from ranks, p from the Wilcoxon test.",
  set_association = "Set score against each phenotype, per window, subjects pooled.",
  set_by_arm = "Within-arm correlations with p < 0.05. Every row is in set_results.rds.",
  set_catalog = "Every tested set with its collection and Training_Interaction NES.",
  set_scores = "The set by sample score matrix the analyses above were computed on.",
  input_manifest = "Which files were read and their checksums.",
  package_versions = "Package versions at the time of the run."
)
read_me <- tibble(
  sheet = names(sheets),
  rows = map_int(sheets, nrow),
  holds = unname(descriptions[names(sheets)])
)
writexl::write_xlsx(
  c(list(read_me = read_me), sheets),
  file.path(out, "05_classify_and_associate_sets.xlsx")
)
saveRDS(
  list(
    catalog = catalog, set_auc = set_auc, set_association = set_association,
    by_arm = by_arm, chance_expectation = chance_expectation,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "set_results.rds"),
  compress = "xz"
)
combined <- file.path(figure_dir, "05_classify_and_associate_sets_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote 05_classify_and_associate_sets.xlsx (", length(sheets) + 1, " sheets) and ",
  length(pages), " figures bundled into ", qpdf::pdf_length(combined), " pages"
)
