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
# Clear last run's figures so the bundle holds only this run's pages.
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
    feature = rownames(spec$positive), n_positive = ncol(spec$positive),
    n_negative = ncol(spec$negative)
  ) |>
    mutate(
      auc = map_dbl(seq_along(feature), \(i) {
        as.numeric(pROC::auc(fit_roc(spec$positive[i, ], spec$negative[i, ])))
      }),
      p_value = map_dbl(seq_along(feature), \(i) {
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
  rename(set_id = feature) |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p_value, "BH"), .by = c(task, database))

auc_counts <- set_auc |>
  summarise(
    features = n(), nominal = sum(p_value < 0.05), expected = round(n() * 0.05, 1),
    fdr_sig = sum(fdr < 0.05), max_auc = round(max(auc), 2),
    .by = c(task, task_label, database, n_positive)
  ) |>
  mutate(analysis = "classification", .before = 1)
print(as.data.frame(
  auc_counts |>
    mutate(ratio = round(nominal / expected, 2)) |>
    select(task, database, features, nominal, expected, ratio, fdr_sig) |>
    arrange(task, database)
))


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
    feature = rownames(values), n = length(outcome),
    r = unname(fits["r", ]), p = unname(fits["p", ])
  )
}
set_association <- imap(windows, \(values, window) {
  map(set_names(outcomes), \(name) spearman_by_row(values, outcome_of(name, colnames(values)))) |>
    list_rbind(names_to = "outcome") |>
    mutate(window = window, .before = 1)
}) |>
  list_rbind() |>
  rename(set_id = feature) |>
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
  list_rbind(names_to = "arm") |>
  rename(set_id = feature)

chance_expectation <- bind_rows(
  select(auc_counts, analysis,
    comparison = task_label, database, n = n_positive,
    features, nominal, expected, fdr_sig
  ),
  set_association |>
    summarise(
      features = n(), nominal = sum(p < 0.05), expected = round(n() * 0.05, 1),
      fdr_sig = sum(fdr < 0.05), .by = c(window, outcome, database, n)
    ) |>
    transmute(
      analysis = paste0("association: ", window), comparison = outcome, database, n,
      features, nominal, expected, fdr_sig
    )
) |>
  mutate(ratio = round(nominal / expected, 2)) |>
  relocate(ratio, .after = expected)


# ---- figures -------------------------------------------------------------------------------

# A task can reach nominal p in hundreds of sets, so the ROC and association figures show the
# strongest two from each collection. The hit matrices after them draw every nominal set.
per_database <- 2

save_faceted <- function(plot, name, n_panels, columns, panel = 2.6, header = 2.1,
                         width = 1.2 + columns * panel,
                         height = header + ceiling(n_panels / columns) * panel) {
  # Captions run to several sentences; wrap them to the page width.
  plot$labels$caption <- stringr::str_wrap(plot$labels$caption, width = floor(width * 17))
  figure <- plot +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 7.6, lineheight = 1.3, margin = margin(3, 3, 5, 3)),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40"),
      legend.position = "top"
    )
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 220, bg = "white"
    )
  })
  message("drew ", name, ": ", n_panels, " panels")
  n_panels
}
chance_line <- paste0(
  "Nominal p, uncorrected. Each collection is read against its own expectation: ",
  paste(sprintf(
    "%s %.1f of %d", collection_sizes$database,
    collection_sizes$n * 0.05, collection_sizes$n
  ), collapse = ", "), "."
)

draw_roc_figure <- function(task_name, index) {
  spec <- tasks[[task_name]]
  hits <- set_auc |>
    filter(task == task_name, p_value < 0.05) |>
    slice_min(p_value, n = per_database, by = database, with_ties = FALSE) |>
    arrange(database, p_value)
  if (!nrow(hits)) {
    message("no set reaches nominal p for ", task_name)
    return(0L)
  }
  curves <- pmap(hits, function(set_id, database, pathway, auc, p_value, ...) {
    coordinates <- pROC::coords(fit_roc(spec$positive[set_id, ], spec$negative[set_id, ]), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = paste0(
        database, ": ", enrichVolcano::ev_clean_label(pathway),
        "\nAUC ", sprintf("%.2f", auc), "    p = ", signif(p_value, 2)
      )
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel))) |>
    arrange(panel, fpr, tpr)

  columns <- min(4, nrow(hits))
  save_faceted(
    ggplot(curves, aes(fpr, tpr, colour = direction, fill = direction)) +
      geom_abline(linetype = "22", linewidth = 0.35, colour = "grey60") +
      geom_ribbon(aes(ymin = 0, ymax = tpr), alpha = 0.16, colour = NA) +
      geom_step(linewidth = 0.8, direction = "hv") +
      scale_colour_manual(values = c(higher = "#B2182B", lower = "#2166AC"), guide = "none") +
      scale_fill_manual(values = c(higher = "#B2182B", lower = "#2166AC"), guide = "none") +
      facet_wrap(~panel, ncol = columns) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      labs(
        x = "1 - specificity", y = "sensitivity", title = spec$label,
        subtitle = sprintf(
          "singscore per set, %d %s, AUC with direction fixed", ncol(spec$positive), spec$unit
        ),
        caption = sprintf(
          paste(
            "Strongest %d sets per collection of %d reaching nominal p. Red separates higher in",
            "%s, blue lower. Shading is area under the curve, dashed line chance.",
            "p from the Wilcoxon %s test. %s"
          ),
          per_database, sum(set_auc$task == task_name & set_auc$p_value < 0.05),
          spec$favours, if (spec$paired) "signed-rank" else "rank-sum", chance_line
        )
      ),
    sprintf("01_roc_%02d_%s", index, task_name), nrow(hits), columns
  )
}
# The floor gets no figure. It shows what the method returns when nothing is there, so
# chance_expectation reports it as a number.
drawn_tasks <- setdiff(names(tasks), "Baseline_HRvLR")
roc_drawn <- imap_int(set_names(drawn_tasks), \(task_name, i) {
  draw_roc_figure(task_name, match(task_name, drawn_tasks))
})

draw_association_figure <- function(which_window, index, title) {
  hits <- set_association |>
    filter(window == which_window, p < 0.05) |>
    slice_min(p, n = per_database, by = database, with_ties = FALSE) |>
    arrange(database, p)
  if (!nrow(hits)) {
    message("no set reaches nominal p for ", which_window)
    return(0L)
  }
  values <- windows[[which_window]]
  points <- pmap(hits, function(set_id, database, pathway, outcome, n, r, p, ...) {
    arms <- by_arm |>
      filter(set_id == !!set_id, outcome == !!outcome, window == which_window) |>
      mutate(text = sprintf("%-3s n=%d  r=%+.2f  p=%.3f", arm, n, r, p))
    tibble(
      panel = paste0(
        database, ": ", enrichVolcano::ev_clean_label(pathway), "   vs   ", outcome,
        "\nr = ", sprintf("%+.2f", r), "    p = ", signif(p, 2)
      ),
      d_score = values[set_id, ], d_outcome = outcome_of(outcome, colnames(values)),
      arm = arm_of[colnames(values)], caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    # A subject with no value for that phenotype has nothing to plot.
    filter(!is.na(d_outcome)) |>
    mutate(panel = factor(panel, levels = unique(panel)))

  columns <- min(4, nrow(hits))
  save_faceted(
    ggplot(points, aes(d_score, d_outcome, colour = arm, fill = arm)) +
      geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
      geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.12, linewidth = 0.6) +
      geom_point(size = 1.9, alpha = 0.9) +
      geom_text(
        data = distinct(points, panel, caption), inherit.aes = FALSE,
        aes(x = -Inf, y = Inf, label = caption), family = "mono",
        hjust = -0.05, vjust = 1.25, size = 2.3, lineheight = 1.25, colour = "grey25"
      ) +
      facet_wrap(~panel, scales = "free", ncol = columns) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.3))) +
      scale_colour_manual(values = c(HR = "#2166AC", LR = "#B2182B"), name = NULL) +
      scale_fill_manual(values = c(HR = "#2166AC", LR = "#B2182B"), name = NULL) +
      labs(
        x = paste0("set score, ", title, ", one point per subject"),
        y = "phenotype", title = paste("Set score against phenotype,", title),
        subtitle = "Spearman, subjects pooled across arms",
        caption = sprintf(
          paste(
            "Strongest %d set-outcome pairs per collection of %d reaching nominal p.",
            "Lines are fitted within each arm; the header r is the pooled correlation that",
            "selected the panel. %s"
          ),
          per_database, sum(set_association$window == which_window & set_association$p < 0.05),
          chance_line
        )
      ),
    sprintf("02_association_%02d_%s", index, which_window), nrow(hits), columns,
    panel = 3.1, header = 2.4
  )
}
window_titles <- c(
  training = "training change (T2 - T1)", baseline = "level at T1",
  acute = "acute change (T3 - T2)"
)
iwalk(window_titles, \(title, window) {
  draw_association_figure(window, match(window, names(window_titles)), title)
})

# Every collection with a nominal hit is on each figure, with up to per_database sets.
stopifnot(all(map_lgl(drawn_tasks, \(task_name) {
  available <- filter(set_auc, task == task_name, p_value < 0.05)
  roc_drawn[[task_name]] == sum(pmin(table(available$database), per_database))
})))

chance_figure <- chance_expectation |>
  filter(analysis == "classification") |>
  mutate(comparison = factor(comparison, levels = map_chr(tasks, "label")))
save_faceted(
  ggplot(chance_figure, aes(ratio, database, fill = ratio > 1)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%d of %d", nominal, features)),
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
    theme(panel.grid.major.y = element_blank()),
  "03_chance_by_database", n_distinct(chance_figure$comparison), 2,
  width = 10, height = 9
)

# Every set nominal in at least one column, as a dot matrix, 75 rows to a page. The tasks get one
# matrix, each association window another.
hit_pages <- function(data, columns, title, subtitle, fill_label, midpoint, limits, name) {
  ranked <- data |>
    summarise(n_nominal = sum(p < 0.05), best = min(p), .by = label) |>
    filter(n_nominal > 0) |>
    arrange(desc(n_nominal), best)
  page_of <- ceiling(seq_len(nrow(ranked)) / 75)
  for (page in seq_len(max(page_of, 0))) {
    rows <- ranked$label[page_of == page]
    shown <- data |>
      filter(label %in% rows) |>
      mutate(
        label = factor(label, levels = rev(rows)),
        column = factor(column, levels = columns), nominal = p < 0.05
      )
    figure <- ggplot(shown, aes(column, label)) +
      geom_point(data = filter(shown, !nominal), colour = "grey85", size = 0.6) +
      geom_point(
        data = filter(shown, nominal), aes(size = -log10(p), fill = effect),
        shape = 21, colour = "grey30", stroke = 0.2
      ) +
      geom_point(
        data = filter(shown, fdr < 0.05), aes(size = -log10(p)),
        shape = 21, colour = "black", stroke = 0.9, fill = NA
      ) +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B", midpoint = midpoint, limits = limits
      ) +
      scale_size_continuous(range = c(1, 4.5), limits = c(-log10(0.05), NA)) +
      scale_x_discrete(drop = FALSE) +
      labs(
        title = sprintf("%s (%d/%d)", title, page, max(page_of)),
        subtitle = sprintf(
          "%s; %d sets nominal in at least one column; rows %d-%d", subtitle, nrow(ranked),
          min(which(page_of == page)), max(which(page_of == page))
        ),
        x = NULL, y = NULL, fill = fill_label, size = "-log10 p",
        caption = paste0(
          "Rows: sets at nominal p < 0.05 in at least one column, most columns first, then by ",
          "best p. Filled dot: nominal, fill = ", fill_label, ", size = -log10 p. Black ring: BH ",
          "< 0.05 within collection. Grey speck: not nominal. ",
          "Table: c_data/05_classify_and_associate_sets.xlsx."
        )
      ) +
      theme_minimal(base_size = 9) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1),
        axis.text.y = element_text(size = if (length(rows) > 30) 5.5 else 8),
        plot.caption = element_text(size = 7, colour = "grey45", hjust = 0)
      )
    walk(c("png", "pdf"), \(extension) {
      ggsave(file.path(figure_dir, sprintf("%s_%02d.%s", name, page, extension)), figure,
        width = 11, height = 8.5, dpi = 200, bg = "white"
      )
    })
  }
  max(page_of, 0)
}
set_label <- catalog |>
  transmute(set_id, label = paste0("[", database, "] ", enrichVolcano::ev_clean_label(pathway))) |>
  mutate(label = gsub("\n", " ", label), label = if_else(duplicated(label), set_id, label)) |>
  deframe()
hit_counts <- c(
  classification = hit_pages(
    transmute(set_auc,
      label = set_label[set_id], column = task, effect = auc,
      p = p_value, fdr
    ),
    names(tasks), "Set classification hits", "singscore AUC, Wilcoxon p per task",
    "AUC", 0.5, c(0, 1), "04_hits_classification"
  ),
  imap_int(window_titles, \(title, window) {
    hit_pages(
      set_association |>
        filter(window == .env$window) |>
        transmute(label = set_label[set_id], column = outcome, effect = r, p, fdr),
      outcomes, paste("Set association hits,", title), "Spearman against each phenotype",
      "rho", 0, c(-1, 1), sprintf("05_hits_%s", window)
    )
  })
)
print(hit_counts)


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
  set_auc = arrange(set_auc, p_value),
  set_association = arrange(set_association, p),
  set_by_arm = by_arm,
  set_catalog = set_detail,
  set_scores = rownames_to_column(as.data.frame(set_score), "set_id"),
  input_manifest = manifest,
  package_versions = versions
)
descriptions <- c(
  chance_expectation = "Nominal hits against chance, per collection. Read this first.",
  set_auc = "How well each set separates each task. AUC from ranks, p from the Wilcoxon test.",
  set_association = "Set score against each phenotype, per window, subjects pooled.",
  set_by_arm = "The same correlation computed inside HR and inside LR, descriptive.",
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
  "wrote 05_classify_and_associate_sets.xlsx (", length(sheets) + 1, " sheets) and a ",
  length(pages), "-page figure PDF"
)
