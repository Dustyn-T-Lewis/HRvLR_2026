# Ties each module eigengene to phenotype, the protein level's two ways.
#
# Sample model: every T1 and T2 biopsy carries its own trait values. Per trait, limma fits
# abundance ~ timepoint + between + within with subject blocked, where between is the subject's
# mean trait and within the biopsy's deviation from it. within asks whether a protein rises in the
# biopsies where that person's trait rose; between asks whether people with more of the trait
# carry more of the protein. Arm stays out: it was defined from one of these traits.
#
# Change score: one value per subject and window (training change, T1 level, acute change),
# Spearman against 8 outcomes. Pooled p is the t approximation that cor.test(exact = FALSE) uses;
# within-arm p is exact, since 4 to 8 subjects make the approximation return p = 0 at |rho| = 1.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
  library(grid)
  library(ComplexHeatmap)
})
source(here("R", "helpers.R"))

stage <- here("04_Network", "04_Associate")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  modules = "04_Network/01_Build_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
modules <- readRDS(paths[["modules"]])
values <- modules$eigengenes
meta <- as_tibble(modules$meta)
stopifnot(identical(colnames(values), meta$sample_id))
feature_label <- set_names(rownames(values))
feature_noun <- "modules"
value_label <- "eigengene"
row_groups <- list(modules = modules$module_summary$module)

# One row per subject and biopsy, pre as T1 and post as T2. Size and strength traits only; the
# MyoVision fibre counts stay out.
traits <- c("fcsa_mixed", "fcsa_type1", "fcsa_type2", "mcsa", "leg_press_1rm", "leg_ext_1rm")
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
biopsy <- phenotype |>
  select(subject, arm, matches("_(pre|post)_")) |>
  pivot_longer(-c(subject, arm),
    names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*"
  ) |>
  pivot_wider(names_from = trait) |>
  mutate(timepoint = if_else(when == "pre", "T1", "T2")) |>
  select(subject, arm, timepoint, all_of(traits))
per_subject <- biopsy |>
  pivot_longer(all_of(traits), names_to = "trait") |>
  pivot_wider(names_from = timepoint) |>
  transmute(subject, trait, pct = 100 * (T2 - T1) / T1) |>
  pivot_wider(names_from = trait, values_from = pct, names_prefix = "pct_") |>
  left_join(select(phenotype, subject, comp_hypertrophy, volume_load_total_kg), by = "subject")
outcomes <- c("comp_hypertrophy", "volume_load_total_kg", paste0("pct_", traits))
outcome_of <- \(name, subjects) set_names(per_subject[[name]], per_subject$subject)[subjects]
arm_of <- set_names(phenotype$arm, phenotype$subject)

samples <- meta |>
  filter(timepoint %in% c("T1", "T2")) |>
  select(sample_id, subject, timepoint) |>
  left_join(select(biopsy, -arm), by = c("subject", "timepoint"))
model <- map(set_names(traits), \(trait) {
  s <- samples |>
    filter(!is.na(.data[[trait]])) |>
    mutate(
      between = ave(.data[[trait]], subject),
      within = .data[[trait]] - between
    )
  design <- model.matrix(~ timepoint + between + within, s)
  y <- values[, s$sample_id]
  correlation <- duplicateCorrelation(y, design, block = s$subject)$consensus.correlation
  fit <- lmFit(y, design, block = s$subject, correlation = correlation) |>
    eBayes(robust = TRUE)
  map(c("within", "between"), \(term) {
    topTable(fit, coef = term, number = Inf, sort.by = "none") |>
      as_tibble(rownames = "feature") |>
      transmute(
        trait, term, feature,
        n_samples = nrow(s), effect = logFC, t, p = P.Value, fdr = adj.P.Val,
        qvalue = storey_q(P.Value)
      )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(label = feature_label[feature], .after = feature)

windows <- list(
  training = subject_change(values, meta, "T1", "T2"),
  baseline = subject_values(values, meta, "T1"),
  acute = subject_change(values, meta, "T2", "T3")
)

correlate <- function(x, exact = FALSE) {
  map(set_names(outcomes), \(name) {
    outcome <- outcome_of(name, colnames(x))
    result <- spearman_by_row(x, outcome)
    if (exact) {
      tested <- which(!is.na(result$rho))
      result$p[tested] <- map_dbl(tested, \(i) exact_spearman_p(x[i, ], outcome))
    }
    result
  }) |>
    list_rbind(names_to = "outcome")
}
change_score <- imap(windows, \(x, window) mutate(correlate(x), window = window, .before = 1)) |>
  list_rbind() |>
  filter(!is.na(rho)) |>
  mutate(fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = c(window, outcome)) |>
  mutate(label = feature_label[feature], .after = feature)
by_arm <- map(set_names(c("HR", "LR")), \(arm) {
  imap(windows, \(x, window) {
    kept <- x[, arm_of[colnames(x)] == arm, drop = FALSE]
    mutate(correlate(kept, exact = TRUE), window = window, .before = 1)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "arm") |>
  filter(!is.na(rho))

# comp_hypertrophy defines the arms, so a slope difference on it re-measures the split; it stays
# in the pooled correlation and the two-predictor model. Each arm's slope needs five subjects
# with both the feature and the outcome observed.
arm_difference <- imap(windows, \(x, window) {
  map(set_names(setdiff(outcomes, "comp_hypertrophy")), \(name) {
    predictors <- data.frame(
      arm = factor(arm_of[colnames(x)], c("LR", "HR")), outcome = outcome_of(name, colnames(x))
    )
    seen <- !is.na(x) & rep(!is.na(predictors$outcome), each = nrow(x))
    per_arm <- sapply(c("HR", "LR"), \(a) rowSums(seen[, predictors$arm %in% a, drop = FALSE]))
    subject_fit(
      x[apply(per_arm, 1, min) >= 5, , drop = FALSE], predictors, ~ arm * outcome,
      "armHR:outcome"
    )
  }) |>
    list_rbind(names_to = "outcome") |>
    mutate(window = window, .before = 1)
}) |>
  list_rbind() |>
  select(-coefficient) |>
  mutate(label = feature_label[feature], .after = feature) |>
  mutate(fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = c(window, outcome))
strength_tests <- c("pct_leg_press_1rm", "pct_leg_ext_1rm")
two_predictor <- imap(windows, \(x, window) {
  map(set_names(strength_tests), \(test) {
    predictors <- data.frame(
      growth = outcome_of("comp_hypertrophy", colnames(x)), strength = outcome_of(test, colnames(x))
    )
    subject_fit(x, predictors, ~ growth + strength, c("growth", "strength"))
  }) |>
    list_rbind(names_to = "strength_test") |>
    mutate(window = window, .before = 1)
}) |>
  list_rbind() |>
  mutate(
    outcome = if_else(coefficient == "growth", "comp_hypertrophy", strength_test),
    .after = coefficient
  ) |>
  mutate(label = feature_label[feature], .after = feature) |>
  mutate(fdr = p.adjust(p, "BH"), qvalue = storey_q(p), .by = c(window, strength_test, coefficient))

counts <- bind_rows(
  summarise(model,
    n_tested = sum(!is.na(p)), n_nominal = sum(p < 0.05, na.rm = TRUE),
    n_fdr_05 = sum(fdr < 0.05, na.rm = TRUE), .by = c(term, trait)
  ) |>
    transmute(analysis = paste("model:", term), comparison = trait, n_tested, n_nominal, n_fdr_05),
  summarise(change_score,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(window, outcome)
  ) |>
    transmute(
      analysis = paste("change score:", window), comparison = outcome, n_tested, n_nominal,
      n_fdr_05
    ),
  summarise(arm_difference,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(window, outcome)
  ) |>
    transmute(
      analysis = paste("arm difference:", window), comparison = outcome,
      n_tested, n_nominal, n_fdr_05
    ),
  summarise(two_predictor,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(window, strength_test, coefficient)
  ) |>
    transmute(
      analysis = paste("two predictor:", window),
      comparison = paste(coefficient, strength_test, sep = ", "),
      n_tested, n_nominal, n_fdr_05
    )
)



# Units come from the phenotype column names.
unit_cols <- grep("_pre_", names(phenotype), value = TRUE)
unit_of <- set_names(sub(".*_pre_", "", unit_cols), sub("_pre_.*", "", unit_cols))
outcome_unit <- c(
  comp_hypertrophy = "score", volume_load_total_kg = "kg",
  set_names(rep("%", length(traits)), paste0("pct_", traits))
)
window_axis <- c(training = "change, T1 to T2", baseline = "T1 level", acute = "change, T2 to T3")
arm_colours <- c(HR = "#2166AC", LR = "#B2182B")

# Overviews: sample-model t by trait, then change-score rho by outcome, one page per window.
t_limit <- cap(model$t)
row_order <- paste(
  "Rows follow module size."
)
overviews <- imap(row_groups, \(rows, group) {
  term_grid <- \(term, value) {
    as_grid(filter(model, term == .env$term), rows, traits, "trait", value)
  }
  model_maps <- heat(
    term_grid("within", \(d) d$t), "Within-person t", c(-t_limit, 0, t_limit),
    "t"
  ) +
    heat(term_grid("between", \(d) d$t), "Between-person t", c(-t_limit, 0, t_limit), NULL) +
    heat(term_grid("within", \(d) status_of(d$p, d$fdr)), "Within: p and BH q") +
    heat(term_grid("between", \(d) status_of(d$p, d$fdr)), "Between: p and BH q",
      row_names = TRUE
    )
  model_page <- heatmap_page(model_maps, sprintf("Sample model overview: %s", group), paste(
    sprintf("All %s %s, one row each.", format(length(rows), big.mark = ","), feature_noun),
    "Moderated t of the within and between terms, abundance ~ timepoint + between + within,",
    sprintf("subject blocked. Colours stop at the 99th percentile of |t|, %.2g.", t_limit),
    row_order, "Grey: not tested. Table: c_data/04_associate.xlsx, model."
  ))
  window_pages <- map(names(windows), \(w) {
    d <- filter(change_score, window == w)
    maps <- heat(
      as_grid(d, rows, outcomes, "outcome", \(d) d$rho),
      "Spearman rho, both arms pooled", c(-1, 0, 1), "rho"
    ) +
      heat(as_grid(d, rows, outcomes, "outcome", \(d) status_of(d$p, d$fdr)), "p and BH q",
        row_names = TRUE
      )
    heatmap_page(maps, sprintf("Change score overview, %s window: %s", w, group), paste(
      sprintf("All %s %s, one row each.", format(length(rows), big.mark = ","), feature_noun),
      sprintf(
        "Spearman rho between the %s (%s) and each outcome.", value_label, window_axis[[w]]
      ),
      row_order, "Grey: not tested. Table: c_data/04_associate.xlsx, change_score."
    ))
  })
  c(list(model_page), window_pages)
}) |>
  list_flatten()

# Grids: every p < 0.05 hit, strongest first (p, then |t| or |rho|, then feature ID).
trait_samples <- map(set_names(traits), \(trait) filter(samples, !is.na(.data[[trait]])))
seen_in <- \(feature, trait) !is.na(values[feature, trait_samples[[trait]]$sample_id])
model_hits <- model |>
  filter(p < 0.05) |>
  mutate(
    n_biopsies = map2_int(feature, trait, \(f, tr) sum(seen_in(f, tr))),
    n_subjects = map2_int(feature, trait, \(f, tr) {
      n_distinct(trait_samples[[tr]]$subject[seen_in(f, tr)])
    }),
    name = label,
    panel = sprintf(
      "%s\n%s %.3g  p %.2g  q %.2g\n%d biopsies, %d subjects", name, term, effect, p, fdr,
      n_biopsies, n_subjects
    )
  ) |>
  arrange(match(trait, traits), desc(term), p, desc(abs(t)), feature) |>
  paginate(c("trait", "term"))
model_page <- function(hits) {
  trait <- hits$trait[1]
  term <- hits$term[1]
  info <- counts |>
    filter(analysis == paste("model:", term), comparison == trait) |>
    summarise(across(c(n_tested, n_nominal, n_fdr_05), sum))
  points <- hits |>
    select(feature, panel) |>
    cross_join(select(samples, sample_id, subject, timepoint, trait_value = all_of(trait))) |>
    mutate(
      value = values[cbind(feature, sample_id)], arm = arm_of[subject],
      panel = factor(panel, slots(hits$panel))
    ) |>
    filter(!is.na(value), !is.na(trait_value))
  ggplot(points, aes(trait_value, value, colour = arm)) +
    geom_line(
      aes(group = subject),
      data = \(d) filter(d, n() == 2, .by = c(panel, subject)), alpha = 0.6
    ) +
    geom_point(aes(shape = timepoint), size = 1.5) +
    facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free", drop = FALSE) +
    scale_colour_manual(values = arm_colours, name = NULL) +
    scale_shape_manual(values = c(T1 = 1, T2 = 16), name = NULL) +
    labs(
      title = sprintf("%s, %s-person term", trait, term),
      subtitle = counts_line(info, feature_noun, hits),
      x = sprintf("%s at the biopsy (%s)", trait, unit_of[[trait]]), y = value_label,
      caption = stringr::str_wrap(width = 160, paste(
        "limma, abundance ~ timepoint + between + within, subject blocked, on the T1 and T2",
        "biopsies with the trait measured. Lines join each subject's two biopsies: the within",
        "term reads their slopes, the between term the spread of subjects. Header: effect per",
        "trait unit, p, BH q within trait and term,",
        "and the biopsies and subjects where the feature was observed.",
        "Table: c_data/04_associate.xlsx, model."
      ))
    ) +
    page_theme
}

arm_text <- by_arm |>
  mutate(text = sprintf("%-3s n=%d rho=%+.2f p=%.3f", arm, n_subjects, rho, p)) |>
  summarise(caption = paste(text, collapse = "\n"), .by = c(window, outcome, feature))
change_hits <- change_score |>
  filter(p < 0.05) |>
  mutate(
    name = label,
    panel = sprintf(
      "%s vs %s (%s)\nrho %+.2f  p %.2g  q %.2g  n %d", name, outcome, outcome_unit[outcome],
      rho, p, fdr, n_subjects
    )
  ) |>
  arrange(match(window, names(windows)), p, desc(abs(rho)), feature, outcome) |>
  paginate("window")
page_captions <- c(
  "change score" = paste(
    "One point per subject. Lines are a least-squares fit per arm, drawn as a guide; the",
    "tests are Spearman. Header: rho and p pooled over both arms,",
    "BH q within window, outcome%s, and the subjects used.",
    "Corner: each arm's rho and exact p. Arm is the median split of",
    "comp_hypertrophy, so pooled correlations with it and with fibre area change partly",
    "carry the group split. Table: c_data/04_associate.xlsx, change_score and",
    "change_score_by_arm."
  ),
  "arm difference" = paste(
    "One point per subject. Lines are a least-squares fit per arm. Header: HR's slope minus",
    "LR's per outcome unit from limma, feature ~ arm * outcome, with p, BH q within window,",
    "outcome%s, and the subjects used. Corner: each arm's Spearman rho and exact p.",
    "Table: c_data/04_associate.xlsx, arm_difference."
  ),
  "two predictor" = paste(
    "One point per subject. Header: the growth or strength coefficient from limma, feature ~",
    "comp_hypertrophy + strength gain, each adjusted for the other, with p, BH q within",
    "window, strength test, term%s, and the subjects used. The x axis shows the term's own",
    "outcome; the drawing does not hold the other predictor fixed. Corner: each arm's",
    "Spearman rho and exact p. Table: c_data/04_associate.xlsx, two_predictor."
  )
)
page_captions[] <- sprintf(page_captions, "")
change_page <- function(hits, kind = "change score") {
  w <- hits$window[1]
  info <- counts |>
    filter(analysis == paste0(kind, ": ", w)) |>
    summarise(across(c(n_tested, n_nominal, n_fdr_05), sum))
  points <- hits |>
    left_join(arm_text, by = c("window", "outcome", "feature")) |>
    select(feature, outcome, panel, caption) |>
    cross_join(tibble(subject = colnames(windows[[w]]))) |>
    mutate(
      value = windows[[w]][cbind(feature, subject)],
      outcome_value = map2_dbl(outcome, subject, outcome_of), arm = arm_of[subject],
      panel = factor(panel, slots(hits$panel))
    ) |>
    filter(!is.na(value), !is.na(outcome_value))
  ggplot(points, aes(value, outcome_value, colour = arm)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
    geom_point(size = 1.3) +
    geom_text(
      data = distinct(points, panel, caption), aes(-Inf, Inf, label = caption),
      inherit.aes = FALSE, hjust = -0.03, vjust = 1.1, size = 1.9, family = "mono"
    ) +
    facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free", drop = FALSE) +
    scale_colour_manual(values = arm_colours, name = NULL) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
    labs(
      title = sprintf("%s, %s window", stringr::str_to_sentence(kind), w),
      subtitle = counts_line(info, "feature-outcome pairs", hits),
      x = sprintf("%s, %s", value_label, window_axis[[w]]), y = "outcome (unit in header)",
      caption = stringr::str_wrap(page_captions[[kind]], 160)
    ) +
    page_theme
}
model_pages <- map(split(model_hits, model_hits$pdf_page), model_page)
change_pages <- map(split(change_hits, change_hits$pdf_page), change_page)
slope_hits <- bind_rows(
  arm_difference |>
    filter(p < 0.05) |>
    mutate(
      kind = "arm difference",
      header = sprintf("slope HR-LR %+.3g", effect)
    ),
  two_predictor |>
    filter(p < 0.05) |>
    mutate(
      kind = "two predictor",
      header = sprintf(
        "%s term, %s model %+.3g", coefficient, sub("pct_", "", strength_test), effect
      )
    )
) |>
  mutate(
    name = label,
    panel = sprintf(
      "%s vs %s (%s)\n%s  p %.2g  q %.2g  n %d", name, outcome, outcome_unit[outcome], header,
      p, fdr, n_subjects
    )
  ) |>
  arrange(
    match(kind, c("arm difference", "two predictor")), match(window, names(windows)), p,
    desc(abs(t)), feature, outcome
  ) |>
  paginate(c("kind", "window"))
slope_pages <- map(split(slope_hits, slope_hits$pdf_page), \(hits) change_page(hits, hits$kind[1]))

plot_index <- bind_rows(
  tibble(
    analysis = "overview", label = map_chr(overviews, "title"), pdf_page = seq_along(overviews)
  ),
  transmute(model_hits,
    analysis = "model", trait, term, feature, label,
    pdf_page = pdf_page + length(overviews), position
  ),
  transmute(change_hits,
    analysis = "change score", window, outcome, feature, label,
    pdf_page = pdf_page + length(overviews) + length(model_pages), position
  ),
  transmute(slope_hits,
    analysis = kind, window, outcome, strength_test, coefficient, feature, label,
    pdf_page = pdf_page + length(overviews) + length(model_pages) + length(change_pages), position
  )
)
stopifnot(
  nrow(model_hits) == sum(model$p < 0.05, na.rm = TRUE),
  nrow(change_hits) == sum(change_score$p < 0.05),
  !anyDuplicated(select(filter(plot_index, !is.na(feature)), -label, -pdf_page, -position)),
  nrow(slope_hits) == sum(arm_difference$p < 0.05) + sum(two_predictor$p < 0.05),
  max(plot_index$pdf_page) ==
    length(overviews) + length(model_pages) + length(change_pages) + length(slope_pages)
)

pdf(file.path(figures, "04_associate_figures.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(c(model_pages, change_pages, slope_pages), print)
invisible(dev.off())

write_workbook(
  list(
    counts = counts, model = arrange(model, trait, term, p),
    change_score = arrange(change_score, window, outcome, p),
    change_score_by_arm = by_arm,
    arm_difference = arrange(arm_difference, window, outcome, p),
    two_predictor = arrange(two_predictor, window, strength_test, coefficient, p),
    plot_index = plot_index
  ),
  c(
    "Tested, p < 0.05 and BH < 0.05 counts per term and trait, and per window and outcome.",
    "Sample model: effect per trait unit, t, p, BH and Storey q per module, trait and term.",
    "Change score: pooled Spearman rho, p, BH and Storey q per module, window and outcome.",
    "Change score within each arm, exact p, every module.",
    "Arm difference: HR minus LR slope, t, p, BH and Storey q per module, window and outcome.",
    "Two predictor: growth and strength terms, t, p, BH and Storey q per module, window and test.",
    "PDF page of each overview, and page and panel position of every p < 0.05 row."
  ),
  file.path(out, "04_associate.xlsx"),
  inputs, paths
)
