# Ties each protein to phenotype two ways.
#
# Sample model: every T1 and T2 biopsy carries its own trait values. Per trait, limma fits
# abundance ~ timepoint + between + within with subject blocked, where between is the subject's
# mean trait and within the biopsy's deviation from it. within asks whether a protein rises in the
# biopsies where that person's trait rose; between asks whether people with more of the trait
# carry more of the protein. Arm stays out: it was defined from one of these traits. A protein
# enters only when two thirds of the trait's biopsies saw it, the change score's floor.
#
# Change score: one value per subject and window (training change, T1 level, acute change),
# Spearman against 20 outcomes. Pooled p is the t approximation that cor.test(exact = FALSE) uses;
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

stage <- here("02_Differential_Expression", "04_Associate")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  proteins = "01_Preprocess/02_Normalization/c_data/DAList_normalized.rds",
  imputed = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
dal <- readRDS(paths[["proteins"]])
values <- as.matrix(dal$data)
meta <- as_tibble(dal$metadata)
feature_label <- dal$annotation |>
  mutate(label = if_else(
    is.na(gene) | duplicated(gene) | duplicated(gene, fromLast = TRUE),
    paste0(coalesce(gene, "?"), " (", uniprot_id, ")"), gene
  )) |>
  with(set_names(label, uniprot_id))
feature_noun <- "protein"
value_label <- "log2 abundance"
# Display order for the overviews: average-linkage clustering on 1 - Pearson r of the imputed
# matrix, shared by every protein overview.
imputed <- as.matrix(readRDS(paths[["imputed"]])$data)
row_groups <- list(
  proteins = rownames(imputed)[hclust(as.dist(1 - cor(t(imputed))), "average")$order]
)

# One row per subject and biopsy, pre as T1 and post as T2. type1_share is the MyoVision type I
# fibre count as a percentage of the mixed count.
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
biopsy <- phenotype |>
  select(subject, arm, matches("_(pre|post)_")) |>
  pivot_longer(-c(subject, arm),
    names_to = c("trait", "when"), names_pattern = "(.*)_(pre|post)_.*"
  ) |>
  pivot_wider(names_from = trait) |>
  mutate(
    timepoint = if_else(when == "pre", "T1", "T2"),
    type1_share = 100 * fibres_type1 / fibres_mixed
  ) |>
  select(-when)
traits <- setdiff(names(biopsy), c("subject", "arm", "timepoint"))
per_subject <- biopsy |>
  pivot_longer(all_of(traits), names_to = "trait") |>
  pivot_wider(names_from = timepoint) |>
  transmute(subject, trait, d = T2 - T1, pct = 100 * (T2 - T1) / T1) |>
  pivot_wider(names_from = trait, values_from = c(d, pct)) |>
  left_join(select(phenotype, subject, comp_hypertrophy, volume_load_total_kg), by = "subject")
outcomes <- c(
  "comp_hypertrophy", "volume_load_total_kg", paste0("d_", traits), paste0("pct_", traits)
)
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
  # The change score's floor applies here too: a protein seen in half the biopsies can fit a
  # large t through a handful of points.
  y <- values[, s$sample_id]
  y <- y[rowSums(!is.na(y)) >= ceiling(2 * ncol(y) / 3), ]
  correlation <- duplicateCorrelation(y, design, block = s$subject)$consensus.correlation
  fit <- lmFit(y, design, block = s$subject, correlation = correlation) |>
    eBayes(robust = TRUE)
  map(c("within", "between"), \(term) {
    topTable(fit, coef = term, number = Inf, sort.by = "none") |>
      as_tibble(rownames = "feature") |>
      transmute(
        trait, term, feature,
        n_samples = nrow(s), effect = logFC, t, p = P.Value, fdr = adj.P.Val
      )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(label = feature_label[feature], .after = feature)

by_subject <- function(timepoint) {
  rows <- filter(meta, timepoint == .env$timepoint)
  set_names(as.data.frame(values[, rows$sample_id, drop = FALSE]), rows$subject) |> as.matrix()
}
change <- function(from, to) {
  a <- by_subject(from)
  b <- by_subject(to)
  both <- intersect(colnames(a), colnames(b))
  b[, both, drop = FALSE] - a[, both, drop = FALSE]
}
windows <- list(
  training = change("T1", "T2"), baseline = by_subject("T1"), acute = change("T2", "T3")
)

# cor.test(exact = FALSE)'s p, step for step: S statistic, rho back from S, then the t tail.
spearman_p <- function(rho, n) {
  s <- (n^3 - n) * (1 - rho) / 6
  r <- 1 - s / ((n * (n^2 - 1)) / 6)
  t <- r / sqrt((1 - r^2) / (n - 2))
  pmin(2 * if_else(s > (n^3 - n) / 6, pt(t, n - 2), pt(t, n - 2, lower.tail = FALSE)), 1)
}
# A protein is tested only when two thirds of the window's subjects saw it.
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
  mutate(fdr = p.adjust(p, "BH"), .by = c(window, outcome)) |>
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
    )
)

# Units come from the phenotype column names; type1_share is a percentage.
unit_cols <- grep("_pre_", names(phenotype), value = TRUE)
unit_of <- c(set_names(sub(".*_pre_", "", unit_cols), sub("_pre_.*", "", unit_cols)),
  type1_share = "%"
)
outcome_unit <- c(
  comp_hypertrophy = "score", volume_load_total_kg = "kg",
  set_names(unit_of[traits], paste0("d_", traits)),
  set_names(rep("%", length(traits)), paste0("pct_", traits))
)
window_axis <- c(training = "change, T1 to T2", baseline = "T1 level", acute = "change, T2 to T3")
arm_colours <- c(HR = "#2166AC", LR = "#B2182B")
# Blank slots keep a short last page on the same 3 x 4 panel size.
slots <- \(panels) c(panels, strrep(" ", seq_len(12 - length(panels))))
page_theme <- theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 7), plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0), legend.position = "top"
  )
# Twelve panels to a page, strongest first within each group: p, then |t| or |rho|, then ID.
paginate <- function(hits, by) {
  hits |>
    mutate(
      page = (row_number() - 1) %/% 12 + 1, n_pages = max(page),
      position = (row_number() - 1) %% 12 + 1, .by = all_of(by)
    ) |>
    mutate(pdf_page = consecutive_id(pick(all_of(by)), page))
}

trait_samples <- map(set_names(traits), \(trait) filter(samples, !is.na(.data[[trait]])))
seen_in <- \(feature, trait) !is.na(values[feature, trait_samples[[trait]]$sample_id])
model_hits <- model |>
  filter(p < 0.05) |>
  mutate(
    n_biopsies = map2_int(feature, trait, \(f, tr) sum(seen_in(f, tr))),
    n_subjects = map2_int(feature, trait, \(f, tr) {
      n_distinct(trait_samples[[tr]]$subject[seen_in(f, tr)])
    }),
    panel = sprintf(
      "%s\n%s %.3g  p %.2g  q %.2g\n%d biopsies, %d subjects", label, term,
      effect, p, fdr, n_biopsies, n_subjects
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
      subtitle = sprintf(
        "%d of %d %ss at p < 0.05, %d at BH < 0.05. Page %d of %d.",
        info$n_nominal, info$n_tested, feature_noun, info$n_fdr_05, hits$page[1],
        hits$n_pages[1]
      ),
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
  mutate(panel = sprintf(
    "%s vs %s (%s)\nrho %+.2f  p %.2g  q %.2g  n %d", label, outcome,
    outcome_unit[outcome], rho, p, fdr, n_subjects
  )) |>
  arrange(match(window, names(windows)), p, desc(abs(rho)), feature, outcome) |>
  paginate("window")
change_page <- function(hits) {
  w <- hits$window[1]
  info <- counts |>
    filter(analysis == paste("change score:", w)) |>
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
      title = sprintf("Change score, %s window", w),
      subtitle = sprintf(
        "%d of %d %s-outcome pairs at p < 0.05, %d at BH < 0.05. Page %d of %d.",
        info$n_nominal, info$n_tested, feature_noun, info$n_fdr_05, hits$page[1],
        hits$n_pages[1]
      ),
      x = sprintf("%s, %s", value_label, window_axis[[w]]), y = "outcome (unit in header)",
      caption = stringr::str_wrap(width = 160, paste(
        "One point per subject. Lines are a least-squares fit per arm, drawn as a guide; the",
        "tests are Spearman. Header: rho and p pooled over both arms,",
        "BH q within window and outcome, and the subjects used.",
        "Corner: each arm's rho and exact p. Arm is the median split of",
        "comp_hypertrophy, so pooled correlations with it and with fibre area change partly",
        "carry the group split. Table: c_data/04_associate.xlsx, change_score and",
        "change_score_by_arm."
      ))
    ) +
    page_theme
}
model_pages <- map(split(model_hits, model_hits$pdf_page), model_page)
change_pages <- map(split(change_hits, change_hits$pdf_page), change_page)

# Overviews: every feature by trait (sample model t) and by outcome (change-score rho, one page
# per window). Each effect map sits beside a status map; grey cells were not tested.
ht_opt$message <- FALSE
status_of <- \(p, fdr) {
  case_when(fdr < 0.05 ~ "BH < 0.05", p < 0.05 ~ "p < 0.05", !is.na(p) ~ "not significant")
}
status_colours <- c("not significant" = "white", "p < 0.05" = "#35978F", "BH < 0.05" = "black")
status_legend <- Legend(
  labels = c(names(status_colours), "not tested"), title = "status",
  legend_gp = gpar(fill = c(status_colours, "grey80"))
)
diverging <- c("#5E3C99", "white", "#E66101")
cell_grid <- function(d, rows, cols, col_var, value) {
  m <- matrix(NA, length(rows), length(cols), dimnames = list(rows, cols))
  d <- filter(d, feature %in% rows, .data[[col_var]] %in% cols)
  m[cbind(d$feature, d[[col_var]])] <- value(d)
  m
}
body <- function(m, name, col, title, legend_at = NULL, legend = TRUE, names = FALSE) {
  rows <- rownames(m)
  Heatmap(m,
    name = name, col = col, na_col = "grey80", heatmap_legend_param = list(at = legend_at),
    show_heatmap_legend = legend, cluster_rows = FALSE, cluster_columns = FALSE,
    show_row_names = names && length(rows) <= 80,
    row_labels = gsub("\n", " ", feature_label[rows]),
    row_names_gp = gpar(fontsize = if (length(rows) <= 40) 7 else 5),
    column_names_gp = gpar(fontsize = 7), column_title = title,
    column_title_gp = gpar(fontsize = 8), use_raster = length(rows) > 200
  )
}
overview <- function(maps, title, caption) {
  list(title = title, draw = function() {
    draw(maps,
      column_title = title, column_title_gp = gpar(fontsize = 11, fontface = "bold"),
      merge_legend = TRUE, heatmap_legend_list = list(status_legend),
      padding = unit(c(14, 4, 4, 4), "mm")
    )
    grid.text(stringr::str_wrap(caption, 190),
      x = unit(4, "mm"), y = unit(3, "mm"),
      just = c("left", "bottom"), gp = gpar(fontsize = 7)
    )
  })
}
t_cap <- unname(quantile(abs(model$t), 0.99, na.rm = TRUE))
model_overview <- function(rows, group) {
  t_colours <- circlize::colorRamp2(c(-t_cap, 0, t_cap), diverging)
  term_grid <- \(term, value) {
    cell_grid(filter(model, term == .env$term), rows, traits, "trait", value)
  }
  maps <- body(term_grid("within", \(d) d$t), "t", t_colours, "Within-person t",
    legend_at = signif(c(-t_cap, 0, t_cap), 2)
  ) +
    body(term_grid("between", \(d) d$t), "t, between", t_colours, "Between-person t",
      legend = FALSE
    ) +
    body(term_grid("within", \(d) status_of(d$p, d$fdr)), "status, within", status_colours,
      "Within: p and BH q",
      legend = FALSE
    ) +
    body(term_grid("between", \(d) status_of(d$p, d$fdr)), "status, between", status_colours,
      "Between: p and BH q",
      legend = FALSE, names = TRUE
    )
  overview(maps, sprintf("Sample model overview: %s", group), paste(
    sprintf("All %s %ss, one row each.", format(length(rows), big.mark = ","), feature_noun),
    "Moderated t of the within and between terms, abundance ~ timepoint + between + within,",
    sprintf("subject blocked. Colours stop at the 99th percentile of |t|, %.2g.", t_cap),
    "Rows are ordered by average-linkage clustering (1 - Pearson r) of the",
    "imputed matrix, for display only.",
    "Grey: not tested. Table: c_data/04_associate.xlsx, model."
  ))
}
change_overview <- function(rows, group, window) {
  d <- filter(change_score, window == .env$window)
  maps <- body(cell_grid(d, rows, outcomes, "outcome", \(d) d$rho), "rho",
    circlize::colorRamp2(c(-1, 0, 1), diverging), "Spearman rho, both arms pooled",
    legend_at = c(-1, 0, 1)
  ) +
    body(cell_grid(d, rows, outcomes, "outcome", \(d) status_of(d$p, d$fdr)), "status",
      status_colours, "p and BH q",
      legend = FALSE, names = TRUE
    )
  overview(maps, sprintf("Change score overview, %s window: %s", window, group), paste(
    sprintf("All %s %ss, one row each.", format(length(rows), big.mark = ","), feature_noun),
    sprintf(
      "Spearman rho between the %s (%s) and each outcome.", value_label, window_axis[[window]]
    ),
    "Rows are ordered by average-linkage clustering (1 - Pearson r) of the",
    "imputed matrix, for display only.",
    "Grey: not tested. Table: c_data/04_associate.xlsx, change_score."
  ))
}
overviews <- imap(row_groups, \(rows, group) {
  c(list(model_overview(rows, group)), map(names(windows), \(w) change_overview(rows, group, w)))
}) |>
  list_flatten()

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
  )
)
stopifnot(
  !anyDuplicated(select(filter(plot_index, analysis != "overview"), -label, -pdf_page, -position)),
  max(plot_index$pdf_page) == length(overviews) + length(model_pages) + length(change_pages),
  nrow(anti_join(filter(model, p < 0.05), model_hits, by = c("trait", "term", "feature"))) == 0,
  nrow(anti_join(filter(change_score, p < 0.05), change_hits,
    by = c("window", "outcome", "feature")
  )) == 0,
  nrow(model_hits) == sum(model$p < 0.05, na.rm = TRUE),
  nrow(change_hits) == sum(change_score$p < 0.05)
)

saveRDS(list(change_score = change_score, model = model), file.path(out, "associate.rds"),
  compress = "xz"
)
pdf(file.path(figures, "04_associate_figures.pdf"), width = 11, height = 8.5)
walk(overviews, \(page) page$draw())
walk(c(model_pages, change_pages), print)
invisible(dev.off())

# The within-arm table runs to hundreds of thousands of rows; the workbook keeps p < 0.05.
sheets <- list(
  counts = counts,
  model = arrange(model, trait, term, p),
  change_score = arrange(change_score, window, outcome, p),
  change_score_by_arm = filter(by_arm, p < 0.05),
  plot_index = plot_index,
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Tested, p < 0.05 and BH < 0.05 counts per model term and trait, and per window and outcome.",
  "Sample model: effect per trait unit, t, p and BH FDR per protein, trait and term.",
  "Change score: pooled Spearman rho, p and BH FDR per protein, window and outcome.",
  "Change score within each arm, exact p, rows with p < 0.05.",
  "PDF page of each overview, and page and panel position of every p < 0.05 row.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "04_associate.xlsx"))
