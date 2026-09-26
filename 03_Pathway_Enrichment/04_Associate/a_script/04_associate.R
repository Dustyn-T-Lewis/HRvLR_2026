# Ties each gene set's singscore to phenotype, the protein level's two ways. Each collection is its
# own screen: BH and the chance count run within collection.
#
# Sample model: every T1 and T2 biopsy carries its own trait values. Per trait, limma fits
# abundance ~ timepoint + between + within with subject blocked, where between is the subject's
# mean trait and within the biopsy's deviation from it. within asks whether a protein rises in the
# biopsies where that person's trait rose; between asks whether people with more of the trait
# carry more of the protein. Arm stays out: it was defined from one of these traits.
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
})

stage <- here("03_Pathway_Enrichment", "04_Associate")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  singscore = "03_Pathway_Enrichment/01_Scores/c_data/singscore.rds",
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  proteins = "01_Preprocess/03_Imputation/c_data/DAList_imputed.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
catalog <- readRDS(paths[["gene_sets"]])$set_catalog |> filter(qualifies)
values <- readRDS(paths[["singscore"]])$scores[catalog$set_id, ]
meta <- as_tibble(readRDS(paths[["proteins"]])$metadata)
stopifnot(identical(colnames(values), meta$sample_id))
feature_label <- with(catalog, set_names(enrichVolcano::ev_clean_label(pathway), set_id))
collection_of <- with(catalog, set_names(collection, set_id))
feature_noun <- "set"
value_label <- "singscore"

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
  y <- values[, s$sample_id]
  correlation <- duplicateCorrelation(y, design, block = s$subject)$consensus.correlation
  fit <- lmFit(y, design, block = s$subject, correlation = correlation) |>
    eBayes(robust = TRUE)
  map(c("within", "between"), \(term) {
    topTable(fit, coef = term, number = Inf, sort.by = "none") |>
      as_tibble(rownames = "feature") |>
      transmute(trait, term, feature, n_samples = nrow(s), effect = logFC, t, p = P.Value)
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(collection = collection_of[feature], label = feature_label[feature], .after = feature) |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(trait, term, collection))

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
  mutate(collection = collection_of[feature], label = feature_label[feature], .after = feature) |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(window, outcome, collection))
by_arm <- map(set_names(c("HR", "LR")), \(arm) {
  imap(windows, \(x, window) {
    kept <- x[, arm_of[colnames(x)] == arm, drop = FALSE]
    mutate(correlate(kept, exact = TRUE), window = window, .before = 1)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "arm") |>
  filter(!is.na(rho))

chance_expectation <- bind_rows(
  summarise(model,
    n_tested = sum(!is.na(p)), n_nominal = sum(p < 0.05, na.rm = TRUE),
    n_fdr_05 = sum(fdr < 0.05, na.rm = TRUE), .by = c(term, trait, collection)
  ) |>
    transmute(
      analysis = paste("model:", term), comparison = trait, collection, n_tested, n_nominal,
      n_fdr_05
    ),
  summarise(change_score,
    n_tested = n(), n_nominal = sum(p < 0.05), n_fdr_05 = sum(fdr < 0.05),
    .by = c(window, outcome, collection)
  ) |>
    transmute(
      analysis = paste("change score:", window), comparison = outcome, collection, n_tested,
      n_nominal, n_fdr_05
    )
) |>
  mutate(n_expected = 0.05 * n_tested, ratio = round(n_nominal / n_expected, 2), .after = n_nominal)

# Twelve panels to a page; split() pages the data so each page builds only its own panels.
paged <- function(data, draw) {
  pages <- split(data, (as.integer(data$panel) - 1) %/% 12)
  imap(unname(pages), \(page, i) {
    draw(droplevels(page)) +
      facet_wrap(~panel, nrow = 3, ncol = 4, scales = "free") +
      labs(tag = sprintf("page %d of %d", i, length(pages))) +
      theme_minimal(base_size = 9) +
      theme(
        strip.text = element_text(size = 7), plot.caption = element_text(hjust = 0),
        plot.tag.position = "topright", plot.tag = element_text(size = 8), legend.position = "top"
      )
  })
}
arm_colours <- c(HR = "#2166AC", LR = "#B2182B")

model_pages <- model |>
  filter(p < 0.05) |>
  arrange(match(trait, traits), desc(term), p) |>
  group_split(trait, term, .keep = TRUE) |>
  map(\(hits) {
    info <- chance_expectation |>
      filter(analysis == paste("model:", hits$term[1]), comparison == hits$trait[1]) |>
      summarise(across(c(n_tested, n_nominal, n_expected), sum))
    points <- hits |>
      mutate(panel = sprintf(
        "%s: %s\n%s %.3g  p %.2g  q %.2g", collection, label, term, effect, p, fdr
      )) |>
      select(feature, panel) |>
      cross_join(
        select(samples, sample_id, subject, timepoint, trait_value = all_of(hits$trait[1]))
      ) |>
      mutate(
        value = values[cbind(feature, sample_id)], arm = arm_of[subject],
        panel = factor(panel, unique(panel))
      ) |>
      filter(!is.na(value), !is.na(trait_value))
    paged(points, \(page) {
      ggplot(page, aes(trait_value, value, colour = arm)) +
        geom_line(
          aes(group = subject),
          data = \(d) filter(d, n() == 2, .by = c(panel, subject)), alpha = 0.6
        ) +
        geom_point(aes(shape = timepoint), size = 1.5) +
        scale_colour_manual(values = arm_colours, name = NULL) +
        scale_shape_manual(values = c(T1 = 1, T2 = 16), name = NULL) +
        labs(
          title = sprintf(
            "%s against %s, %s-person term", value_label, hits$trait[1], hits$term[1]
          ),
          x = paste(hits$trait[1], "at the biopsy"), y = value_label,
          subtitle = paste(
            "One line per subject from T1 to T2: within is the slope of the lines,",
            "between the spread of subjects"
          ),
          caption = sprintf(
            paste(
              "%d of %d %ss reach p < 0.05 on this term; chance predicts %.0f.",
              "limma on %d biopsies, subject blocked; q from BH within trait, term and collection.",
              "Table: c_data/04_associate.xlsx, model."
            ),
            info$n_nominal, info$n_tested, feature_noun, info$n_expected, hits$n_samples[1]
          )
        )
    })
  }) |>
  list_flatten()

arm_text <- by_arm |>
  mutate(text = sprintf("%-3s n=%d rho=%+.2f p=%.3f", arm, n_subjects, rho, p)) |>
  summarise(caption = paste(text, collapse = "\n"), .by = c(window, outcome, feature))
change_pages <- change_score |>
  filter(p < 0.05) |>
  arrange(match(window, names(windows)), p) |>
  group_split(window, .keep = TRUE) |>
  map(\(hits) {
    w <- hits$window[1]
    n_tested <- sum(filter(chance_expectation, analysis == paste("change score:", w))$n_tested)
    points <- hits |>
      left_join(arm_text, by = c("window", "outcome", "feature")) |>
      mutate(panel = sprintf(
        "%s: %s vs %s\nrho %+.2f  p %.2g  q %.2g", collection, label, outcome, rho, p, fdr
      )) |>
      select(feature, outcome, panel, caption) |>
      cross_join(tibble(subject = colnames(windows[[w]]))) |>
      mutate(
        value = windows[[w]][cbind(feature, subject)],
        outcome_value = map2_dbl(outcome, subject, outcome_of), arm = arm_of[subject],
        panel = factor(panel, unique(panel))
      ) |>
      filter(!is.na(value), !is.na(outcome_value))
    paged(points, \(page) {
      ggpubr::ggscatter(page, "value", "outcome_value",
        color = "arm", palette = arm_colours, add = "reg.line", size = 1.3
      ) +
        geom_text(
          data = distinct(page, panel, caption), aes(-Inf, Inf, label = caption),
          inherit.aes = FALSE, hjust = -0.03, vjust = 1.1, size = 1.9, family = "mono"
        ) +
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.45))) +
        labs(
          title = sprintf("%s against phenotype, %s window", value_label, w),
          x = sprintf("%s, %s window, one point per subject", value_label, w), y = "outcome",
          colour = NULL,
          caption = stringr::str_wrap(width = 160, sprintf(
            paste(
              "%d of %d %s-outcome pairs reach pooled p < 0.05; chance predicts %.0f. Header:",
              "pooled Spearman, q from BH within window, outcome and collection. Corner: each",
              "arm's rho and exact p. Table: c_data/04_associate.xlsx, change_score and",
              "change_score_by_arm."
            ),
            nrow(hits), n_tested, feature_noun, 0.05 * n_tested
          ))
        )
    })
  }) |>
  list_flatten()

chance_figure <- chance_expectation |>
  summarise(across(c(n_nominal, n_expected), sum), .by = c(analysis, comparison)) |>
  mutate(
    ratio = n_nominal / n_expected,
    comparison = factor(comparison, rev(unique(c(traits, outcomes))))
  ) |>
  ggplot(aes(ratio, comparison)) +
  geom_vline(xintercept = 1, colour = "grey40") +
  geom_col(fill = "grey60", width = 0.7) +
  facet_wrap(~analysis, nrow = 1, scales = "free_y") +
  labs(
    title = "Nominal associations against chance", x = "nominal hits / chance expectation",
    y = NULL, caption = sprintf(
      "Collections pooled; chance is 5%% of %ss tested. Per collection: c_data/04_associate.xlsx.",
      feature_noun
    )
  ) +
  theme_minimal(base_size = 8) +
  theme(plot.caption = element_text(hjust = 0))

stopifnot(
  sum(map_int(model_pages, \(g) nlevels(g$data$panel))) == sum(model$p < 0.05, na.rm = TRUE),
  sum(map_int(change_pages, \(g) nlevels(g$data$panel))) == sum(change_score$p < 0.05)
)

pdf(file.path(figures, "04_associate_figures.pdf"), width = 11, height = 8.5)
walk(c(list(chance_figure), model_pages, change_pages), print)
invisible(dev.off())

# The within-arm table runs to hundreds of thousands of rows; the workbook keeps p < 0.05.
sheets <- list(
  chance_expectation = chance_expectation,
  model = arrange(model, trait, term, p),
  change_score = arrange(change_score, window, outcome, p),
  change_score_by_arm = filter(by_arm, p < 0.05),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
)
read_me <- tibble(sheet = names(sheets), holds = c(
  "Nominal hits against chance, per model term and trait and per change-score window and outcome.",
  "Sample model: effect per trait unit, t, p and BH FDR per set, trait and term.",
  "Change score: pooled Spearman rho, p and BH FDR per set, window and outcome.",
  "Change score within each arm, exact p, rows with p < 0.05.",
  "Files read, with md5.",
  "Packages loaded at run time."
))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "04_associate.xlsx"))
