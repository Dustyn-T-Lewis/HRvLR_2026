# The chance line for the candidate counts. Each shuffle hands every subject another subject's
# phenotype row, arm and traits together, so trait correlations and arm sizes hold and only the
# link between a person's proteome and their response breaks. The tests behind the candidate
# lists are rerun on the shuffled labels and hits are counted with the same list rules, so a list
# can be read against what pure chance produces. Context only: nothing is filtered.
#
# Shuffle 0 is the real labelling run through this same code, so observed and permuted counts
# share one code path. Two simplifications keep 200 shuffles tractable: within-arm Spearman p uses
# the t approximation rather than the exact permutation p, and set scores (singscore) and module
# eigengenes are taken as computed once, not rebuilt per shuffle.

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
source(here("R", "helpers.R"))

stage <- here("05_Summary", "01_Permutation")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  proteins = "01_Preprocess/03_Quantification/c_data/proteins.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  singscore = "03_Pathway_Enrichment/01_Scores/c_data/singscore.rds",
  modules = "04_Network/01_Build_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
n_shuffles <- 200
cores <- min(8, parallel::detectCores())

quantified <- readRDS(paths[["proteins"]])
proteins <- quantified$proteins
observed <- proteins$other$n.observations > 0
d <- readRDS(paths[["design"]])
contrast_matrix <- d$dal$design$contrast_matrix
key_contrasts <- c(hr_lr_contrasts, "Baseline_HRvLR")
gs <- readRDS(paths[["gene_sets"]])
selected <- filter(gs$protein_map, selected)
set_rows <- map(gs$sets, \(genes) {
  match(selected$uniprot_id[match(genes, selected$gene)], rownames(proteins$E))
})
catalog <- filter(gs$set_catalog, qualifies)
set_scores <- readRDS(paths[["singscore"]])$scores[catalog$set_id, ]
modules <- readRDS(paths[["modules"]])
module_members <- filter(modules$membership, module != "grey")
module_rows <- split(match(module_members$uniprot_id, rownames(proteins$E)), module_members$module)

base_meta <- quantified$metadata
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
traits <- c("fcsa_mixed", "fcsa_type1", "fcsa_type2", "mcsa", "leg_press_1rm", "leg_ext_1rm")
percent <- map(set_names(traits), \(trait) {
  pre <- phenotype[[grep(paste0("^", trait, "_pre"), names(phenotype))]]
  post <- phenotype[[grep(paste0("^", trait, "_post"), names(phenotype))]]
  100 * (post - pre) / pre
})
responses <- tibble(
  subject = phenotype$subject, arm = phenotype$arm,
  comp_hypertrophy = phenotype$comp_hypertrophy,
  volume_load_total_kg = phenotype$volume_load_total_kg,
  !!!set_names(percent, paste0("pct_", traits))
)
outcomes <- setdiff(names(responses), c("subject", "arm"))
strength_tests <- c("pct_leg_press_1rm", "pct_leg_ext_1rm")

# One labelling's evidence: every test that can place a feature on a list, as a long table of
# level, feature, contrast, task, window and p, in the shape on_lists() reads.
evidence_for <- function(responses) {
  meta <- base_meta |>
    select(-arm) |>
    left_join(select(responses, subject, arm), by = "subject") |>
    mutate(group = factor(paste(arm, timepoint, sep = "_"), levels(base_meta$group)))
  design <- model.matrix(~ 0 + group + blood, meta)
  colnames(design) <- sub("^group", "", colnames(design))
  arm_of <- set_names(responses$arm, responses$subject)
  corrected <- blood_corrected(proteins$E, meta$blood, model.matrix(~ 0 + group, meta))
  levels <- list(
    proteins = replace(corrected, !observed, NA), sets = set_scores, modules = modules$eigengenes
  )

  protein_fit <- fit_limpa(proteins, design, contrast_matrix[, key_contrasts], meta$subject)
  contrast_rows_of <- function(rows, prefix = "") {
    mutate(rows, contrast = sub(paste0("^", prefix), "", test))
  }
  fgsea_rows <- map(key_contrasts, \(contrast) {
    stats <- set_names(protein_fit$t[selected$uniprot_id, contrast], selected$gene)
    set.seed(1)
    fgsea::fgsea(gs$sets, stats, minSize = 15, maxSize = 500) |>
      as_tibble() |>
      transmute(feature = pathway, test = paste("fgsea", contrast), p = pval)
  }) |>
    list_rbind()
  contrasts <- bind_rows(
    proteins = contrast_rows(protein_fit) |> contrast_rows_of(),
    sets = bind_rows(
      fgsea_rows,
      fry_rows(proteins$E, set_rows, design, contrast_matrix[, key_contrasts], meta$subject)
    ) |>
      contrast_rows_of("(fgsea|fry) "),
    modules = bind_rows(
      fit_contrasts(modules$eigengenes, design, contrast_matrix[, key_contrasts], meta$subject) |>
        contrast_rows("eigengene "),
      fry_rows(proteins$E, module_rows, design, contrast_matrix[, key_contrasts], meta$subject)
    ) |>
      contrast_rows_of("(eigengene|fry) "),
    .id = "level"
  )

  per_level <- imap(levels, \(values, level) {
    by_time <- map(set_names(c("T1", "T2", "T3")), \(tp) subject_values(values, meta, tp))
    windows <- list(
      training = subject_change(values, meta, "T1", "T2"), baseline = by_time$T1,
      acute = subject_change(values, meta, "T2", "T3")
    )
    tasks <- list(
      Baseline_HRvLR = by_time$T1, Trained_HRvLR = by_time$T2,
      Training_change_HRvLR = windows$training, Acute_change_HRvLR = windows$acute
    )
    classify <- imap(tasks, \(x, task) {
      hr <- x[, arm_of[colnames(x)] == "HR", drop = FALSE]
      lr <- x[, arm_of[colnames(x)] == "LR", drop = FALSE]
      tested <- rowSums(!is.na(hr)) >= 3 & rowSums(!is.na(lr)) >= 3
      tibble(
        feature = rownames(x)[tested], task,
        p = suppressWarnings(matrixTests::row_wilcoxon_twosample(hr[tested, ], lr[tested, ])$pvalue)
      )
    }) |>
      list_rbind()
    associate <- imap(windows, \(x, window) {
      outcome_of <- \(name) set_names(responses[[name]], responses$subject)[colnames(x)]
      arm <- arm_of[colnames(x)]
      pooled <- map(outcomes, \(name) {
        spearman_by_row(x, outcome_of(name)) |> mutate(test = paste("pooled", name))
      })
      within <- map(c("HR", "LR"), \(a) {
        map(outcomes, \(name) {
          spearman_by_row(x[, arm == a, drop = FALSE], outcome_of(name)[arm == a]) |>
            mutate(test = paste(a, name))
        })
      }) |>
        list_flatten()
      slopes <- map(setdiff(outcomes, "comp_hypertrophy"), \(name) {
        predictors <- data.frame(arm = factor(arm, c("LR", "HR")), outcome = outcome_of(name))
        seen <- !is.na(x) & rep(!is.na(predictors$outcome), each = nrow(x))
        per_arm <- sapply(c("HR", "LR"), \(a) rowSums(seen[, arm %in% a, drop = FALSE]))
        subject_fit(
          x[apply(per_arm, 1, min) >= 5, , drop = FALSE], predictors, ~ arm * outcome,
          "armHR:outcome"
        ) |>
          mutate(test = paste("slope", name))
      })
      two_predictor <- map(strength_tests, \(test) {
        predictors <- data.frame(
          growth = outcome_of("comp_hypertrophy"), strength = outcome_of(test)
        )
        subject_fit(x, predictors, ~ growth + strength, c("growth", "strength")) |>
          mutate(test = paste(coefficient, test))
      })
      bind_rows(c(pooled, within, slopes, two_predictor)) |>
        transmute(feature, window, test, p)
    }) |>
      list_rbind()
    ancova_rows <- ancova(values, meta) |>
      transmute(feature, window = paste0("ancova_", window), test = paste("ANCOVA", window), p)
    bind_rows(classify, associate, ancova_rows)
  }) |>
    list_rbind(names_to = "level")

  bind_rows(contrasts, per_level) |>
    filter(!is.na(p))
}

# Hits per feature on each list for one labelling.
hit_counts <- function(responses) {
  evidence_for(responses) |>
    on_lists() |>
    summarise(hits = sum(p < 0.05), .by = c(level, list, feature))
}

observed_hits <- hit_counts(responses)
shuffled <- parallel::mclapply(seq_len(n_shuffles), \(i) {
  set.seed(i)
  donor <- sample(nrow(responses))
  hit_counts(mutate(responses[donor, ], subject = responses$subject))
}, mc.cores = cores)
stopifnot(length(shuffled) == n_shuffles, !any(map_lgl(shuffled, inherits, "try-error")))
permuted_hits <- list_rbind(shuffled, names_to = "shuffle")

thresholds <- 1:8
at_least <- function(hits) {
  map(thresholds, \(k) summarise(hits, k = k, n_features = sum(hits >= k), .by = c(level, list))) |>
    list_rbind()
}
permuted_counts <- map(split(permuted_hits, permuted_hits$shuffle), at_least) |>
  list_rbind(names_to = "shuffle")
chance_line <- at_least(observed_hits) |>
  rename(observed = n_features) |>
  left_join(
    summarise(permuted_counts,
      permuted_median = median(n_features), permuted_95 = quantile(n_features, 0.95),
      .by = c(level, list, k)
    ),
    by = c("level", "list", "k")
  ) |>
  left_join(
    permuted_counts |>
      left_join(
        rename(at_least(observed_hits), observed = n_features),
        by = c("level", "list", "k")
      ) |>
      summarise(share_at_or_above = mean(n_features >= observed), .by = c(level, list, k)),
    by = c("level", "list", "k")
  )

# A feature's score: the share of shuffles in which it reached at least its observed hit count.
feature_scores <- observed_hits |>
  rename(observed_hits = hits) |>
  left_join(
    permuted_hits |>
      left_join(rename(observed_hits, observed_hits = hits), by = c("level", "list", "feature")) |>
      summarise(
        permutation_score = (sum(hits >= observed_hits) + 1) / (n_shuffles + 1),
        .by = c(level, list, feature)
      ),
    by = c("level", "list", "feature")
  ) |>
  mutate(permutation_score = coalesce(permutation_score, 1 / (n_shuffles + 1))) |>
  arrange(level, list, permutation_score)

chance_figure <- chance_line |>
  ggplot(aes(k)) +
  geom_ribbon(aes(ymin = 0, ymax = permuted_95), fill = "grey85") +
  geom_line(aes(y = permuted_median), colour = "grey45", linetype = 2) +
  geom_line(aes(y = observed), colour = "firebrick") +
  geom_point(aes(y = observed), colour = "firebrick", size = 1.2) +
  facet_grid(level ~ list, scales = "free_y") +
  scale_y_continuous(trans = "log1p") +
  labs(
    title = "Candidate counts against the chance line",
    x = "at least k nominal hits", y = "features (log scale)",
    caption = stringr::str_wrap(width = 170, paste(
      sprintf("Grey band: up to the 95th percentile over %d shuffles of the", n_shuffles),
      "subjects' phenotype rows; dashed: median; red: the real labelling through the same code.",
      "Red above the band means more features reach k hits than chance gives. Context only;",
      "nothing is filtered. Table: c_data/01_permutation.xlsx, chance_line."
    ))
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.caption = element_text(hjust = 0))
pdf(file.path(figures, "01_permutation_figures.pdf"), width = 11, height = 8.5)
print(chance_figure)
invisible(dev.off())

write_workbook(
  list(chance_line = chance_line, feature_scores = feature_scores),
  c(
    paste(
      "Per level, list and k: features with at least k hits under the real labelling, the permuted",
      "median and 95th percentile, and the share of shuffles at or above the real count."
    ),
    paste(
      "Per feature: its hits under the real labelling (this step's code) and the share of",
      "shuffles reaching that many, (b + 1) / (m + 1)."
    )
  ),
  file.path(out, "01_permutation.xlsx"),
  inputs, paths
)
