# Candidate tables: every protein, set and module with at least one p < 0.05 among the tests that
# bear on HR against LR, one table for training and one for the acute bout. Every nominal test
# counts once and nothing is weighted, so hit counts from correlated tests overlap. The same rules
# at T1, before any training, give the reference distribution of hit counts.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(readxl)
  library(writexl)
})
source(here("R", "helpers.R"))

out <- here("05_Integration", "02_Candidates", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  protein_contrasts = "02_Proteins/02_Contrasts/c_data/02_contrasts.xlsx",
  protein_classify = "02_Proteins/03_Classify/c_data/03_classify.xlsx",
  protein_associate = "02_Proteins/04_Associate/c_data/04_associate.xlsx",
  set_tests = "03_Pathways/03_Contrasts/c_data/set_tests.rds",
  set_contrasts = "03_Pathways/03_Contrasts/c_data/03_contrasts.xlsx",
  set_classify = "03_Pathways/04_Classify/c_data/04_classify.xlsx",
  set_associate = "03_Pathways/05_Associate/c_data/05_associate.xlsx",
  gene_sets = "03_Pathways/01_Gene_Sets/c_data/gene_sets.rds",
  modules = "04_Modules/01_Build/c_data/01_build.xlsx",
  module_contrasts = "04_Modules/02_Contrasts/c_data/02_contrasts.xlsx",
  module_classify = "04_Modules/03_Classify/c_data/03_classify.xlsx",
  module_associate = "04_Modules/04_Associate/c_data/04_associate.xlsx",
  permutation = "05_Integration/01_Chance_Line/c_data/01_chance_line.xlsx"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
sheet <- \(input, name) suppressMessages(read_excel(paths[[input]], name))

# One long table per level: feature, family, test, effect, p with its BH and Storey q (each from
# its own screen), and the window or contrast that places it on a list. Within-arm correlations
# keep only their p < 0.05 rows upstream, so they carry no correction.
associations <- function(input) {
  bind_rows(
    sheet(input, "change_score") |>
      transmute(feature, window, test = paste("pooled rho", outcome), effect = rho, p, fdr, qvalue),
    sheet(input, "change_score_by_arm") |>
      transmute(feature, window, test = paste(arm, "rho", outcome), effect = rho, p),
    sheet(input, "arm_difference") |>
      transmute(feature, window, test = paste("slope HR-LR", outcome), effect, p, fdr, qvalue),
    sheet(input, "two_predictor") |>
      transmute(
        feature, window,
        test = paste(coefficient, "term,", sub("pct_", "", strength_test), "model"), effect, p,
        fdr, qvalue
      )
  ) |>
    mutate(family = "associate")
}
classify <- function(input) {
  sheet(input, "auc") |>
    filter(!paired) |>
    transmute(
      feature, task,
      test = paste("classify", task), effect = auc, p, fdr, qvalue, family = "classify"
    )
}
ancova_rows <- function(input) {
  sheet(input, "ancova") |>
    transmute(feature,
      test = paste("ANCOVA", window), window = paste0("ancova_", window), effect, p, fdr, qvalue,
      family = "contrast"
    )
}

protein_de <- sheet("protein_contrasts", "DEP_matrix") |>
  select(feature = uniprot_id, matches("_(logFC|p|fdr|qvalue)$")) |>
  pivot_longer(-feature,
    names_to = c("contrast", ".value"), names_pattern = "(.*)_(logFC|p|fdr|qvalue)$"
  ) |>
  transmute(feature, contrast, test = contrast, effect = logFC, p, fdr, qvalue, family = "contrast")
# cameraPR corroborates fgsea and is not counted as a hit of its own.
set_de <- readRDS(paths[["set_tests"]])$set_tests |>
  filter(method != "cameraPR") |>
  transmute(
    feature = set_id, contrast, test = paste(method, contrast),
    effect = if_else(method == "fgsea", nes, if_else(direction == "Up", 1, -1)), p, fdr, qvalue,
    family = "contrast"
  )
module_de <- bind_rows(
  sheet("module_contrasts", "eigengene_tests") |>
    transmute(
      feature = module, contrast,
      test = paste("eigengene", contrast), effect = logFC, p, fdr, qvalue
    ),
  sheet("module_contrasts", "fry_tests") |>
    transmute(
      feature = module, contrast, test = paste("fry", contrast),
      effect = if_else(direction == "Up", 1, -1), p, fdr, qvalue
    )
) |>
  mutate(family = "contrast")

evidence <- list(
  proteins = bind_rows(
    protein_de, ancova_rows("protein_contrasts"), classify("protein_classify"),
    associations("protein_associate")
  ),
  sets = bind_rows(
    set_de, ancova_rows("set_contrasts"), classify("set_classify"), associations("set_associate")
  ),
  modules = bind_rows(
    module_de, ancova_rows("module_contrasts"), classify("module_classify"),
    associations("module_associate")
  )
) |>
  list_rbind(names_to = "level") |>
  filter(!is.na(p))

on_list <- on_lists(evidence)

protein_names <- sheet("modules", "membership") |>
  transmute(feature = uniprot_id, label = gene, module)
set_names_tbl <- readRDS(paths[["gene_sets"]])$set_catalog |>
  filter(qualifies) |>
  transmute(feature = set_id, label = pathway, collection)
labels <- bind_rows(
  mutate(protein_names, level = "proteins"),
  mutate(set_names_tbl, level = "sets"),
  tibble(level = "modules", feature = unique(filter(evidence, level == "modules")$feature)) |>
    mutate(label = feature)
)
baseline_de <- evidence |>
  filter(contrast == "Baseline_HRvLR", test %in% c(
    "Baseline_HRvLR", "fgsea Baseline_HRvLR",
    "eigengene Baseline_HRvLR"
  )) |>
  select(level, feature, baseline_effect = effect, baseline_p = p)
unadjusted <- sheet("protein_contrasts", "unadjusted") |>
  transmute(
    level = "proteins", feature = uniprot_id,
    unadjusted_p_training = Training_Interaction_p, unadjusted_p_acute = Acute_Interaction_p
  )

# A contrast hit is stable when it stays at p < 0.05 with any one subject left out.
stability <- list(
  proteins = "protein_contrasts", sets = "set_contrasts", modules = "module_contrasts"
) |>
  map(\(input) sheet(input, "leave_one_out")) |>
  list_rbind(names_to = "level") |>
  transmute(level, feature, test, stable = n_nominal == n_reruns)
hits <- filter(on_list, p < 0.05) |>
  left_join(stability, by = c("level", "feature", "test"))
candidates <- hits |>
  arrange(p) |>
  summarise(
    n_hits = n(), n_contrast = sum(family == "contrast"), n_classify = sum(family == "classify"),
    n_associate = sum(family == "associate"), n_stable = sum(stable, na.rm = TRUE),
    n_bh_20 = sum(fdr < 0.2, na.rm = TRUE), n_storey_20 = sum(qvalue < 0.2, na.rm = TRUE),
    min_p = min(p), min_bh = min(fdr, na.rm = TRUE), min_storey = min(qvalue, na.rm = TRUE),
    hits = paste(
      if_else(
        is.na(fdr), sprintf("%s (%+.3g, p %.2g)", test, effect, p),
        sprintf("%s (%+.3g, p %.2g, BH %.2g, q %.2g)", test, effect, p, fdr, qvalue)
      ),
      collapse = "; "
    ),
    .by = c(level, list, feature)
  ) |>
  mutate(on_both = all(c("training", "acute") %in% list), .by = c(level, feature)) |>
  left_join(labels, by = c("level", "feature")) |>
  left_join(baseline_de, by = c("level", "feature")) |>
  left_join(unadjusted, by = c("level", "feature")) |>
  left_join(
    select(sheet("permutation", "feature_scores"), level, list, feature, permutation_score),
    by = c("level", "list", "feature")
  ) |>
  arrange(level, list, desc(n_hits), min_p, feature)
stopifnot(!anyNA(candidates$label[candidates$level != "modules"]))

tested <- on_list |>
  summarise(n_tests = n(), .by = c(level, list, feature)) |>
  summarise(
    features = n(), median_tests = median(n_tests), max_tests = max(n_tests),
    .by = c(level, list)
  )
distribution <- on_list |>
  summarise(n_hits = sum(p < 0.05), .by = c(level, list, feature)) |>
  mutate(n_hits = pmin(n_hits, 6)) |>
  count(level, list, n_hits) |>
  mutate(n_hits = if_else(n_hits == 6, "6+", as.character(n_hits))) |>
  pivot_wider(names_from = n_hits, values_from = n, values_fill = 0, names_prefix = "hits_") |>
  select(level, list, num_range("hits_", 0:5), any_of("hits_6+")) |>
  left_join(tested, by = c("level", "list"))

table_of <- function(level, list) {
  d <- filter(candidates, .data$level == .env$level, .data$list == .env$list)
  cols <- c(
    "feature", "label", if (level == "proteins") "module", if (level == "sets") "collection",
    "n_hits", "permutation_score", "n_contrast", "n_stable", "n_classify", "n_associate",
    "n_bh_20", "n_storey_20",
    "min_p", "min_bh", "min_storey", "on_both",
    "baseline_effect",
    "baseline_p", if (level == "proteins") paste0("unadjusted_p_", list), "hits"
  )
  select(d, all_of(cols))
}
tables <- expand_grid(level = c("proteins", "sets", "modules"), list = c("training", "acute"))
sheets <- set_names(
  pmap(tables, table_of),
  paste(tables$level, tables$list, sep = "_")
)
holds <- sprintf(
  "%s with at least one p < 0.05 on the %s list, by hit count then smallest p.",
  stringr::str_to_sentence(tables$level), tables$list
)

sheets <- c(sheets, list(
  hit_distribution = distribution,
  evidence = hits |>
    select(level, list, feature, family, test, effect, p, fdr, qvalue, stable) |>
    arrange(level, list, feature, p)
))
holds <- c(
  holds,
  paste(
    "Features by hit count per level and list, with tests per feature. The baseline list applies",
    "the same rules at T1, before training, as the reference."
  ),
  paste(
    "Every p < 0.05 behind the tables, one row per test; stable marks contrast hits that hold",
    "with any one subject left out."
  )
)
write_workbook(sheets, holds, file.path(out, "02_candidates.xlsx"), inputs, paths,
  note = paste(
    "Tests are counted as they are and overlap: the interaction, the training-change",
    "classification and the comp_hypertrophy association all read each subject's training change,",
    "and arm is the median split of comp_hypertrophy. Effects: logFC for protein contrasts, NES",
    "for fgsea, +1 or -1 for fry direction, eigengene logFC, ANCOVA arm effect, AUC, Spearman rho,",
    "slope difference or model coefficient. Contrast rows come from the blood-adjusted fit;",
    "BH and q are exploratory corrections within each test's own screen (Storey's q with a",
    "bootstrap pi0, equal to BH in screens under 100 tests); n_bh_20 and n_storey_20 count hits",
    "below 0.20.",
    "unadjusted_p is the protein interaction p before adjustment. n_stable counts contrast hits",
    "that stay at p < 0.05 in every leave-one-subject-out rerun."
  )
)
