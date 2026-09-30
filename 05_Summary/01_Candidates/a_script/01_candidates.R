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

out <- here("05_Summary", "01_Candidates", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  protein_contrasts = "02_Differential_Expression/02_Contrasts/c_data/02_contrasts.xlsx",
  protein_classify = "02_Differential_Expression/03_Classify/c_data/03_classify.xlsx",
  protein_associate = "02_Differential_Expression/04_Associate/c_data/04_associate.xlsx",
  set_tests = "03_Pathway_Enrichment/02_Contrasts/c_data/set_tests.rds",
  set_contrasts = "03_Pathway_Enrichment/02_Contrasts/c_data/02_contrasts.xlsx",
  set_classify = "03_Pathway_Enrichment/03_Classify/c_data/03_classify.xlsx",
  set_associate = "03_Pathway_Enrichment/04_Associate/c_data/04_associate.xlsx",
  gene_sets = "03_Pathway_Enrichment/00_Gene_Sets/c_data/gene_sets.rds",
  modules = "04_Network/01_Build_Modules/c_data/01_build_modules.xlsx",
  module_contrasts = "04_Network/02_Contrasts/c_data/02_contrasts.xlsx",
  module_classify = "04_Network/03_Classify/c_data/03_classify.xlsx",
  module_associate = "04_Network/04_Associate/c_data/04_associate.xlsx"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
sheet <- \(input, name) suppressMessages(read_excel(paths[[input]], name))

# Which tests feed which list. Acute_Interaction feeds both, by decision (docs/decisions.md).
lists <- list(
  training = list(
    contrasts = c("Training_Interaction", "Trained_HRvLR", "Acute_Interaction"),
    ancova = "training", tasks = c("Training_change_HRvLR", "Trained_HRvLR"), window = "training"
  ),
  acute = list(
    contrasts = c("Acute_Interaction", "Acute_HRvLR"),
    ancova = "acute", tasks = "Acute_change_HRvLR", window = "acute"
  ),
  baseline = list(
    contrasts = "Baseline_HRvLR", ancova = character(0), tasks = "Baseline_HRvLR",
    window = "baseline"
  )
)

# One long table per level: feature, family, test, effect, p, and the window or contrast that
# places it on a list.
associations <- function(input) {
  bind_rows(
    sheet(input, "change_score") |>
      transmute(feature, window, test = paste("pooled rho", outcome), effect = rho, p),
    sheet(input, "change_score_by_arm") |>
      transmute(feature, window, test = paste(arm, "rho", outcome), effect = rho, p),
    sheet(input, "arm_difference") |>
      transmute(feature, window, test = paste("slope HR-LR", outcome), effect, p),
    sheet(input, "two_predictor") |>
      transmute(
        feature, window,
        test = paste(coefficient, "term,", sub("pct_", "", strength_test), "model"), effect, p
      )
  ) |>
    mutate(family = "associate")
}
classify <- function(input) {
  sheet(input, "auc") |>
    filter(!paired) |>
    transmute(feature, task, test = paste("classify", task), effect = auc, p, family = "classify")
}
ancova <- function(input) {
  sheet(input, "ancova") |>
    transmute(feature,
      window = paste0("ancova_", window), test = paste("ANCOVA", window), effect, p,
      family = "contrast"
    )
}

protein_de <- sheet("protein_contrasts", "DEP_matrix") |>
  select(feature = uniprot_id, ends_with("_logFC"), ends_with("_p")) |>
  pivot_longer(-feature, names_to = c("contrast", ".value"), names_pattern = "(.*)_(logFC|p)$") |>
  transmute(feature, contrast, test = contrast, effect = logFC, p, family = "contrast")
set_de <- readRDS(paths[["set_tests"]])$set_tests |>
  transmute(
    feature = set_id, contrast, test = paste(method, contrast),
    effect = if_else(method == "fgsea", nes, if_else(direction == "Up", 1, -1)), p,
    family = "contrast"
  )
module_de <- bind_rows(
  sheet("module_contrasts", "eigengene_tests") |>
    transmute(feature = module, contrast, test = paste("eigengene", contrast), effect = logFC, p),
  sheet("module_contrasts", "fry_tests") |>
    transmute(
      feature = module, contrast, test = paste("fry", contrast),
      effect = if_else(direction == "Up", 1, -1), p
    )
) |>
  mutate(family = "contrast")

evidence <- list(
  proteins = bind_rows(
    protein_de, ancova("protein_contrasts"), classify("protein_classify"),
    associations("protein_associate")
  ),
  sets = bind_rows(
    set_de, ancova("set_contrasts"), classify("set_classify"), associations("set_associate")
  ),
  modules = bind_rows(
    module_de, ancova("module_contrasts"), classify("module_classify"),
    associations("module_associate")
  )
) |>
  list_rbind(names_to = "level") |>
  filter(!is.na(p))

on_list <- imap(lists, \(rule, name) {
  evidence |>
    filter(
      contrast %in% rule$contrasts | window %in% c(paste0("ancova_", rule$ancova), rule$window) |
        task %in% rule$tasks
    ) |>
    mutate(list = name)
}) |>
  list_rbind()

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

hits <- filter(on_list, p < 0.05)
candidates <- hits |>
  arrange(p) |>
  summarise(
    n_hits = n(), n_contrast = sum(family == "contrast"), n_classify = sum(family == "classify"),
    n_associate = sum(family == "associate"), min_p = min(p),
    hits = paste(sprintf("%s (%+.3g, p %.2g)", test, effect, p), collapse = "; "),
    .by = c(level, list, feature)
  ) |>
  mutate(on_both = all(c("training", "acute") %in% list), .by = c(level, feature)) |>
  left_join(labels, by = c("level", "feature")) |>
  left_join(baseline_de, by = c("level", "feature")) |>
  left_join(unadjusted, by = c("level", "feature")) |>
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
    "n_hits", "n_contrast", "n_classify", "n_associate", "min_p", "on_both", "baseline_effect",
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
  evidence = filter(on_list, p < 0.05) |>
    select(level, list, feature, family, test, effect, p) |>
    arrange(level, list, feature, p),
  input_manifest = tibble(
    input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
  ),
  package_versions = sessioninfo::package_info("loaded", dependencies = FALSE) |>
    as_tibble() |>
    select(package, version = loadedversion, source)
))
holds <- c(
  holds,
  paste(
    "Features by hit count per level and list, with tests per feature. The baseline list applies",
    "the same rules at T1, before training, as the reference."
  ),
  "Every p < 0.05 behind the tables, one row per test.",
  "Files read, with md5.",
  "Packages loaded at run time."
)
read_me <- tibble(sheet = names(sheets), holds = holds) |>
  add_row(sheet = "note", holds = paste(
    "Tests are counted as they are and overlap: the interaction, the training-change",
    "classification and the comp_hypertrophy association all read each subject's training change,",
    "and arm is the median split of comp_hypertrophy. Effects: logFC for protein contrasts, NES",
    "for fgsea, +1 or -1 for fry direction, eigengene logFC, ANCOVA arm effect, AUC, Spearman rho,",
    "slope difference or model coefficient. Contrast rows come from the blood-adjusted fit;",
    "unadjusted_p is the protein interaction p before adjustment."
  ))
write_xlsx(c(list(read_me = read_me), sheets), file.path(out, "01_candidates.xlsx"))
