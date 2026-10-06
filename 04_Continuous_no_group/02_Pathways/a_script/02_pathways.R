# Cross-level continuous analyses, every subject pooled and each trait continuous. The per-level
# associations live in each 04_Associate; this step adds what needs more than one level:
#   pathways     fgsea and cameraPR on the protein ranking by each trait (the limma t saved by
#                04_Continuous_no_group/01_Associate/01_Proteins), against 100 phenotype shuffles: sets at FDR < 0.05
#                per screen, and per set the share of shuffles with as large an |NES|
#   spls         sparse PLS of each level against all traits of a pairing at once (first-variate
#                correlation), against 200 shuffles of the trait rows
#   membership   per module, the correlation of protein module membership (kME) with |t|: does
#                the module's core carry the trait link? Against the same 100 shuffles
#   string       STRINGdb interactions among each screen's top 50 proteins, against the top 50 of
#                20 shuffles
#   composition  the growth x strength interaction refit with the connective-tissue (magenta)
#                and non-fibre (yellow) module changes as covariates, per level
# The shuffled protein t is saved for 05_Clustering.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
  library(STRINGdb)
})
source(here("R", "helpers.R"))

stage <- here("04_Continuous_no_group", "02_Pathways")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  corrected = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  associate = "04_Continuous_no_group/01_Associate/01_Proteins/c_data/associate.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds",
  singscore = "01_Trunk/05_Scores/c_data/singscore.rds",
  modules = "01_Trunk/06_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv",
  string_links = "00_Input/downloads/9606.protein.links.v12.0.min700.txt.gz",
  string_aliases = "00_Input/downloads/9606.protein.aliases.v12.0.txt.gz",
  string_info = "00_Input/downloads/9606.protein.info.v12.0.txt.gz"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
n_shuffles <- 100
cores <- min(8, parallel::detectCores())

dal <- readRDS(paths[["corrected"]])
meta <- as_tibble(dal$metadata)
gs <- readRDS(paths[["gene_sets"]])
catalog <- filter(gs$set_catalog, qualifies)
modules <- readRDS(paths[["modules"]])
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
pairings <- phenotype_tables(phenotype)$pairings
feature_label <- c(
  with(dal$annotation, set_names(gene, uniprot_id)),
  with(catalog, set_names(pathway, set_id)),
  set_names(rownames(modules$eigengenes))
)
inputs_of <- function(values) {
  list(
    training = subject_change(values, meta, "T1", "T2"),
    baseline = subject_values(values, meta, "T1"),
    acute = subject_change(values, meta, "T2", "T3"),
    baseline_level = subject_values(values, meta, "T1"),
    post_level = subject_values(values, meta, "T2")
  )
}
levels_of <- list(
  proteins = inputs_of(dal$data),
  sets = inputs_of(readRDS(paths[["singscore"]])$scores[catalog$set_id, ]),
  modules = inputs_of(modules$eigengenes)
)
screens <- imap(pairings, \(tt, window) tibble(window, outcome = setdiff(names(tt), "subject"))) |>
  list_rbind()
perm <- \(null_values, observed) (sum(null_values >= observed) + 1) / (length(null_values) + 1)

# Observed protein t from 01_Associate/01_Proteins; the shuffles refit the same limma model.
observed_t <- readRDS(paths[["associate"]])$change_score |>
  select(window, outcome, feature, t, p)
protein_t <- function(tables) {
  pmap(screens, \(window, outcome) {
    x <- levels_of$proteins[[window]]
    subject_fit(
      x, data.frame(trait = pairing_value(tables, window, outcome, colnames(x))),
      ~trait, "trait"
    ) |>
      transmute(window, outcome, feature, t, p)
  }) |>
    list_rbind()
}
shuffled_t <- parallel::mclapply(seq_len(n_shuffles), \(i) {
  protein_t(shuffle_pairings(pairings, i))
}, mc.cores = cores)
saveRDS(shuffled_t, file.path(out, "shuffled_t.rds"), compress = "xz")

# Pathways from the protein ranking.
selected <- filter(gs$protein_map, selected)
gene_of <- set_names(selected$gene, selected$uniprot_id)
rank_tests <- function(d) {
  d <- filter(d, feature %in% selected$uniprot_id)
  stats <- set_names(d$t, gene_of[d$feature])
  set.seed(1)
  fg <- fgsea::fgsea(gs$sets, stats, minSize = 15, maxSize = 500) |>
    as_tibble() |>
    transmute(set_id = pathway, size, nes = NES, p = pval, fdr = padj)
  index <- map(gs$sets, \(genes) which(names(stats) %in% genes))
  camera <- cameraPR(stats, index[lengths(index) >= 15]) |>
    as_tibble(rownames = "set_id") |>
    transmute(set_id, camera_direction = Direction, p_camera = PValue)
  left_join(fg, camera, by = "set_id") |>
    mutate(
      direction = if_else(camera_direction == "Up", 1, -1),
      camera_agrees = fdr < 0.05 & p_camera < 0.05 & sign(nes) == direction
    ) |>
    select(-direction)
}
by_screen <- \(d) split(d, ~ window + outcome, drop = TRUE, sep = "|")
pathways <- by_screen(observed_t) |>
  imap(\(d, key) mutate(rank_tests(d), screen = key)) |>
  list_rbind()
pathway_null <- parallel::mclapply(shuffled_t, \(t_long) {
  imap(by_screen(t_long), \(d, key) {
    rank_tests(d) |> transmute(screen = key, set_id, nes, significant = fdr < 0.05)
  }) |>
    list_rbind()
}, mc.cores = cores) |>
  list_rbind(names_to = "shuffle")
pathways <- pathways |>
  left_join(summarise(pathway_null, abs_null = list(abs(nes)), .by = c(screen, set_id)),
    by = c("screen", "set_id")
  ) |>
  mutate(set_score = map2_dbl(abs_null, abs(nes), perm), label = unname(feature_label[set_id])) |>
  select(-abs_null) |>
  separate_wider_delim(screen, "|", names = c("window", "outcome"), cols_remove = FALSE) |>
  arrange(window, outcome, p)
pathway_chance <- pathways |>
  summarise(
    sets_fdr_05 = sum(fdr < 0.05, na.rm = TRUE),
    camera_agreed = sum(camera_agrees, na.rm = TRUE), .by = c(window, outcome, screen)
  ) |>
  left_join(summarise(pathway_null, n = sum(significant, na.rm = TRUE), .by = c(shuffle, screen)) |>
    summarise(
      null_median = median(n), null_95 = quantile(n, 0.95), null_n = list(n),
      .by = screen
    ), by = "screen") |>
  mutate(perm_p = map2_dbl(null_n, sets_fdr_05, perm)) |>
  select(-null_n)

# Sparse PLS per level and pairing.
spls_r <- function(x, y) {
  fit <- mixOmics::spls(x, y, ncomp = 2, keepX = c(10, 10), mode = "regression")
  list(r = abs(cor(fit$variates$X[, 1], fit$variates$Y[, 1])), fit = fit)
}
spls_table <- imap(levels_of, \(by_window, level) {
  imap(by_window, \(x, window) {
    tt <- pairings[[window]]
    y <- as.matrix(select(tt, -subject))
    rownames(y) <- tt$subject
    subjects <- intersect(colnames(x), rownames(y)[complete.cases(y)])
    xm <- t(x[, subjects, drop = FALSE])
    xm <- xm[, apply(xm, 2, sd) > 0, drop = FALSE]
    ym <- y[subjects, , drop = FALSE]
    observed <- spls_r(xm, ym)
    null_r <- map_dbl(seq_len(200), \(i) {
      set.seed(i)
      spls_r(xm, ym[sample(nrow(ym)), , drop = FALSE])$r
    })
    loadings_x <- observed$fit$loadings$X[, 1]
    loadings_y <- observed$fit$loadings$Y[, 1]
    tibble(level, window,
      n_subjects = length(subjects), n_features = ncol(xm), r = observed$r,
      null_median = median(null_r), null_95 = quantile(null_r, 0.95),
      perm_p = perm(null_r, observed$r),
      features = paste(
        feature_label[names(sort(abs(loadings_x[loadings_x != 0]), decreasing = TRUE))],
        collapse = ", "
      ),
      trait_loadings = paste(sprintf("%s %+.2f", names(loadings_y), loadings_y), collapse = "; ")
    )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(perm_p, "BH"), screen = paste(level, window))

# Module membership against the trait link.
kme <- modules$membership |>
  filter(module != "grey") |>
  select(feature = uniprot_id, module, kme)
membership_of <- function(t_long) {
  inner_join(t_long, kme, by = "feature", relationship = "many-to-many") |>
    summarise(r = cor(kme, abs(t)), n_members = n(), .by = c(window, outcome, module))
}
membership <- membership_of(observed_t) |>
  left_join(
    list_rbind(map(shuffled_t, membership_of)) |>
      summarise(null_r = list(r), .by = c(window, outcome, module)),
    by = c("window", "outcome", "module")
  ) |>
  mutate(perm_p = map2_dbl(null_r, r, perm), fdr = p.adjust(perm_p, "BH")) |>
  select(-null_r) |>
  arrange(perm_p)

# STRING density of each screen's top 50 proteins.
string_db <- STRINGdb$new(
  version = "12.0", species = 9606, score_threshold = 700L,
  network_type = "full", input_directory = here("00_Input", "downloads")
)
mapped <- string_db$map(data.frame(uniprot_id = rownames(dal$data)), "uniprot_id",
  removeUnmappedRows = TRUE
)
string_db$set_background(mapped$STRING_id)
string_of <- set_names(mapped$STRING_id, mapped$uniprot_id)
edges_of <- \(features) string_db$get_ppi_enrichment(unique(na.omit(string_of[features])))$edges
top50 <- \(t_long) slice_min(t_long, p, n = 50, by = c(window, outcome), with_ties = FALSE)
string_table <- top50(observed_t) |>
  summarise(edges = edges_of(feature), .by = c(window, outcome)) |>
  left_join(map(shuffled_t[1:20], \(t_long) {
    top50(t_long) |> summarise(edges = edges_of(feature), .by = c(window, outcome))
  }) |>
    list_rbind() |>
    summarise(
      null_median = median(edges), null_95 = quantile(edges, 0.95), null_e = list(edges),
      .by = c(window, outcome)
    ), by = c("window", "outcome")) |>
  mutate(perm_p = map2_dbl(null_e, edges, perm), screen = paste(window, outcome)) |>
  select(-null_e)

# Composition: the interaction with and without the magenta and yellow module changes.
composition_covariates <- levels_of$modules[c("training", "baseline", "acute")] |>
  map(\(x) t(x[c("magenta", "yellow"), , drop = FALSE]))
composition <- imap(levels_of, \(by_window, level) {
  imap(by_window[c("training", "baseline", "acute")], \(x, window) {
    map(c("pct_leg_press_1rm", "pct_leg_ext_1rm"), \(test) {
      predictors <- data.frame(
        growth = pairing_value(pairings, window, "comp_hypertrophy", colnames(x)),
        strength = pairing_value(pairings, window, test, colnames(x)),
        composition_covariates[[window]][colnames(x), , drop = FALSE]
      )
      plain <- subject_fit(x, predictors, ~ growth * strength, "growth:strength")
      adjusted <- subject_fit(
        x, predictors, ~ growth * strength + magenta + yellow,
        "growth:strength"
      )
      left_join(select(plain, feature, effect, p),
        select(adjusted, feature, effect_adjusted = effect, p_adjusted = p),
        by = "feature"
      ) |>
        mutate(level, window, strength_test = test, .before = 1)
    }) |>
      list_rbind()
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(
    fdr = p.adjust(p, "BH"), fdr_adjusted = p.adjust(p_adjusted, "BH"),
    .by = c(level, window, strength_test)
  ) |>
  mutate(label = unname(feature_label[feature])) |>
  arrange(level, window, strength_test, p)

write_workbook(
  list(
    pathway_chance = pathway_chance, pathways = pathways, spls = spls_table,
    membership = membership, string = string_table, composition = composition
  ),
  c(
    "Per screen: sets at fgsea FDR < 0.05 and those cameraPR corroborates, against 100 shuffles.",
    paste(
      "Per screen and set: fgsea NES, p, FDR; cameraPR direction and p; set_score, the share",
      "of 100 shuffles with |NES| as large."
    ),
    "Per level and pairing: sPLS first-variate correlation against 200 trait-row shuffles.",
    "Per module and screen: cor(kME, |t|) among members against 100 shuffles.",
    "Per screen: STRING edges among the top 50 proteins against the top 50 of 20 shuffles.",
    paste(
      "Per feature: growth x strength interaction effect and p, then the same with the",
      "magenta and yellow module changes as covariates, BH per level, window and strength test."
    )
  ),
  file.path(out, "04_continuous.xlsx"),
  inputs, paths
)

pdf(file.path(figures, "04_continuous_figures.pdf"), width = 11, height = 8.5, useDingbats = TRUE)
print(plot_chance(
  pathway_chance, "sets_fdr_05",
  "Sets at fgsea FDR < 0.05 from each trait's protein ranking, against 100 phenotype shuffles",
  paste(
    "Correlated proteins make set tests fire on any ranking, so the comparison is with the",
    "grey bar (shuffles' 95th percentile; tick, median), not with zero. Table:",
    "c_data/04_continuous.xlsx, pathway_chance."
  )
))
print(plot_chance(spls_table, "r", "Sparse PLS: first-variate correlation against 200 shuffles",
  "Each row a level and pairing; all traits of the pairing at once. Table: spls.",
  facet = ~level
))
print(plot_chance(
  string_table, "edges",
  "STRING interactions among each screen's top 50 proteins, against 20 shuffles",
  "Background: the measured proteins. Table: string."
))
print(
  ggplot(composition, aes(-log10(p), -log10(p_adjusted))) +
    geom_abline(linetype = 2, colour = "grey60") +
    geom_point(alpha = 0.3, size = 0.6) +
    facet_grid(level ~ paste(window, strength_test)) +
    labs(
      x = "-log10 p, growth x strength", y = "-log10 p, adjusted for composition",
      title = "Composition check: the interaction with and without the magenta and yellow changes",
      caption = "Points on the diagonal do not depend on the composition covariates."
    ) +
    page_theme
)
invisible(dev.off())
