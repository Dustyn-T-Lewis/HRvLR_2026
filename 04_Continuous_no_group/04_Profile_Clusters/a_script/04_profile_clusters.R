# Clustering, with every structure judged against what shuffled data give.
#   profiles     proteins clustered by their association profile, the limma t against every trait
#                in every pairing (04_Continuous_no_group/01_Associate/01_Proteins); average linkage on 1 - Pearson r,
#                cut by dynamicTreeCut; mean silhouette against profiles from 100 phenotype shuffles
#                (04_Continuous_no_group/02_Pathways)
#   trait        per trait, the proteins at p < 0.05 grouped by co-variation across people, each
#                group annotated (direction, module, top pathway, STRING); per trait family, each
#                protein's t across the pairings clustered by PAM (k 2 to 6)
#   people       the subjects clustered by their proteome, HR/LR unused: PAM, k 2 to 4, on module
#                eigengenes and on the 500 most variable proteins, per input; silhouette against
#                200 shuffles of each feature across people; agreement with HR/LR (adjusted Rand
#                index, Fisher) and with each trait (Kruskal-Wallis)
#   groups       the subjects clustered within each pathway's and module's proteins (average
#                linkage, k = 2), per input; the best silhouette and the best agreement with HR/LR
#                judged against the best over 100 shuffled searches
#   blocks       two-way hierarchical blocks of the profile matrix and biclust's plaid model

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(writexl)
  library(cluster)
  library(STRINGdb)
})
source(here("R", "helpers.R"))

stage <- here("04_Continuous_no_group", "04_Profile_Clusters")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
for (path in c(out, figures)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  corrected = "01_Trunk/03_Quantify/c_data/DAList_corrected.rds",
  associate = "04_Continuous_no_group/01_Associate/01_Proteins/c_data/associate.rds",
  gene_sets = "01_Trunk/04_Gene_Sets/c_data/gene_sets.rds",
  modules = "01_Trunk/06_Modules/c_data/modules.rds",
  shuffled_t = "04_Continuous_no_group/02_Pathways/c_data/shuffled_t.rds",
  phenotype = "00_Input/phenotype.csv",
  string_links = "00_Input/downloads/9606.protein.links.v12.0.min700.txt.gz",
  string_aliases = "00_Input/downloads/9606.protein.aliases.v12.0.txt.gz",
  string_info = "00_Input/downloads/9606.protein.info.v12.0.txt.gz"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))
cores <- min(8, parallel::detectCores())

dal <- readRDS(paths[["corrected"]])
meta <- as_tibble(dal$metadata)
gs <- readRDS(paths[["gene_sets"]])
modules <- readRDS(paths[["modules"]])
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
pairings <- phenotype_tables(phenotype)$pairings
arm_of <- set_names(phenotype$arm, phenotype$subject)
gene_of <- with(filter(gs$protein_map, selected), set_names(gene, uniprot_id))
feature_label <- c(
  with(dal$annotation, set_names(gene, uniprot_id)),
  with(gs$set_catalog, set_names(pathway, set_id)),
  set_names(rownames(modules$eigengenes))
)
module_of <- with(modules$membership, set_names(module, uniprot_id))
observed_t <- readRDS(paths[["associate"]])$change_score |> select(window, outcome, feature, t, p)
shuffled_t <- readRDS(paths[["shuffled_t"]])
perm <- \(null_values, observed) (sum(null_values >= observed) + 1) / (length(null_values) + 1)
mean_silhouette <- \(groups, d) {
  if (length(unique(groups)) < 2) {
    return(NA_real_)
  }
  mean(silhouette(groups, d)[, "sil_width"])
}
inputs_of <- function(values) {
  list(
    T1 = subject_values(values, meta, "T1"), training = subject_change(values, meta, "T1", "T2"),
    T2 = subject_values(values, meta, "T2"), acute = subject_change(values, meta, "T2", "T3"),
    T3 = subject_values(values, meta, "T3")
  )
}
protein_inputs <- inputs_of(dal$data)

term2gene <- tibble(set_id = rep(names(gs$sets), lengths(gs$sets)), gene = unlist(gs$sets))
top_set <- function(members, universe) {
  res <- clusterProfiler::enricher(na.omit(unname(gene_of[members])),
    universe = na.omit(unname(gene_of[universe])), TERM2GENE = term2gene, minGSSize = 15,
    maxGSSize = 500, pvalueCutoff = 1, qvalueCutoff = 1
  )
  if (is.null(res) || nrow(res@result) == 0) {
    return(tibble(top_set = NA_character_, top_set_fdr = NA_real_))
  }
  top <- slice_min(as_tibble(res@result), pvalue, n = 1, with_ties = FALSE)
  tibble(top_set = unname(feature_label[top$ID]), top_set_fdr = top$p.adjust)
}
string_db <- STRINGdb$new(
  version = "12.0", species = 9606, score_threshold = 700L, network_type = "full",
  input_directory = here("00_Input", "downloads")
)
mapped <- string_db$map(data.frame(uniprot_id = rownames(dal$data)), "uniprot_id",
  removeUnmappedRows = TRUE
)
string_db$set_background(mapped$STRING_id)
string_of <- set_names(mapped$STRING_id, mapped$uniprot_id)
string_density <- function(members) {
  ids <- unique(na.omit(string_of[members]))
  if (length(ids) < 5) {
    return(tibble(edges = NA_real_, expected = NA_real_, p_string = NA_real_))
  }
  e <- string_db$get_ppi_enrichment(ids)
  tibble(edges = e$edges, expected = e$lambda, p_string = e$enrichment)
}

# Profiles.
profile_matrix <- function(t_long) {
  t_long |>
    transmute(feature, screen = paste(window, outcome), t) |>
    pivot_wider(names_from = screen, values_from = t) |>
    column_to_rownames("feature") |>
    as.matrix() |>
    na.omit()
}
cut_profiles <- function(m) {
  d <- as.dist(1 - cor(t(m)))
  groups <- dynamicTreeCut::cutreeDynamic(hclust(d, "average"),
    distM = as.matrix(d),
    minClusterSize = 20, deepSplit = 2, verbose = 0
  )
  kept <- groups != 0
  list(
    groups = set_names(groups, rownames(m)), n_groups = n_distinct(groups[kept]),
    silhouette = mean_silhouette(groups[kept], as.dist(as.matrix(d)[kept, kept]))
  )
}
profiles <- profile_matrix(observed_t)
profile_cut <- cut_profiles(profiles)
profile_null <- parallel::mclapply(shuffled_t, \(t_long) {
  cut_profiles(profile_matrix(t_long))$silhouette
}, mc.cores = cores) |>
  unlist()
profile_groups <- tibble(feature = names(profile_cut$groups), group = profile_cut$groups) |>
  filter(group != 0) |>
  mutate(label = unname(feature_label[feature]), module = unname(module_of[feature]))
profile_means <- as_tibble(profiles[profile_groups$feature, ], rownames = "feature") |>
  left_join(select(profile_groups, feature, group), by = "feature") |>
  pivot_longer(-c(feature, group), names_to = "screen", values_to = "t") |>
  summarise(mean_t = mean(t), .by = c(group, screen))
profile_summary <- profile_groups |>
  summarise(
    n_proteins = n(), members = list(feature),
    main_module = names(which.max(table(module))), main_module_share = max(table(module)) / n(),
    .by = group
  ) |>
  mutate(enrichment = map(members, \(m) top_set(m, rownames(profiles)))) |>
  unnest(enrichment) |>
  select(-members) |>
  arrange(group)

# Per trait: co-variation groups among the associated proteins, and profiles across pairings.
covariation_groups <- function(features, x) {
  m <- x[intersect(features, rownames(x)), , drop = FALSE]
  d <- as.dist(1 - cor(t(m)))
  groups <- dynamicTreeCut::cutreeDynamic(hclust(d, "average"),
    distM = as.matrix(d),
    minClusterSize = 8, deepSplit = 2, verbose = 0
  )
  kept <- groups != 0
  list(
    groups = set_names(groups, rownames(m)),
    silhouette = mean_silhouette(groups[kept], as.dist(as.matrix(d)[kept, kept]))
  )
}
pairing_input <- c(training = "training", baseline_level = "T1")
trait_screens <- filter(distinct(observed_t, window, outcome), window %in% names(pairing_input))
trait_groups <- pmap(trait_screens, \(window, outcome) {
  d <- filter(observed_t, .data$window == .env$window, .data$outcome == .env$outcome)
  hits <- filter(d, p < 0.05)
  x <- protein_inputs[[pairing_input[[window]]]]
  observed <- covariation_groups(hits$feature, x)
  null <- map_dbl(shuffled_t[1:50], \(s) {
    top <- filter(s, .data$window == .env$window, .data$outcome == .env$outcome) |>
      slice_min(p, n = nrow(hits), with_ties = FALSE)
    covariation_groups(top$feature, x)$silhouette
  })
  annotated <- tibble(feature = names(observed$groups), group = observed$groups) |>
    filter(group != 0) |>
    left_join(select(hits, feature, t), by = "feature") |>
    summarise(
      n_proteins = n(), mean_t = mean(t), share_positive = mean(t > 0),
      members = list(feature), main_module = names(which.max(table(module_of[feature]))),
      proteins = paste(head(feature_label[feature][order(-abs(t))], 15), collapse = ", "),
      .by = group
    ) |>
    mutate(
      enrichment = map(members, \(m) top_set(m, d$feature)),
      string = map(members, string_density)
    ) |>
    unnest(c(enrichment, string)) |>
    select(-members) |>
    mutate(window = window, outcome = outcome, .before = 1)
  list(
    summary = tibble(window, outcome,
      proteins = nrow(hits), groups = nrow(annotated),
      silhouette = observed$silhouette, null_median = median(null, na.rm = TRUE),
      null_95 = quantile(null, 0.95, na.rm = TRUE),
      perm_p = perm(na.omit(null), observed$silhouette)
    ),
    groups = annotated
  )
})
family_of <- \(outcome) sub("^pct_", "", outcome)
families <- c("fcsa_mixed", "fcsa_type1", "fcsa_type2", "mcsa", "leg_press_1rm", "leg_ext_1rm")
family_profiles <- \(t_long, family) {
  t_long |>
    filter(family_of(outcome) == family) |>
    select(feature, window, t) |>
    pivot_wider(names_from = window, values_from = t) |>
    column_to_rownames("feature") |>
    as.matrix() |>
    na.omit()
}
best_pam <- function(m) {
  d <- dist(m)
  fits <- map(2:6, \(k) pam(d, k, diss = TRUE))
  best <- which.max(map_dbl(fits, \(f) f$silinfo$avg.width))
  list(
    k = best + 1, silhouette = fits[[best]]$silinfo$avg.width,
    groups = set_names(fits[[best]]$clustering, rownames(m))
  )
}
trait_profiles <- map(set_names(families), \(family) {
  m <- family_profiles(observed_t, family)
  observed <- best_pam(m)
  null <- unlist(parallel::mclapply(shuffled_t,
    \(s) best_pam(family_profiles(s, family))$silhouette,
    mc.cores = cores
  ))
  means <- as_tibble(m, rownames = "feature") |>
    mutate(group = observed$groups[feature]) |>
    pivot_longer(-c(feature, group), names_to = "window", values_to = "t") |>
    summarise(mean_t = mean(t), n_proteins = n_distinct(feature), .by = c(group, window)) |>
    mutate(family = family, .before = 1)
  annotation <- tibble(feature = names(observed$groups), group = observed$groups) |>
    summarise(
      n_proteins = n(), members = list(feature),
      main_module = names(which.max(table(module_of[feature]))), .by = group
    ) |>
    mutate(enrichment = map(members, \(mm) top_set(mm, rownames(m)))) |>
    unnest(enrichment) |>
    select(-members) |>
    mutate(family = family, .before = 1)
  list(
    summary = tibble(family,
      proteins = nrow(m), k = observed$k, silhouette = observed$silhouette,
      null_median = median(null), null_95 = quantile(null, 0.95),
      perm_p = perm(null, observed$silhouette)
    ),
    means = means, annotation = annotation
  )
})

# People.
person_distance <- \(x) as.dist(1 - cor(x))
shuffle_features <- \(x, seed) {
  set.seed(seed)
  t(apply(x, 1, sample))
}
pam_people <- function(x) {
  d <- person_distance(x)
  map(2:4, \(k) {
    fit <- pam(d, k, diss = TRUE)
    tibble(k = k, silhouette = fit$silinfo$avg.width, groups = list(fit$clustering))
  }) |>
    list_rbind()
}
trait_p <- function(groups) {
  traits <- pairings$training
  map_dbl(set_names(setdiff(names(traits), "subject")), \(tr) {
    v <- set_names(traits[[tr]], traits$subject)[names(groups)]
    ok <- !is.na(v)
    if (n_distinct(groups[ok]) < 2) {
      return(NA_real_)
    }
    kruskal.test(v[ok], groups[ok])$p.value
  })
}
people_inputs <- list(
  modules = inputs_of(modules$eigengenes),
  proteins = map(protein_inputs, \(x) x[order(-apply(x, 1, var))[1:500], , drop = FALSE])
)
people <- imap(people_inputs, \(by_input, features) {
  imap(by_input, \(x, input) {
    null <- map(1:200, \(i) select(pam_people(shuffle_features(x, i)), k, silhouette)) |>
      list_rbind() |>
      summarise(
        null_median = median(silhouette), null_95 = quantile(silhouette, 0.95),
        null_s = list(silhouette), .by = k
      )
    pam_people(x) |>
      left_join(null, by = "k") |>
      mutate(
        perm_p = map2_dbl(null_s, silhouette, perm),
        sizes = map_chr(groups, \(g) paste(table(g), collapse = "/")),
        hr_lr_ari = map_dbl(groups, \(g) mclust::adjustedRandIndex(g, arm_of[names(g)])),
        hr_lr_fisher_p = map_dbl(groups, \(g) fisher.test(table(g, arm_of[names(g)]))$p.value),
        min_trait_p = map_dbl(groups, \(g) min(trait_p(g), na.rm = TRUE)),
        best_trait = map_chr(groups, \(g) names(which.min(trait_p(g)))),
        members = map_chr(groups, \(g) paste(names(g), g, sep = ":", collapse = " ")),
        features = features, input = input
      ) |>
      select(-null_s, -groups)
  }) |>
    list_rbind()
}) |>
  list_rbind()

# Groups: each pathway's and module's proteins, people cut into two clusters.
groups <- c(
  map(gs$sets, \(genes) names(gene_of)[gene_of %in% genes]),
  split(modules$membership$uniprot_id, modules$membership$module)
)
groups <- groups[lengths(groups) >= 10]
split_people <- function(x, members) {
  m <- x[intersect(members, rownames(x)), , drop = FALSE]
  d <- person_distance(m)
  g <- cutree(hclust(d, "average"), k = 2)
  list(groups = g, silhouette = mean_silhouette(g, d))
}
group_scan <- imap(protein_inputs, \(x, input) {
  imap(groups, \(members, group) {
    s <- split_people(x, members)
    tibble(input, group,
      silhouette = s$silhouette,
      ari = mclust::adjustedRandIndex(s$groups, arm_of[names(s$groups)]),
      min_trait_p = suppressWarnings(min(trait_p(s$groups), na.rm = TRUE)),
      split = list(s$groups)
    )
  }) |>
    list_rbind()
}) |>
  list_rbind() |>
  mutate(label = coalesce(unname(feature_label[group]), group))
group_null <- parallel::mclapply(1:100, \(i) {
  set.seed(i)
  shuffled_arm <- set_names(sample(arm_of), names(arm_of))
  imap(protein_inputs, \(x, input) {
    xs <- shuffle_features(x, i)
    real <- filter(group_scan, .data$input == .env$input)
    tibble(input,
      max_silhouette = max(map_dbl(groups, \(m) split_people(xs, m)$silhouette), na.rm = TRUE),
      max_ari = max(map_dbl(real$split, \(g) mclust::adjustedRandIndex(g, shuffled_arm[names(g)])))
    )
  }) |>
    list_rbind()
}, mc.cores = cores) |>
  list_rbind()
group_summary <- group_scan |>
  summarise(
    best_silhouette = max(silhouette, na.rm = TRUE),
    best_silhouette_group = label[which.max(silhouette)], best_ari = max(ari, na.rm = TRUE),
    best_ari_group = label[which.max(ari)], best_trait_p = min(min_trait_p, na.rm = TRUE),
    best_trait_group = label[which.min(min_trait_p)], .by = input
  ) |>
  left_join(
    summarise(group_null, null_s = list(max_silhouette), null_a = list(max_ari), .by = input),
    by = "input"
  ) |>
  mutate(
    perm_p_silhouette = map2_dbl(null_s, best_silhouette, perm),
    perm_p_ari = map2_dbl(null_a, best_ari, perm),
    null_95_silhouette = map_dbl(null_s, quantile, 0.95),
    null_95_ari = map_dbl(null_a, quantile, 0.95),
    best_trait_p_bonferroni = pmin(1, best_trait_p * length(groups) * 8)
  ) |>
  select(-null_s, -null_a)

# Blocks.
column_cut <- \(m) cutree(hclust(as.dist(1 - cor(m)), "average"), k = 4)
block_means <- function(m, rows, cols) {
  expand_grid(row_group = unique(rows[rows != 0]), column_group = unique(cols)) |>
    mutate(mean_abs_t = map2_dbl(row_group, column_group, \(r, c) {
      mean(abs(m[names(rows)[rows == r], cols == c, drop = FALSE]))
    }))
}
column_groups <- column_cut(profiles)
blocks <- block_means(profiles, profile_cut$groups, column_groups)
block_null <- unlist(parallel::mclapply(shuffled_t, \(t_long) {
  m <- profile_matrix(t_long)
  max(block_means(m, cut_profiles(m)$groups, column_cut(m))$mean_abs_t)
}, mc.cores = cores))
blocks <- blocks |>
  mutate(
    perm_p = map_dbl(mean_abs_t, \(v) perm(block_null, v)),
    screens = map_chr(column_group, \(c) {
      paste(colnames(profiles)[column_groups == c], collapse = "; ")
    })
  )
plaid <- function(m) {
  set.seed(1)
  biclust::biclust(m,
    method = biclust::BCPlaid(), verbose = FALSE, back.fit = 2, shuffle = 3,
    fit.model = y ~ m + a + b, iter.startup = 5, iter.layer = 30
  )
}
plaid_fit <- plaid(profiles)
plaid_table <- map(seq_len(plaid_fit@Number), \(b) {
  rows <- rownames(profiles)[plaid_fit@RowxNumber[, b]]
  tibble(
    bicluster = b, n_proteins = length(rows), n_screens = sum(plaid_fit@NumberxCol[b, ]),
    screens = paste(colnames(profiles)[plaid_fit@NumberxCol[b, ]], collapse = "; "),
    main_module = names(which.max(table(module_of[rows]))),
    proteins = paste(head(feature_label[rows], 25), collapse = ", ")
  )
}) |>
  list_rbind()
plaid_null <- map_int(shuffled_t[1:20], \(t_long) plaid(profile_matrix(t_long))@Number)

summary_table <- bind_rows(
  tibble(
    analysis = "profiles", screen = "all pairings", statistic = "silhouette",
    value = profile_cut$silhouette, null_median = median(profile_null, na.rm = TRUE),
    null_95 = quantile(profile_null, 0.95, na.rm = TRUE),
    perm_p = perm(na.omit(profile_null), profile_cut$silhouette)
  ),
  list_rbind(map(trait_groups, "summary")) |>
    transmute(
      analysis = "trait groups", screen = paste(window, outcome), statistic = "silhouette",
      value = silhouette, null_median, null_95, perm_p
    ),
  list_rbind(map(trait_profiles, "summary")) |>
    transmute(
      analysis = "trait profiles", screen = family, statistic = "silhouette",
      value = silhouette, null_median, null_95, perm_p
    ),
  people |>
    slice_max(silhouette, n = 1, by = c(features, input)) |>
    transmute(
      analysis = "people", screen = paste(features, input, "k", k),
      statistic = "silhouette", value = silhouette, null_median, null_95, perm_p
    ),
  group_summary |> transmute(
    analysis = "groups", screen = input,
    statistic = "best silhouette", value = best_silhouette, null_95 = null_95_silhouette,
    perm_p = perm_p_silhouette
  ),
  group_summary |> transmute(
    analysis = "groups", screen = input, statistic = "best HR/LR ARI",
    value = best_ari, null_95 = null_95_ari, perm_p = perm_p_ari
  ),
  tibble(
    analysis = "blocks", screen = "two-way", statistic = "largest block mean |t|",
    value = max(blocks$mean_abs_t), null_median = median(block_null),
    null_95 = quantile(block_null, 0.95), perm_p = perm(block_null, max(blocks$mean_abs_t))
  ),
  tibble(
    analysis = "blocks", screen = "plaid", statistic = "biclusters",
    value = plaid_fit@Number, null_median = median(plaid_null),
    null_95 = quantile(plaid_null, 0.95), perm_p = perm(plaid_null, plaid_fit@Number)
  )
)

write_workbook(
  list(
    summary = summary_table, profile_groups = profile_summary, profile_members = profile_groups,
    profile_means = profile_means,
    trait_groups = list_rbind(map(trait_groups, "groups")),
    trait_profiles = list_rbind(map(trait_profiles, "annotation")),
    trait_profile_means = list_rbind(map(trait_profiles, "means")),
    people = people, groups = select(group_scan, -split), group_summary = group_summary,
    blocks = blocks, plaid = plaid_table
  ),
  c(
    "Per analysis and screen: the statistic, its shuffle null (median, 95th percentile) and p.",
    "Profile groups: size, main co-expression module, top pathway (over-representation).",
    "Profile groups: every grouped protein.",
    "Profile groups: mean t per group and screen.",
    paste(
      "Per trait, groups of co-varying associated proteins: size, direction, module, pathway,",
      "STRING."
    ),
    "Per trait family, protein clusters across pairings: size, module, top pathway.",
    "Per trait family: mean t per cluster and pairing.",
    paste(
      "People: PAM per feature set, input and k; silhouette against feature shuffles; group",
      "sizes; agreement with HR/LR and the traits; members."
    ),
    "Per pathway or module and input: silhouette, ARI with HR/LR, smallest trait p.",
    "Per input: the best group per score against the best over 100 shuffled searches.",
    "Two-way blocks: mean |t| against the largest block of 100 shuffled matrices.",
    "Plaid biclusters of the profile matrix."
  ),
  file.path(out, "05_clustering.xlsx"),
  inputs, paths
)

pdf(file.path(figures, "05_clustering_figures.pdf"), width = 11, height = 8.5, useDingbats = TRUE)
print(plot_chance(summary_table, "value",
  "Every clustering against its shuffle null",
  paste(
    "Silhouette, agreement or block statistic per analysis; grey bar to the shuffles' 95th",
    "percentile, tick at their median. Table: c_data/05_clustering.xlsx, summary."
  ),
  facet = ~analysis
))
print(
  ggplot(profile_means, aes(screen, factor(group), fill = mean_t)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", high = "#B2182B", name = "mean t") +
    labs(
      x = NULL, y = "profile group", title = "Protein groups by association profile",
      caption = "Table: profile_groups, profile_means."
    ) +
    page_theme +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
)
print(
  list_rbind(map(trait_profiles, "means")) |>
    ggplot(aes(window, mean_t, group = group, colour = factor(group))) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_line() +
    geom_point(aes(size = n_proteins)) +
    facet_wrap(~family) +
    labs(
      x = NULL, y = "mean t", colour = "group",
      title = "Per trait: protein clusters by association across pairings",
      caption = "Table: trait_profiles, trait_profile_means."
    ) +
    page_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
)
print(
  group_scan |>
    ggplot(aes(silhouette, ari)) +
    geom_point(alpha = 0.25, size = 0.6) +
    geom_vline(aes(xintercept = null_95_silhouette),
      data = group_summary, linetype = 2,
      colour = "firebrick"
    ) +
    geom_hline(aes(yintercept = null_95_ari),
      data = group_summary, linetype = 2,
      colour = "firebrick"
    ) +
    facet_wrap(~input) +
    labs(
      x = "silhouette of the two-way split of people", y = "adjusted Rand index with HR/LR",
      title = "People split within every pathway and module",
      caption = stringr::str_wrap(paste(
        "Dashed lines: 95th percentile of the best score over 100",
        "shuffled searches; a point beyond a line beats the search's own chance. Table: groups."
      ), 160)
    ) +
    page_theme
)
invisible(dev.off())
