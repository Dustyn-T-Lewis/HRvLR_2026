# People clustered on their muscle-size response, before any proteome result is opened. PAM on the
# z-scored % changes over training. K diagnostics: the gap statistic for K = 1 to 3 with the
# one-standard-error rule of Tibshirani, Walther & Hastie (2001), mean silhouette for K = 2 and 3,
# group sizes, and bootstrap Jaccard stability per cluster (fpc::clusterboot; Hennig 2007).
# imposed_k fixes K by hand; both the data-selected and the imposed K are recorded with the reason.
# If K = 1 is selected and nothing is imposed, the family is not grouped and this arc stops here.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
  library(cluster)
})
source(here("R", "helpers.R"))

family <- "Muscle_size"
traits <- c("pct_fcsa_mixed", "pct_fcsa_type1", "pct_fcsa_type2", "pct_mcsa")
imposed_k <- 2L
imposed_reason <- paste0(
  "Gap statistic selects K = 1; K = 2 imposed by decision (2026-10-06). It reproduc",
  "es the HR/LR split exactly, so this arc repeats that comparison under a size-clu",
  "ster label."
)
k_max <- 3
n_boot <- 500

stage <- here("03_Cluster_groups", family, "01_Clustering")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(phenotype = "00_Input/phenotype.csv")
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
change <- phenotype_tables(phenotype)$pairings$training
complete <- change |>
  select(subject, all_of(traits)) |>
  drop_na()
x <- scale(as.matrix(select(complete, -subject)))
rownames(x) <- complete$subject

pam1 <- \(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
set.seed(1)
gap <- clusGap(x, FUNcluster = pam1, K.max = k_max, B = n_boot)
gap_table <- as_tibble(gap$Tab) |> mutate(k = seq_len(n()), .before = 1)
k_gap <- maxSE(gap_table$gap, gap_table$SE.sim, method = "Tibs2001SEmax")

solutions <- map(set_names(2:k_max), \(k) pam(x, k))
stability <- map(set_names(2:k_max), \(k) {
  boot <- fpc::clusterboot(x,
    B = n_boot, clustermethod = fpc::claraCBI, k = k, usepam = TRUE,
    seed = 1, count = FALSE
  )
  tibble(k = k, cluster = seq_along(boot$bootmean), jaccard = boot$bootmean)
}) |>
  list_rbind()
diagnostics <- gap_table |>
  left_join(
    imap(solutions, \(fit, k) {
      tibble(
        k = as.integer(k), mean_silhouette = fit$silinfo$avg.width,
        sizes = paste(table(fit$clustering), collapse = "/")
      )
    }) |>
      list_rbind(),
    by = "k"
  ) |>
  left_join(
    summarise(stability,
      min_jaccard = min(jaccard),
      jaccard = paste(round(jaccard, 2), collapse = "/"), .by = k
    ),
    by = "k"
  )

k_used <- coalesce(imposed_k, k_gap)
choice <- tibble(
  family,
  n_subjects = nrow(x), traits = paste(traits, collapse = ", "),
  k_gap = k_gap, k_imposed = imposed_k, imposed_reason, k_used,
  grouped = k_used > 1
)
assignments <- complete |>
  mutate(
    k2 = solutions[["2"]]$clustering, k3 = solutions[["3"]]$clustering,
    cluster = if (k_used > 1) solutions[[as.character(k_used)]]$clustering else 1L
  ) |>
  right_join(select(change, subject), by = "subject") |>
  arrange(cluster, subject)

# The trunk's inputs, relabelled for this arc: arm and group carry the clusters, named Hi and Lo by
# the mean of the family's first trait, and the samples of any unassigned subject drop out. The
# arc's stages read these copies, so they run the same code as 02_HRvLR_groups on new labels.
if (k_used == 2) {
  label_of <- assignments |>
    filter(!is.na(cluster)) |>
    summarise(lead = mean(.data[[traits[1]]]), .by = cluster) |>
    mutate(label = if_else(lead == max(lead), "Hi", "Lo")) |>
    select(cluster, label)
  group_of <- assignments |>
    inner_join(label_of, by = "cluster") |>
    with(set_names(label, subject))
  relabel <- function(meta) {
    meta <- meta[meta$subject %in% names(group_of), , drop = FALSE]
    meta$arm <- unname(group_of[meta$subject])
    meta$group <- paste(meta$arm, meta$timepoint, sep = "_")
    meta
  }
  trunk <- \(path) readRDS(here("01_Trunk", path))
  inputs_dir <- file.path(out, "inputs")
  dir.create(inputs_dir, showWarnings = FALSE)
  proteins <- trunk("03_Quantification/c_data/proteins.rds")
  proteins$metadata <- relabel(proteins$metadata)
  proteins$proteins <- proteins$proteins[, proteins$metadata$sample_id]
  saveRDS(proteins, file.path(inputs_dir, "proteins.rds"), compress = "xz")
  walk(c("DAList_corrected", "DAList_limpa"), \(name) {
    d <- trunk(file.path("03_Quantification/c_data", paste0(name, ".rds")))
    d$metadata <- relabel(d$metadata)
    d$data <- d$data[, d$metadata$sample_id]
    d$observed <- d$observed[, d$metadata$sample_id]
    saveRDS(d, file.path(inputs_dir, paste0(name, ".rds")), compress = "xz")
  })
  scores <- trunk("05_Set_Scores/c_data/singscore.rds")
  scores$scores <- scores$scores[, proteins$metadata$sample_id]
  saveRDS(scores, file.path(inputs_dir, "singscore.rds"), compress = "xz")
  modules <- trunk("06_Module_Building/c_data/modules.rds")
  modules$meta <- relabel(modules$meta)
  modules$eigengenes <- modules$eigengenes[, modules$meta$sample_id]
  saveRDS(modules, file.path(inputs_dir, "modules.rds"), compress = "xz")
  saveRDS(
    list(group_of = group_of, labels = label_of, family = family, traits = traits),
    file.path(inputs_dir, "groups.rds")
  )
}

pdf(file.path(figures, "clustering_report.pdf"), width = 11, height = 8.5)
print(
  ggplot(gap_table, aes(k, gap)) +
    geom_pointrange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
    geom_line() +
    geom_vline(xintercept = k_gap, linetype = "dashed") +
    scale_x_continuous(breaks = seq_len(k_max)) +
    labs(
      title = sprintf("%s: gap statistic, PAM on %d subjects", family, nrow(x)),
      subtitle = sprintf(
        "One-SE rule (Tibs2001SEmax) selects K = %d; K used = %d. %d bootstrap reference sets",
        k_gap, k_used, n_boot
      ),
      x = "K", y = "Gap (± SE)"
    ) +
    theme_bw(base_size = 10)
)
pcs <- prcomp(x)
print(
  as_tibble(pcs$x[, 1:2], rownames = "subject") |>
    left_join(select(assignments, subject, k2, k3), by = "subject") |>
    pivot_longer(c(k2, k3), names_to = "solution", values_to = "cluster") |>
    ggplot(aes(PC1, PC2, colour = factor(cluster), label = subject)) +
    geom_point(size = 2) +
    geom_text(size = 2.2, vjust = -0.8, show.legend = FALSE) +
    facet_wrap(~solution) +
    labs(
      title = sprintf("%s: subjects on the first two components of the scaled traits", family),
      subtitle = paste(traits, collapse = ", "), colour = "Cluster"
    ) +
    theme_bw(base_size = 10)
)
invisible(dev.off())

# Each person on the first two components of the scaled traits: the space the clusters were found
# in, for the arc's Figure 1.
space <- as_tibble(pcs$x[, 1:2], rownames = "subject") |>
  left_join(select(assignments, subject, cluster), by = "subject") |>
  mutate(
    arm = set_names(phenotype$arm, phenotype$subject)[subject],
    group = if (k_used == 2) unname(group_of[subject]) else NA_character_,
    .after = subject
  )
space_variance <- tibble(
  component = c("PC1", "PC2"), variance_explained = summary(pcs)$importance[2, 1:2]
)

write_workbook(
  list(
    choice = choice, diagnostics = diagnostics, stability = stability,
    assignments = assignments, space = space, space_variance = space_variance
  ),
  c(
    "K selected by the gap statistic, any imposed K and its reason, and the K used downstream.",
    "Per K: gap, SE, mean silhouette, group sizes, bootstrap Jaccard per cluster.",
    "Per K and cluster: mean bootstrap Jaccard (fpc::clusterboot, PAM).",
    "Per subject: the traits, the K = 2 and K = 3 clusters, and the cluster used downstream.",
    "Per subject: HR/LR arm, group, cluster, and PC1 and PC2 of the scaled traits.",
    "Share of variance on PC1 and PC2 of the scaled traits."
  ),
  file.path(out, "clustering.xlsx"), inputs, paths
)
