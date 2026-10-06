# People clustered on their proteome, no phenotype or group label used, separately at each input:
# T1, T2 and T3 levels, training change (T2 - T1) and acute change (T3 - T2). Features are the
# trunk's module eigengenes, z-scored per input. The K rule matches the phenotype clusters: PAM,
# the gap statistic for K = 1 to 3 with the one-standard-error rule (Tibshirani, Walther & Hastie
# 2001), mean silhouette for K = 2 and 3, group sizes and bootstrap Jaccard stability per cluster
# (fpc::clusterboot). imposed_k fixes K for an input by hand, with its reason. Agreement with the
# HR/LR split is reported as a description (adjusted Rand index), never used.

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

inputs_used <- c("T1", "T2", "T3", "training", "acute")
imposed_k <- set_names(rep(2L, 5), inputs_used)
imposed_reason <- set_names(rep(paste(
  "Gap statistic selects K = 1 at every input, which is reported as the result; K = 2 imposed",
  "by decision (2026-10-06) so the exploratory arc can compare the two halves the proteome",
  "suggests. Read each split with its bootstrap Jaccard."
), 5), inputs_used)
k_max <- 3
n_boot <- 500

stage <- here("03_Cluster_groups", "Proteome", "01_Cluster")
out <- file.path(stage, "c_data")
figures <- file.path(stage, "b_reports")
walk(c(out, figures), dir.create, recursive = TRUE, showWarnings = FALSE)
inputs <- c(
  modules = "01_Trunk/06_Modules/c_data/modules.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
stopifnot(file.exists(paths))

modules <- readRDS(paths[["modules"]])
meta <- as_tibble(modules$meta)
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
arm_of <- set_names(phenotype$arm, phenotype$subject)
values <- modules$eigengenes
by_input <- list(
  T1 = subject_values(values, meta, "T1"), T2 = subject_values(values, meta, "T2"),
  T3 = subject_values(values, meta, "T3"),
  training = subject_change(values, meta, "T1", "T2"),
  acute = subject_change(values, meta, "T2", "T3")
)

pam1 <- \(x, k) list(cluster = pam(x, k, cluster.only = TRUE))
per_input <- imap(by_input, \(m, input) {
  x <- scale(t(m))
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
          sizes = paste(table(fit$clustering), collapse = "/"),
          hr_lr_ari = mclust::adjustedRandIndex(fit$clustering, arm_of[rownames(x)])
        )
      }) |>
        list_rbind(),
      by = "k"
    ) |>
    left_join(
      summarise(stability,
        min_jaccard = min(jaccard), jaccard = paste(round(jaccard, 2), collapse = "/"), .by = k
      ),
      by = "k"
    ) |>
    mutate(input = input, .before = 1)
  k_used <- coalesce(as.integer(imposed_k[[input]]), k_gap)
  list(
    diagnostics = diagnostics,
    stability = mutate(stability, input = input, .before = 1),
    choice = tibble(
      input,
      n_subjects = nrow(x), n_features = ncol(x), k_gap, k_imposed = imposed_k[[input]],
      imposed_reason = imposed_reason[[input]], k_used, grouped = k_used > 1
    ),
    assignments = tibble(
      input,
      subject = rownames(x), arm = unname(arm_of[rownames(x)]),
      k2 = solutions[["2"]]$clustering, k3 = solutions[["3"]]$clustering,
      cluster = if (k_used > 1) solutions[[as.character(k_used)]]$clustering else 1L
    ),
    gap = mutate(gap_table, input = input, k_gap = k_gap)
  )
})
results <- map(set_names(names(per_input[[1]])), \(sheet) list_rbind(map(per_input, sheet)))

pdf(file.path(figures, "01_cluster_figures.pdf"), width = 11, height = 8.5)
print(
  ggplot(results$gap, aes(k, gap)) +
    geom_pointrange(aes(ymin = gap - SE.sim, ymax = gap + SE.sim)) +
    geom_line() +
    geom_vline(aes(xintercept = k_gap), linetype = "dashed") +
    facet_wrap(~input, nrow = 1) +
    scale_x_continuous(breaks = seq_len(k_max)) +
    labs(
      title = "Proteome clusters: gap statistic per input, PAM on z-scored module eigengenes",
      subtitle = sprintf(
        "Dashed: K selected by the one-SE rule (Tibs2001SEmax); %d bootstrap reference sets",
        n_boot
      ),
      x = "K", y = "Gap (± SE)"
    ) +
    theme_bw(base_size = 10)
)
invisible(dev.off())

write_workbook(
  results[c("choice", "diagnostics", "stability", "assignments")],
  c(
    "Per input: subjects, features, K selected by the gap statistic, any imposed K and reason.",
    "Per input and K: gap, SE, mean silhouette, sizes, HR/LR adjusted Rand index, Jaccard.",
    "Per input, K and cluster: mean bootstrap Jaccard (fpc::clusterboot, PAM).",
    "Per input and subject: HR/LR (for description), the K = 2 and K = 3 clusters, cluster used."
  ),
  file.path(out, "01_cluster.xlsx"), inputs, paths
)
