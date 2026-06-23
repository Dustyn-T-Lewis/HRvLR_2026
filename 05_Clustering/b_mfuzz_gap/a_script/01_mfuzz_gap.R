# Fuzzy c-means clustering of the HR-LR gap trajectory across T1/T2/T3.
# Helpers below follow Mfuzz / Schwämmle & Jensen (2010) Bioinformatics 26:2923.
# e1071::cmeans replaces Mfuzz's Fclust (which requires tcltk/XQuartz).

pacman::p_load(here, e1071, ggplot2)
source(here("05_Clustering", "_shared_inputs.R"))
i <- load_clustering_inputs()

# Z-score each protein across timepoints; drop zero-variance rows
gap_z <- t(scale(t(i$gap)))
bad <- apply(gap_z, 1, function(x) any(is.nan(x)) | any(is.na(x)))
gap_z <- gap_z[!bad, ]
message("Zero-variance proteins dropped: ", sum(bad))

# Fuzzifier m from Schwämmle & Jensen 2010 eq. 1; D = 3 (timepoints)
n <- nrow(gap_z)
D <- 3L
m <- 1 +
  (1418 / n + 22.05) * D^(-2) +
  (12.33 / n + 0.243) * D^(-0.0406 * log(n) - 0.1134)
message("Fuzzifier m (mestimate): ", round(m, 4))

# Min inter-centroid distance (Dmin) for c = 2:8; cap at 4–6 for 3 timepoints
set.seed(42)
dmin <- vapply(2:8, function(c) {
  fit <- cmeans(gap_z, centers = c, m = m, iter.max = 200)
  d <- as.matrix(dist(fit$centers))
  diag(d) <- Inf
  min(d)
}, numeric(1))
names(dmin) <- as.character(2:8)

# c chosen from the printed Dmin curve (manual, Mfuzz-style),
# capped to the 4–6 range justified by 3 timepoints
chosen_c <- 4L
message(
  "Dmin values c=2:8: ",
  paste(round(dmin, 3), collapse = " "),
  " | chosen c=", chosen_c
)

set.seed(42)
fit <- cmeans(gap_z, centers = chosen_c, m = m, iter.max = 200)

max_mem <- apply(fit$membership, 1, max)
hard_cluster <- apply(fit$membership, 1, which.max)
n_below <- sum(max_mem < 0.5)
message("Proteins below 0.5 core cutoff: ", n_below, " / ", nrow(gap_z))

# Seed stability: ARI (label-invariant, no permutation matching needed)
# Uses mclust::adjustedRandIndex when available; falls back to pair-counting ARI
ari_fn <- if (requireNamespace("mclust", quietly = TRUE)) {
  mclust::adjustedRandIndex
} else {
  function(a, b) {
    tab <- table(a, b)
    n_pts <- sum(tab)
    sum_a <- choose(rowSums(tab), 2)
    sum_b <- choose(colSums(tab), 2)
    sum_ab <- sum(choose(tab, 2))
    expected <- sum(sum_a) * sum(sum_b) / choose(n_pts, 2)
    max_index <- (sum(sum_a) + sum(sum_b)) / 2
    (sum_ab - expected) / (max_index - expected)
  }
}

all_hard <- lapply(seq_len(25), function(s) {
  set.seed(s)
  f <- cmeans(gap_z, centers = chosen_c, m = m, iter.max = 200)
  apply(f$membership, 1, which.max)
})
pairs <- combn(25L, 2L)
ag <- apply(pairs, 2, function(p) ari_fn(all_hard[[p[1]]], all_hard[[p[2]]]))
message(
  "Seed stability (25 seeds, ", ncol(pairs), " pairs): ",
  "median ARI = ", round(median(ag), 4),
  ", min ARI = ", round(min(ag), 4)
)

# ME = PC1 of per-sample abundances for each cluster's proteins;
# sign-fixed so positive ME = higher mean abundance.
# scale. = FALSE because proteins are already on the log2 scale.
compute_me <- function(proteins, abund) {
  mat <- abund[proteins, , drop = FALSE]
  pc <- prcomp(t(mat), center = TRUE, scale. = FALSE)
  me <- pc$x[, 1]
  if (cor(me, colMeans(mat)) < 0) me <- -me
  me
}

clust_ids <- sort(unique(hard_cluster))
eigen_rows <- lapply(clust_ids, function(cl) {
  prots <- names(hard_cluster)[hard_cluster == cl]
  me <- compute_me(prots, i$abund)
  # Guard: PC1 inherits sample names from colnames(abund); must match meta
  stopifnot(identical(names(me), i$meta$sample_id))
  data.frame(
    sample_id = names(me),
    subject = i$meta$subject,
    group_arm = i$meta$group_arm,
    timepoint = i$meta$timepoint,
    group_id = cl,
    ME = me
  )
})
eigengene <- do.call(rbind, eigen_rows)

membership <- data.frame(
  protein_id = names(hard_cluster),
  group_id = hard_cluster,
  membership_weight = max_mem
)

write.csv(
  membership,
  here("05_Clustering", "b_mfuzz_gap", "c_data", "membership.csv"),
  row.names = FALSE
)
write.csv(
  eigengene,
  here("05_Clustering", "b_mfuzz_gap", "c_data", "eigengene.csv"),
  row.names = FALSE
)

# Gap trajectory plots (membership-shaded, Mfuzz-style)
tp_labels <- c("T1" = "Baseline", "T2" = "Trained", "T3" = "Acute")

gap_df <- as.data.frame(gap_z)
gap_df$protein_id <- rownames(gap_z)
gap_df$cluster <- factor(hard_cluster)
gap_df$mem <- max_mem

gap_long <- reshape(
  gap_df,
  direction = "long",
  varying = c("T1", "T2", "T3"),
  v.names = "z",
  timevar = "tp",
  times = c("T1", "T2", "T3"),
  idvar = "protein_id"
)
gap_long$tp_label <- factor(
  tp_labels[gap_long$tp],
  levels = c("Baseline", "Trained", "Acute")
)

clust_mean <- aggregate(z ~ cluster + tp_label, data = gap_long, FUN = mean)

pdf(
  here("05_Clustering", "b_mfuzz_gap", "b_reports", "gap_clusters.pdf"),
  width = 9, height = 7
)
p <- ggplot(gap_long, aes(x = tp_label, y = z, group = protein_id)) +
  geom_line(aes(colour = mem), alpha = 0.4, linewidth = 0.3) +
  geom_line(
    data = clust_mean, aes(x = tp_label, y = z, group = cluster),
    colour = "black", linewidth = 1.2, inherit.aes = FALSE
  ) +
  scale_colour_gradient(
    low = "grey80", high = "steelblue",
    name = "Membership"
  ) +
  facet_wrap(~cluster, labeller = label_both) +
  labs(
    x = NULL, y = "z-scored HR-LR gap",
    title = "Fuzzy c-means gap-trajectory clusters (z-scored HR-LR gap)"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
print(p)
dev.off()

# Group-mean abundance profiles per cluster
abund_long <- as.data.frame(t(i$abund))
abund_long$sample_id <- rownames(abund_long)
abund_long <- merge(abund_long, i$meta, by = "sample_id")

shape_rows <- lapply(clust_ids, function(cl) {
  prots <- names(hard_cluster)[hard_cluster == cl]
  sub <- abund_long[, c("sample_id", "group_arm", "timepoint", prots)]
  sub$mean_abund <- rowMeans(sub[, prots, drop = FALSE])
  agg <- aggregate(
    mean_abund ~ group_arm + timepoint,
    data = sub,
    FUN = mean
  )
  agg$cluster <- cl
  agg
})
shape_df <- do.call(rbind, shape_rows)
shape_df$timepoint <- factor(
  shape_df$timepoint,
  levels = c("Baseline", "Trained", "Acute")
)

pdf(
  here("05_Clustering", "b_mfuzz_gap", "b_reports", "shape_profiles.pdf"),
  width = 9, height = 7
)
p2 <- ggplot(
  shape_df,
  aes(x = timepoint, y = mean_abund, colour = group_arm, group = group_arm)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  facet_wrap(~cluster, labeller = label_both, scales = "free_y") +
  labs(
    x = NULL, y = "Mean log2-abundance",
    colour = "Group",
    title = "Cluster shape profiles: HR vs LR"
  ) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank())
print(p2)
dev.off()

message(
  "Done. membership: ", nrow(membership), " rows | ",
  "eigengene: ", nrow(eigengene), " rows (",
  chosen_c, " clusters × 45 samples)"
)
