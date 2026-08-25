# Is there any group structure to find, anywhere?
#
# The project asked this repeatedly before settling on a continuous design, and
# the answer decides whether that design is a preference or a necessity. It
# needs asking before any partition is drawn, not after: at 15 subjects with
# more features than observations, every clustering method returns groups, so a
# cluster count on its own carries no information.
#
# Five methods, chosen because they disagree about what a cluster is and
# because two of them can return "no clusters" as an answer:
#
#   dip test        unimodality of the leading projection. Adolfsson, Ackerman
#                   & Brownstein (2019) find this the most reliable
#                   clusterability check; the null is one group, which is the
#                   question. Hartigan & Hartigan (1985).
#   gap statistic   compares within-cluster dispersion against a uniform
#                   reference and can select k = 1. Tibshirani, Walther &
#                   Hastie (2001).
#   silhouette      magnitude, not just the arg-max. Below 0.5 is "weak, could
#                   be artificial" on the Kaufman-Rousseeuw scale.
#   mclust BIC      with a simulated single-Gaussian null, because BIC alone is
#                   not interpretable here. See the caveat below.
#   consensus PAC   subsampling stability against shuffled data. Senbabaoglu et
#                   al. (2014) showed consensus clustering splits unimodal null
#                   data into apparently stable clusters, so the shuffled
#                   comparison is the whole point.
#
# The mclust column is reported and not trusted. At p = 1902 and K = 2 even its
# most constrained family estimates thousands of covariance parameters from 16
# observations, and mclust drops unfittable models from the BIC table silently
# rather than warning, so the curve spans only the degenerate-but-estimable
# subset. Bouveyron & Brunet-Saumard (2014) is the review; this is why the
# column exists to be contradicted by the other four.

pacman::p_load(
  here, dplyr, purrr, tibble, readr, cluster, diptest, MASS, openxlsx
)

source(here("functions", "association.R"))

OUT_DIR <- here("03_Features", "c_data")
N_NULL <- 199
GAP_B <- 100
PAC_B <- 200
K_MAX <- 5

pheno <- phenotype_table()
features <- feature_matrices()
proteins <- features$proteins[stats::complete.cases(features$proteins), ]

phenotype_only <- pheno |>
  dplyr::select(-"subject", -"group_arm") |>
  as.matrix() |>
  (\(m) m[stats::complete.cases(m), , drop = FALSE])()

spaces <- list(
  `subjects, baseline modules` = t(subject_window(features$modules, "T1")),
  `subjects, T2 modules` = t(subject_window(features$modules, "T2")),
  `subjects, baseline proteins` = t(subject_window(proteins, "T1")),
  `all 45 samples, proteins` = t(proteins),
  `subjects, ten phenotypes` = phenotype_only
)

# Consensus stability without the package: subsample, cluster, and record how
# often each pair lands together. PAC is the share of pair-consensus values
# stuck in the ambiguous middle, so low means a clean split.
pac <- function(x, k, b = PAC_B, frac = 0.8, lo = 0.1, hi = 0.9) {
  n <- nrow(x)
  co <- matrix(0, n, n)
  seen <- matrix(0, n, n)
  for (i in seq_len(b)) {
    idx <- sample(n, round(frac * n))
    cl <- cluster::pam(stats::dist(x[idx, , drop = FALSE]), k)$clustering
    seen[idx, idx] <- seen[idx, idx] + 1
    together <- outer(cl, cl, "==")
    co[idx, idx] <- co[idx, idx] + together
  }
  m <- co / pmax(seen, 1)
  v <- m[upper.tri(m)]
  mean(v > lo & v < hi)
}

best_g <- function(z) {
  withr::local_package("mclust")
  f <- try(mclust::Mclust(z, G = 1:K_MAX, verbose = FALSE), silent = TRUE)
  if (inherits(f, "try-error") || is.null(f)) NA_integer_ else f$G
}

assess <- function(name, x) {
  x <- scale(x)
  pc <- stats::prcomp(x)$x
  d <- stats::dist(x)

  dip <- diptest::dip.test(pc[, 1])

  gap <- cluster::clusGap(x,
    FUN = stats::kmeans, K.max = K_MAX, B = GAP_B,
    nstart = 20
  )
  gap_k <- cluster::maxSE(
    gap$Tab[, "gap"], gap$Tab[, "SE.sim"],
    method = "Tibs2001SEmax"
  )

  sil <- vapply(2:K_MAX, function(g) {
    mean(cluster::silhouette(cluster::pam(d, g))[, 3])
  }, numeric(1))

  scores <- pc[, 1:2, drop = FALSE]
  observed_g <- best_g(scores)
  null_g <- vapply(seq_len(N_NULL), function(i) {
    best_g(MASS::mvrnorm(nrow(scores), colMeans(scores), stats::cov(scores)))
  }, numeric(1))

  tibble(
    space = name, n = nrow(x), p = ncol(x),
    dip_D = unname(dip$statistic), dip_p = dip$p.value,
    gap_k = gap_k,
    silhouette_best_k = (2:K_MAX)[which.max(sil)],
    silhouette_max = max(sil),
    mclust_g = observed_g,
    mclust_null_g_above_1 = mean(null_g > 1, na.rm = TRUE),
    pac_k2 = pac(x, 2),
    pac_k2_shuffled = mean(replicate(10, pac(apply(x, 2, sample), 2)))
  )
}

set.seed(42)
clusterability <- imap_dfr(spaces, function(x, nm) assess(nm, x))

# The composite is the one variable the original HR/LR label was cut from, and
# it is one-dimensional, so it needs no projection and gets its own row.
composite_dip <- diptest::dip.test(pheno$comp_hypertrophy)
composite <- tibble(
  variable = "comp_hypertrophy",
  n = length(pheno$comp_hypertrophy),
  dip_D = unname(composite_dip$statistic),
  dip_p = composite_dip$p.value
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write.xlsx(
  list(clusterability = clusterability, composite_modality = composite),
  file.path(OUT_DIR, "06_clusterability.xlsx")
)

print(as.data.frame(clusterability |> dplyr::select(
  space, n, p, dip_p, gap_k, silhouette_max, mclust_g, mclust_null_g_above_1
)), digits = 3)
print(as.data.frame(composite), digits = 3)
message(
  "\n", sum(clusterability$dip_p < 0.05), " of ", nrow(clusterability),
  " feature spaces reject unimodality; ",
  sum(clusterability$gap_k > 1), " have a gap-statistic k above 1"
)
