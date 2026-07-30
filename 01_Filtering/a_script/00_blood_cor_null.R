# Why blood_cor stopped gating removal.
#
# The cut was 0.45, justified as a rounding above a permutation null recorded as
# 0.43. This recomputes that null on the current input and reports it by
# observation count, which is where it fails: a protein seen in 15 samples and
# one seen in 48 do not draw |rho| from the same distribution, so no
# single cut is calibrated for both. Removal moved to the curated blood list in
# filter_config.R. Kept as the evidence for that move, and rerun if 00_input
# changes.

pacman::p_load(here, readxl, readr, dplyr, tibble)
source(here("functions", "blood_cor_null.R"))
source(here("01_Filtering", "a_script", "filter_config.R"))
cfg <- filter_cfg

raw <- read_excel(here("00_input", "HRvLR_raw.xlsx"))
metadata <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE)

annotation <- raw[, c("uniprot_id", "gene")]
intensity <- data.matrix(raw[, metadata$Col_ID])

log_int <- log2(intensity)
log_int[!is.finite(log_int)] <- NA
blood_index <- colMeans(
  log_int[annotation$gene %in% cfg$blood_anchor, , drop = FALSE],
  na.rm = TRUE
)

# 01_run_filtering.R sets blood_cor to NA below this many observations, because
# a protein seen in three samples reaches |rho| = 1 on any permutation. The null
# has to be drawn from the same proteins that can actually receive a call, or it
# saturates at 1 and reports that nothing is ever significant.
testable <- rowSums(!is.na(log_int)) >= cfg$miss_min_reps * 3
log_int <- log_int[testable, ]

N_PERM <- 300
quantiles <- c(0.99, 0.999, 0.9999)
null_cuts <- vapply(quantiles, function(p) {
  blood_cor_null_quantile(log_int, blood_index, n_perm = N_PERM, prob = p)
}, numeric(1))

cat(sprintf(
  "%d testable proteins of %d, x %d samples, %d permutations\n\n",
  nrow(log_int), length(testable), ncol(log_int), N_PERM
))
print(as.data.frame(tibble(
  quantile = quantiles, null_abs_rho = round(null_cuts, 4)
)))
# The pooled quantile hides the problem, so report it by observation count too.
strata <- list(c(15, 25), c(26, 40), c(41, 48))
n_obs <- rowSums(!is.na(log_int))
by_n <- vapply(strata, function(r) {
  sel <- n_obs >= r[1] & n_obs <= r[2]
  blood_cor_null_quantile(log_int[sel, ], blood_index, n_perm = N_PERM)
}, numeric(1))

cat("\n99.9% of the null by observation count:\n")
print(as.data.frame(tibble(
  observations = vapply(strata, function(r) sprintf("%d-%d", r[1], r[2]), ""),
  proteins = vapply(strata, function(r) sum(n_obs >= r[1] & n_obs <= r[2]), 0L),
  null_abs_rho = round(by_n, 3)
)))
cat("\nOne flat cut cannot serve that range; removal is by curated list.\n")
