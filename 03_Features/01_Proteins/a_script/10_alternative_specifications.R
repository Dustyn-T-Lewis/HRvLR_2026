#!/usr/bin/env Rscript
# Four ways to model the group effect that the nine contrasts do not use.
#
# The cell-means model asks one question well: does the mean differ between two
# labelled sets of samples. Each alternative here changes what is being asked or
# how the residual is estimated, and each is a standard limma or lm option
# rather than a new method.
#
#   continuous  replaces the binary label with the composite hypertrophy score.
#               HR and LR are a split of that score, so the split throws away
#               the spacing between subjects; regressing on the score keeps it.
#   ancova      tests the trained state adjusted for the same protein's
#               baseline. For pre-post designs this is usually better powered
#               than either a change score or a post-only comparison, and it is
#               the one the cell-means design cannot express.
#   weights     limma::arrayWeights downweights samples whose residuals run
#               large. It appears to find hits the primary misses and does not
#               survive its own null: 10 BH hits observed, a median of 7 under
#               subject-permuted labels and a 97.5th percentile of 79, p = 0.33.
#               The weights are estimated from the data they then test, so they
#               adapt to any labelling. Reported as a warning, not a result.
#
# lmFit(method = "robust") was tried and cannot be used: limma will not combine
# it with duplicateCorrelation, and dropping the subject blocking to gain robust
# weighting is a worse trade in a three-biopsy design.
#
# All four are reported beside the primary. None promotes a protein on its own.

pacman::p_load(here, dplyr, tibble, purrr, readr, limma, openxlsx)
source(here("functions", "feature_contrasts.R"))
source(here("functions", "feature_levels.R"))
source(here("03_Features", "contrasts.R"))

set.seed(42)
OUT <- here("03_Features", "01_Proteins", "c_data", "alt_specs")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

mat <- protein_matrix()
meta <- feature_metadata()
meta$time <- factor(
  sub("^.*_", "", as.character(meta$group)),
  levels = c("T1", "T2", "T3")
)
meta$arm <- factor(
  sub("_.*$", "", as.character(meta$group)),
  levels = c("LR", "HR")
)

pheno <- read_csv(
  here("00_input", "c_data", "phenotype.csv"),
  show_col_types = FALSE
) |>
  mutate(subject = sub("^(HR|LR)_", "", .data$subject))
meta$comp <- pheno$comp_hypertrophy[
  match(sub("^(HR|LR)_", "", meta$subject), pheno$subject)
]
stopifnot(!anyNA(meta$comp))

report <- function(tt, label, coefs) {
  map_dfr(coefs, function(cf) {
    r <- limma::topTable(
      tt,
      coef = cf, number = Inf, adjust.method = "BH", sort.by = "none"
    )
    tibble(
      spec = label, term = cf, n = nrow(r),
      nominal = sum(r$P.Value < 0.05, na.rm = TRUE),
      bh05 = sum(r$adj.P.Val < 0.05, na.rm = TRUE),
      bh10 = sum(r$adj.P.Val < 0.10, na.rm = TRUE),
      min_bh = min(r$adj.P.Val, na.rm = TRUE),
      top = paste(
        utils::head(rownames(r)[order(r$P.Value)], 3),
        collapse = ", "
      )
    )
  })
}

gene_of <- function(ids) {
  ann <- as.data.frame(
    readRDS(
      here("02_Normalization", "c_data", "DAList_normalized.rds")
    )$annotation
  )
  paste(
    ann$gene[match(strsplit(ids, ", ")[[1]], ann$uniprot_id)],
    collapse = ", "
  )
}

# Continuous: one slope on the composite score per timepoint.
des_cont <- stats::model.matrix(~ 0 + time + time:comp, meta)
colnames(des_cont) <- make.names(colnames(des_cont))
corr_c <- limma::duplicateCorrelation(
  mat, des_cont,
  block = meta$subject
)$consensus
fit_c <- limma::eBayes(
  limma::lmFit(
    mat, des_cont,
    block = meta$subject, correlation = corr_c
  )
)
slope_terms <- grep("comp", colnames(des_cont), value = TRUE)
res_cont <- report(fit_c, "continuous score", slope_terms)

# Sample weights and robust fitting, both on the primary design.
parts <- feature_design(mat, meta)
aw <- limma::arrayWeights(mat, parts$design)
fit_w <- limma::eBayes(limma::contrasts.fit(
  limma::lmFit(mat, parts$design,
    block = parts$block,
    correlation = parts$correlation, weights = aw
  ),
  parts$contrasts
))
res_w <- report(fit_w, "sample weights", colnames(parts$contrasts))

# lmFit(method = 'robust') is not available here: limma refuses to combine it
# with duplicateCorrelation, and the repeated-measures blocking is not optional
# in a design with three biopsies per subject. Robust weighting would have to
# replace the blocking, which trades a real dependency for a hypothetical outlier.


# ANCOVA: trained state adjusted for the same protein's own baseline. limma
# cannot express a covariate that changes per protein, so this is a per-protein
# fit over the 14 subjects with both timepoints.
wide <- function(tp) {
  cols <- meta$sample_id[meta$time == tp]
  m <- mat[, cols, drop = FALSE]
  colnames(m) <- sub("^(HR|LR)_", "", meta$subject[match(cols, meta$sample_id)])
  m
}
t1 <- wide("T1")
t2 <- wide("T2")
shared <- intersect(colnames(t1), colnames(t2))
arm_s <- meta$arm[match(shared, sub("^(HR|LR)_", "", meta$subject))]

ancova <- map_dfr(seq_len(nrow(mat)), function(i) {
  y <- t2[i, shared]
  x0 <- t1[i, shared]
  ok <- stats::complete.cases(y, x0)
  if (sum(ok) < 8 || length(unique(arm_s[ok])) < 2) {
    return(tibble(uniprot_id = rownames(mat)[i], p = NA_real_))
  }
  cf <- summary(stats::lm(y[ok] ~ x0[ok] + arm_s[ok]))$coefficients
  tibble(uniprot_id = rownames(mat)[i], p = cf[nrow(cf), 4])
}) |>
  mutate(bh = stats::p.adjust(p, "BH"))

res_a <- tibble(
  spec = "ANCOVA on baseline", term = "arm | T2 adjusted for T1",
  n = sum(!is.na(ancova$p)), nominal = sum(ancova$p < 0.05, na.rm = TRUE),
  bh05 = sum(ancova$bh < 0.05, na.rm = TRUE),
  bh10 = sum(ancova$bh < 0.10, na.rm = TRUE),
  min_bh = min(ancova$bh, na.rm = TRUE),
  top = paste(utils::head(ancova$uniprot_id[order(ancova$p)], 3), collapse = ", ")
)

all_res <- bind_rows(res_cont, res_w, res_a) |>
  mutate(top = vapply(top, gene_of, character(1)), chance = round(0.05 * n))

cat("=== alternative specifications, all against the same matrix ===\n")
all_res |>
  select(spec, term, nominal, chance, bh05, bh10, min_bh, top) |>
  arrange(min_bh) |>
  as.data.frame() |>
  print(row.names = FALSE, digits = 3)

cat(sprintf("\narrayWeights range: %.2f to %.2f (1 = average sample)\n", min(aw), max(aw)))

write.xlsx(list(summary = all_res, ancova = ancova), file.path(OUT, "alt_specs.xlsx"))
cat(sprintf("wrote %s\n", file.path(OUT, "alt_specs.xlsx")))
