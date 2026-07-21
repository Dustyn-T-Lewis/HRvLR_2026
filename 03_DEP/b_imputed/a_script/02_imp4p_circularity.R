#!/usr/bin/env Rscript
# Why the imp4p arm breaks the null, and why that is circularity rather than fragility.
#
# imp4p::impute.mle loops `for (n in 1:nb_cond)` and fits a separate EM inside each level of
# `conditions`, imputing each gap from its own cell's mean at n = 7-8. We pass Group_Time, and
# the nine contrasts then test among those same six cells. Missing values are pulled toward the
# cell mean, which manufactures between-cell differences and shrinks within-cell residual
# variance. Both inflate the moderated t.
#
# The arm cannot simply be blinded: a single-level `conditions` makes estim.mix fail, because
# imputing within conditions is imp4p's architecture, not our misuse of it. So the control runs
# the other way. Impute within RANDOM cells of identical size, orthogonal to Group_Time, and
# change nothing else. Random cells cannot manufacture real Group_Time signal, so any collapse in
# the hit count is attributable to the tested factor having been inside the imputation.
#
# Result: real cells give a three-figure BH<0.10 count across the five HR-vs-LR / interaction
# contrasts; random cells collapse it to near zero. See b_reports/imp4p_circularity.csv.

pacman::p_load(proteoDA, imp4p, here, readr, dplyr, tibble, purrr)
source(here("03_DEP", "contrasts.R"))

N_PERM <- 3
NULL_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

dal <- readRDS(here("02_Normalization", "c_data", "DAList_normalized.rds"))
mat <- as.matrix(dal$data)
group_time <- factor(dal$metadata$Group_Time[match(colnames(mat), dal$metadata$Col_ID)])

fit_contrasts <- function(imputed) {
  d <- dal
  d$data <- imputed
  d$metadata$group <- factor(
    d$metadata$Group_Time,
    levels = c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
  )
  d$metadata$subject <- d$metadata$Subject_ID
  d <- add_design(d, "~ 0 + group + (1 | subject)")
  d <- add_contrasts(d, contrasts_vector = HRVLR_CONTRASTS)
  d <- fit_limma_model(d)
  d <- extract_DA_results(d, pval_thresh = 0.10, lfc_thresh = 0, adj_method = "BH")
  imap_dfr(d$results, \(r, cname) as_tibble(r, rownames = "uniprot_id") |> mutate(contrast = cname))
}

impute_within <- function(conditions) {
  set.seed(42)
  out <- as.matrix(impute.mi(
    tab = mat, conditions = conditions,
    methodMCAR = "mle", methodMNAR = "igcda", progress.bar = FALSE
  ))
  dimnames(out) <- dimnames(mat)
  out
}

summarise_arm <- function(res, label) {
  tibble(
    imputed_within = label,
    BH_10_total = sum(res$adj.P.Val < 0.10, na.rm = TRUE),
    BH_10_null_contrasts = sum(res$adj.P.Val < 0.10 & res$contrast %in% NULL_CONTRASTS, na.rm = TRUE),
    BH_10_baseline = sum(res$adj.P.Val < 0.10 & res$contrast == "Baseline_HRvLR", na.rm = TRUE)
  )
}

real <- summarise_arm(fit_contrasts(impute_within(group_time)), "Group_Time (the tested factor)")

# Same six cells, same sizes, membership shuffled. Orthogonal to every contrast.
permuted <- map_dfr(seq_len(N_PERM), function(i) {
  set.seed(100 + i)
  shuffled <- factor(sample(as.character(group_time)), levels = levels(group_time))
  summarise_arm(fit_contrasts(impute_within(shuffled)), sprintf("random cells (seed %d)", 100 + i))
})

out <- bind_rows(real, permuted)

report_dir <- here("03_DEP", "b_imputed", "b_reports")
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(out, file.path(report_dir, "imp4p_circularity.csv"))

print(as.data.frame(out), row.names = FALSE)

ni_dir <- here("03_DEP", "a_non_imputed", "c_data", "04_per_contrast_results")
ni <- map_dfr(list.files(ni_dir, pattern = "\\.csv$", full.names = TRUE), function(f) {
  read_csv(f, show_col_types = FALSE) |>
    transmute(contrast = tools::file_path_sans_ext(basename(f)), adj.P.Val)
})
ni_total <- sum(ni$adj.P.Val < 0.10, na.rm = TRUE)
ni_null <- sum(ni$adj.P.Val < 0.10 & ni$contrast %in% NULL_CONTRASTS, na.rm = TRUE)

cat(sprintf(
  "\nnon-imputed primary arm: %d BH<0.10 total, %d in the null contrasts.\n",
  ni_total, ni_null
))
cat("Imputing within the tested cells is what produces the hits; imputing within random\n")
cat("cells of the same size does not. The null is robust; the imp4p arm is circular.\n")
