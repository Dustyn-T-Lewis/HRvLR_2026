# Permutation control for the Pi-score selection, the counterpart to
# 02_imp4p_circularity.R.
#
# 01_run_dep.R selects by Pi < 0.05 and the manuscript quotes the resulting
# counts. Pi carries no error rate, so those counts have no reference until one
# is built: this refits the nine contrasts under permuted HR/LR labels and asks
# how many Pi hits the design produces when there is nothing to find.

pacman::p_load(here, dplyr, openxlsx)
source(here("functions", "pi_permutation.R"))
source(here("functions", "feature_levels.R"))

N_PERM <- 200

res <- pi_permutation_null(protein_matrix(), n_perm = N_PERM, seed = 42)

print(as.data.frame(res), digits = 3)
cat(sprintf(
  "\nObserved total across the 5 HR-vs-LR / interaction contrasts: %d\n",
  sum(res$observed[grepl("HRvLR|Interaction", res$contrast)])
))
cat(sprintf(
  "Permuted median for the same five: %.0f\n",
  sum(res$null_median[grepl("HRvLR|Interaction", res$contrast)])
))

write.xlsx(
  list(pi_permutation = res),
  here("03_DEP", "a_non_imputed", "c_data", "06_pi_permutation.xlsx")
)
