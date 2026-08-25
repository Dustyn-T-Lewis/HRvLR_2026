# Does the refit reproduce V1's committed numbers?
#
# Stage 01 rebuilds the nine contrasts rather than copying V1's workbook, so
# something has to prove the rebuild did not move them. This joins on uniprot
# id per contrast and reports the largest absolute discrepancy in logFC, raw p
# and BH.
#
# A pipeline that cannot reproduce its predecessor's numbers has either found a
# bug or introduced one, and this is the script that says which.

pacman::p_load(here, dplyr, tibble, readr, openxlsx)

source(here("03_Features", "contrasts.R"))

OUT_DIR <- here("03_Features", "01_DEP", "b_reports")
TOL <- 1e-6

V1_BOOK <- file.path(
  dirname(here()), "A_HRvLR_2026", "03_Analysis", "categorical",
  "01_Proteins", "c_data", "05_results.xlsx"
)

if (!file.exists(V1_BOOK)) {
  message("V1 workbook not found; skipping equivalence check")
  quit(save = "no")
}

fitted <- read_csv(
  here("03_Features", "01_DEP", "c_data", "01_dep_results.csv"),
  show_col_types = FALSE
)

equivalence <- bind_rows(lapply(CONTRAST_NAMES, function(ct) {
  ref <- openxlsx::read.xlsx(V1_BOOK, ct)
  j <- inner_join(
    fitted |> filter(.data$contrast == ct),
    ref |> transmute(
      uniprot_id = .data$uniprot_id, ref_lfc = .data$logFC,
      ref_p = .data$P.Value, ref_bh = .data$adj.P.Val
    ),
    by = "uniprot_id"
  )
  worst <- max(
    abs(j$logFC - j$ref_lfc), abs(j$P.Value - j$ref_p),
    abs(j$adj.P.Val - j$ref_bh),
    na.rm = TRUE
  )
  tibble(
    contrast = ct, n_matched = nrow(j),
    max_abs_diff = worst, equivalent = worst < TOL
  )
}))

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
write_csv(equivalence, file.path(OUT_DIR, "v1_equivalence.csv"))

print(as.data.frame(equivalence), row.names = FALSE, digits = 3)
if (all(equivalence$equivalent)) {
  message("\nall nine contrasts reproduce V1 to within ", TOL)
} else {
  stop(
    "contrasts differ from V1: ",
    paste(equivalence$contrast[!equivalence$equivalent], collapse = ", ")
  )
}
