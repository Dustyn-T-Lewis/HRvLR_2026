pacman::p_load(here, dplyr, readr, tibble, purrr)
fdir <- here("04_Figures", "05_test", "modules", "a_script", "functions")
for (f in c("fit.R", "coupling.R", "prediction.R", "honest_refit.R")) source(file.path(fdir, f))

da <- readRDS(here("02_Normalization", "imputation", "c_data", "DAList_imputed_missforest.rds"))
abund <- as.matrix(da$data)
sample_meta <- da$metadata |>
  transmute(
    sample_id = Col_ID, subject = Subject_ID, timepoint = Timepoint,
    group = Group
  )
pheno <- read_csv(here("00_input", "c_data", "phenotype.csv"), show_col_types = FALSE)

fit <- fit_modules(abund, sample_meta)
out <- here("04_Figures", "05_test", "modules", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_csv(fit$membership, file.path(out, "module_assignments.csv"))
fit$eigengenes |>
  tibble::rownames_to_column("sample_id") |>
  write_csv(file.path(out, "eigengenes.csv"))

bind_rows(lapply(c("baseline", "training", "acute"), function(p) {
  module_coupling(fit$eigengenes, sample_meta, pheno, p)
})) |>
  write_csv(file.path(out, "coupling.csv"))

hr_lr <- setNames(
  as.integer(sample_meta$group[match(unique(sample_meta$subject), sample_meta$subject)] == "HR"),
  unique(sample_meta$subject)
)
for (phase in c("training", "acute")) {
  run_honest_refit(abund, sample_meta, hr_lr, phase) |>
    write_csv(file.path(out, paste0("loso_refit_", phase, ".csv")))
}
cat("05_test/modules engine written\n")
