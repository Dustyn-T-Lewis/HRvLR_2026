# Module-level tests: the nine contrasts, classification and phenotype
# association, all on the eigengene matrix.
#
# The contrasts reuse the protein fit's design and subject block, with the
# within-subject correlation re-estimated on the eigengenes, so a module is
# tested exactly as a protein was. Twelve features make BH far weaker a filter
# here than over 1900 proteins; the chance tables are read with that in mind.

pacman::p_load(here, dplyr, tibble, purrr, limma)

source(here("functions", "contrasts.R"))
source(here("functions", "classify.R"))
source(here("functions", "shared_utils.R"))

OUT_DIR <- here("04_Networks", "03_Classify_Associate", "c_data")

modules <- readRDS(here("04_Networks", "01_Modules", "c_data", "modules.rds"))
fit <- readRDS(here(
  "02_Proteins", "01_Differential", "c_data", "01_limma_DAList.rds"
))
me <- modules$eigengenes
design <- fit$design$design_matrix
subject <- fit$metadata$subject
stopifnot(identical(colnames(me), rownames(design)))

correlation <- duplicateCorrelation(me, design, block = subject)$consensus
me_fit <- lmFit(me, design, block = subject, correlation = correlation) |>
  contrasts.fit(fit$design$contrast_matrix) |>
  eBayes()
contrast_res <- map(CONTRAST_NAMES, function(ct) {
  topTable(me_fit, coef = ct, number = Inf, sort.by = "none") |>
    as_tibble(rownames = "module") |>
    transmute(
      contrast = ct, module = .data$module, logFC = .data$logFC,
      t = .data$t, p = .data$P.Value, bh = .data$adj.P.Val
    )
}) |>
  list_rbind()

screens <- c(
  list(contrasts = contrast_res, correlation = tibble(consensus = correlation)),
  screen_level(me, modules$meta)
)

clear_dir(OUT_DIR)
saveRDS(screens, file.path(OUT_DIR, "module_screens.rds"))
openxlsx::write.xlsx(screens, file.path(OUT_DIR, "module_screens.xlsx"))

message(sprintf("eigengene within-subject correlation: %.3f", correlation))
print(
  filter(contrast_res, .data$p < 0.05) |> arrange(.data$p),
  n = Inf
)
