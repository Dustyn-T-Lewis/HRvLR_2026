# F01 supplement (pheno): change score by group across every outcome, from the
# same stat machinery as the divergence forest.
if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))

secondary <- build_secondary(meta)
rpt_pheno <- file.path(F01_RPT, "supp", "pheno")
dir.create(rpt_pheno, recursive = TRUE, showWarnings = FALSE)
save_panel(secondary$plot, file.path(rpt_pheno, "F01_pheno_supp"), 230, 150)
F01_AUDIT[["secondary_delta"]] <- secondary$audit
