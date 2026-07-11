# F01 supplement (hlm): one performance::check_model exemplar (mCSA) and the
# per-model fit table. The six models share structure, so one figure shows the
# assumption pattern; compare_performance validates all six quantitatively (R2,
# ICC, singular flag), and the leave-one-subject table is in F01_lmm_report.qmd.
if (!exists("meta")) source(here::here("04_Figures", "F01", "a_script", "setup.R"))
pacman::p_load(performance, see, patchwork)

models <- lmm_models(meta)
rpt_hlm <- file.path(F01_RPT, "supp", "hlm")
dir.create(rpt_hlm, recursive = TRUE, showWarnings = FALSE)

ggsave(
  file.path(rpt_hlm, "F01_hlm_check_mCSA.png"),
  plot(check_model(models[["mCSA"]], verbose = FALSE)),
  width = 280, height = 360, units = "mm", dpi = 150, bg = "white" # tall so y-titles clear the subtitles
)

F01_AUDIT[["lmm_fit_summary"]] <- as.data.frame(
  performance::compare_performance(models, metrics = c("AIC", "BIC", "R2", "ICC", "RMSE"))
) |>
  dplyr::transmute(
    measure = Name, AIC, BIC,
    r2_conditional = R2_conditional, r2_marginal = R2_marginal, icc = ICC, RMSE
  )

cat("F01 hlm check_model exemplar + fit table written\n")
