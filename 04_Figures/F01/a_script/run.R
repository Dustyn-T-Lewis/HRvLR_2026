# F01 build: render every panel, stitch the composite, write one workbook.
pacman::p_load(here, patchwork, ggplot2, dplyr, tidyr, tibble, openxlsx)

F01_PANELS <- list()
F01_AUDIT <- list()
source(here("04_Figures", "F01", "a_script", "setup.R"))

unlink(setdiff(
  list.files(F01_RPT, full.names = TRUE, recursive = TRUE),
  file.path(F01_RPT, ".gitkeep")
))
dir.create(file.path(F01_RPT, "panels"), showWarnings = FALSE)
dir.create(file.path(F01_RPT, "supp"), showWarnings = FALSE)

panel_dir <- here("04_Figures", "F01", "a_script", "panels")
for (f in c("panel_a_volume", "panel_b_composite", "panel_c_continuum", "panel_d_forest")) {
  source(file.path(panel_dir, paste0(f, ".R")))
}

supp_dir <- here("04_Figures", "F01", "a_script", "supp")
source(file.path(supp_dir, "pheno", "pheno_supp.R"))
source(file.path(supp_dir, "hlm", "hlm_supp.R"))

composite <- (F01_PANELS$volume | F01_PANELS$composite) /
  (F01_PANELS$continuum | F01_PANELS$forest) +
  plot_layout(heights = c(0.85, 1.15)) +
  plot_annotation(
    title = "F01 · Phenotype atlas: matched work, divergent growth",
    subtitle = "Equal work and strength gains; groups split by muscle growth.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "grey30")
    )
  )

ggsave(file.path(F01_RPT, "F01_composite.png"), composite,
  width = 300, height = 250, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(F01_RPT, "F01_composite.pdf"), composite,
  width = 300, height = 250, units = "mm", device = PDF_DEVICE, bg = "white"
)

summary_by_measure <- F01_AUDIT$secondary_delta |>
  left_join(
    F01_AUDIT$lmm_interaction |>
      select(measure, domain,
        lmm_advantage = estimate, lmm_ci_lo = ci_lo, lmm_ci_hi = ci_hi, lmm_p = p
      ),
    by = "measure"
  ) |>
  left_join(F01_AUDIT$lmm_fit_summary |> select(measure, icc), by = "measure") |>
  select(measure, domain,
    group = Group, n, delta_mean, delta_sem,
    lmm_advantage, lmm_ci_lo, lmm_ci_hi, lmm_p, icc
  )

controls <- bind_rows(
  F01_AUDIT$volume_load |> mutate(measure = "Accumulated volume load", .before = 1),
  F01_AUDIT$composite |> mutate(measure = "Composite hypertrophy", .before = 1)
) |>
  select(measure, group = Group, n, mean, sd, sem)

metadata <- tibble::tribble(
  ~field, ~value,
  "figure", "F01 phenotype atlas (HRvLR)",
  "design", "16 subjects (8 HR, 8 LR), repeated measures T1/T2 (T3 not used here)",
  "groups", "HR/LR defined by the composite hypertrophy score; split is descriptive",
  "model", "lmer(value ~ Group*Timepoint + (1|subject)), REML, Satterthwaite df",
  "divergence", "Group x Timepoint interaction; standardized as HR-minus-LR change advantage",
  "control", "accumulated volume load compared by effect size + TOST equivalence",
  "source", "00_input/HRvLR_meta.csv"
)

wb <- createWorkbook()
addWorksheet(wb, "summary_by_measure")
writeData(wb, "summary_by_measure", summary_by_measure)
addWorksheet(wb, "controls_and_axis")
writeData(wb, "controls_and_axis", controls, startRow = 1)
writeData(wb, "controls_and_axis", F01_AUDIT$responder_axis, startRow = nrow(controls) + 3)
addWorksheet(wb, "hlm_fit_summary")
writeData(wb, "hlm_fit_summary", F01_AUDIT$lmm_fit_summary)
addWorksheet(wb, "metadata")
writeData(wb, "metadata", metadata)
saveWorkbook(wb, file.path(F01_DAT, "F01_source_data.xlsx"), overwrite = TRUE)

cat("F01 rebuilt: composite, panels, supplement, workbook written\n")
