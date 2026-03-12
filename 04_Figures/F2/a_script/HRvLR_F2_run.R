#!/usr/bin/env Rscript
# HRvLR_F2_run.R — Run all F2 panels
# Run from project root: Rscript 04_Figures/F2/a_script/HRvLR_F2_run.R
#
# Panels:
#   A–G : Volcano + fGSEA ring (7 contrasts; Hallmark + GO:BP + CP databases)
#   H   : Training scatter (Training_HR vs Training_LR logFC)
#   I   : Acute scatter (Acute_HR vs Acute_LR logFC)

here::i_am(".here")

panels <- c("panel_A", "panel_B", "panel_C", "panel_D",
            "panel_E", "panel_F", "panel_G", "panel_H", "panel_I")

t0 <- proc.time()
for (p in panels) {
  message("\n=== Running ", p, " ===")
  source(here::here("04_Figures/F2/a_script/", paste0(p, ".R")))
}
elapsed <- proc.time() - t0
message(sprintf("\nAll F2 panels complete (%.0fs)", elapsed["elapsed"]))
