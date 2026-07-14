# F04 composite: WGCNA module-trait atlas (A) beside the module-level NES
# scatters (B), stitched with cowplot::ggdraw so Panel A overlaps and hides
# Panel B's left margin. Then the source-data workbook.
pacman::p_load(here, ggplot2, cowplot, png, grid, readr, openxlsx)
source(here("04_Figures", "functions", "style.R"))

RPT_DIR <- here("04_Figures", "F04", "b_reports")
DAT_DIR <- here("04_Figures", "F04", "c_data")
PANEL_DIR <- file.path(RPT_DIR, "panels")

read_panel <- function(name) {
  path <- file.path(PANEL_DIR, paste0(name, ".png"))
  if (!file.exists(path)) stop("Missing panel: ", path, " (run HRvLR_F04_run.R first)")
  rasterGrob(readPNG(path), interpolate = TRUE)
}
pA_grob <- read_panel("panel_a_wgcna_map")
pB_grob <- read_panel("panel_b_nes_scatters")
pB_leg_grob <- read_panel("panel_b_nes_legend")

COMP_W <- 470
COMP_H <- 300
mm2x <- function(mm) mm / COMP_W
mm2y <- function(mm) mm / COMP_H

crop_l <- 0.01
crop_r <- 0.80
crop_b <- 0.17
crop_t <- 0.85
save_w <- COMP_W * (crop_r - crop_l)
save_h <- COMP_H * (crop_t - crop_b)

grid_top <- 0.80
grid_bot <- 0.22
b_x <- 0.43
b_w <- 0.55
b_grid_h <- grid_top - grid_bot

tag_sz <- 14
title_sz <- 12
subtitle_sz <- 9

composite <- ggdraw(xlim = c(crop_l, crop_r), ylim = c(crop_b, crop_t)) +
  theme(plot.background = element_rect(fill = "white", color = NA)) +
  draw_grob(pB_grob,
    x = b_x, y = grid_bot, width = b_w, height = b_grid_h,
    hjust = 0, vjust = 0
  ) +
  draw_grob(pB_leg_grob,
    x = 0.70, y = 0.195, width = 0.24, height = 0.04,
    hjust = 0.5, vjust = 0.5
  ) +
  draw_grob(pA_grob,
    x = 0.01, y = 0, width = 0.60, height = 0.96,
    hjust = 0, vjust = 0
  ) +
  draw_label("A",
    x = mm2x(15), y = mm2y(248),
    size = tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label("WGCNA Module–Trait Associations",
    x = mm2x(58), y = mm2y(248),
    size = title_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label("Baseline (T1) → adaptation | Responder (HR − LR) per timepoint",
    x = mm2x(58), y = mm2y(243),
    size = subtitle_sz, fontface = "bold.italic", colour = "grey40",
    hjust = 0, vjust = 1
  ) +
  draw_label("B",
    x = mm2x(285), y = mm2y(248),
    size = tag_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label("Module–Level NES Scatters",
    x = mm2x(294), y = mm2y(248),
    size = title_sz, fontface = "bold", hjust = 0, vjust = 1
  ) +
  draw_label("fgsea on module-member moderated-t ranks",
    x = mm2x(294), y = mm2y(243),
    size = subtitle_sz, fontface = "bold.italic", colour = "grey40",
    hjust = 0, vjust = 1
  )

ggsave(file.path(RPT_DIR, "F04_composite.png"), composite,
  width = save_w, height = save_h, units = "mm", dpi = 300, bg = "white"
)
ggsave(file.path(RPT_DIR, "F04_composite.pdf"), composite,
  width = save_w, height = save_h, units = "mm", device = PDF_DEVICE, bg = "white"
)

sheet_desc <- c(
  wgcna_membership = "Protein-to-module assignment (WGCNA, full imputed proteome)",
  wgcna_eigengene = "Module eigengenes per sample and timepoint",
  module_prediction = "Baseline (T1) eigengene LOO-CV prediction of adaptation traits",
  module_responder = "Per-module responder (HR vs LR) signal per timepoint",
  module_ora = "Over-representation per module (supplement ORA bars)",
  panel_b_module_fgsea = "Panel B: module-level fgsea NES across contrasts",
  module_phenotype = "Supplement: module eigengene delta vs response trait",
  module_phenotype_lmm = "Supplement: module-phenotype linear mixed model",
  phenotype = "Per-subject phenotype table (T2 composite + T1->T2 deltas)"
)
present <- names(sheet_desc)[file.exists(
  file.path(DAT_DIR, paste0(names(sheet_desc), ".csv"))
)]
overview <- data.frame(
  sheet = substr(present, 1, 31), description = unname(sheet_desc[present]),
  stringsAsFactors = FALSE
)
wb <- createWorkbook()
addWorksheet(wb, "overview")
writeData(wb, "overview", overview)
for (s in present) {
  addWorksheet(wb, substr(s, 1, 31))
  writeData(wb, substr(s, 1, 31), read_csv(
    file.path(DAT_DIR, paste0(s, ".csv")),
    show_col_types = FALSE
  ))
}
saveWorkbook(wb, file.path(DAT_DIR, "F04_source_data.xlsx"), overwrite = TRUE)
message("F04 composite + workbook saved to ", RPT_DIR)
