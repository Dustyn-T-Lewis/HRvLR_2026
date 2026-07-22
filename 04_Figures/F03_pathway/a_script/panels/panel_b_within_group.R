# Lead figure: the four within-group responses - what training and the acute bout do to
# the proteome inside each arm. Arms in columns (HR left, LR right), phases in rows, with a
# faint arm-coloured wash so the two arms read as columns. Seeing HR beside LR is the direct
# view of the responder question; the formal between-responder test is the second figure.
if (!exists("fg")) source(here::here("04_Figures", "F03_pathway", "a_script", "setup.R"))

WITHIN_GROUP_CONTRASTS <- c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR")

panel_within_group <- function(fg, dep, pw) {
  resp <- RING_PALETTES$responses
  hr_bg <- scales::alpha(GROUP_COLORS[["HR"]], 0.05)
  lr_bg <- scales::alpha(GROUP_COLORS[["LR"]], 0.05)
  specs <- list(
    list(contrast = "Training_HR", palette = resp, tag = "A", bg = hr_bg),
    list(contrast = "Training_LR", palette = resp, tag = "B", bg = lr_bg),
    list(contrast = "Acute_HR", palette = resp, tag = "C", bg = hr_bg),
    list(contrast = "Acute_LR", palette = resp, tag = "D", bg = lr_bg)
  )
  build_ring_grid(fg, dep, pw, specs, ncol = 2, byrow = TRUE, legend = "right")
}

within_group <- panel_within_group(fg, dep, pw)
F03_PANELS[["within_group"]] <- within_group$plot
F03_REPORTS[["within_group"]] <- within_group$reports

within_top30 <- panel_top30(fg, WITHIN_GROUP_CONTRASTS)
F03_SUPP[["within_group_top30"]] <- within_top30$plots
F03_TABLES[["within_group_top30"]] <- within_top30$table
