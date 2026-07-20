# Panel B (supplement): the four within-group responses - what training and the acute
# bout do to the proteome inside each arm. These do not test responder status, so they
# do not belong in the main figure. They are the positive control: they show the
# proteome does respond, which is what makes the null in Panel A meaningful rather than
# merely underpowered.
#
# This panel owns its supplement: the top-30 audit for its own four contrasts.
if (!exists("fg")) source(here::here("04_Figures", "F03_pathway", "a_script", "setup.R"))

WITHIN_GROUP_CONTRASTS <- c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR")

panel_within_group <- function(fg, dep, pw) {
  resp <- RING_PALETTES$responses
  specs <- list(
    list(contrast = "Training_HR", palette = resp, tag = "A"),
    list(contrast = "Training_LR", palette = resp, tag = "B"),
    list(contrast = "Acute_HR", palette = resp, tag = "C"),
    list(contrast = "Acute_LR", palette = resp, tag = "D")
  )
  build_ring_grid(fg, dep, pw, specs, ncol = 2, byrow = FALSE)
}

within_group <- panel_within_group(fg, dep, pw)
F03_PANELS[["within_group"]] <- within_group$plot
F03_REPORTS[["within_group"]] <- within_group$reports

within_top30 <- panel_top30(fg, WITHIN_GROUP_CONTRASTS)
F03_SUPP[["within_group_top30"]] <- within_top30$plots
F03_TABLES[["within_group_top30"]] <- within_top30$table
