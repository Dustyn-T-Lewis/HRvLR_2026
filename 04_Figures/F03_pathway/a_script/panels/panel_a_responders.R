# Panel A: the five contrasts that actually test the study's question - does the
# proteome of a high responder differ from a low responder? Three between-responder
# rings (baseline, trained, acute; green/purple) beside the two interaction contrasts
# (orange/purple). The interactions are the lowest-powered cells (8 vs 8) and are read
# as hypothesis-generating, not as findings.
#
# This panel owns its supplement: the un-deduplicated top-30 audit for its own five
# contrasts, so the rings' parsimony can be checked against the full ranked enrichment.
if (!exists("fg")) source(here::here("04_Figures", "F03_pathway", "a_script", "setup.R"))

RESPONDER_CONTRASTS <- c(
  "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
  "Training_Interaction", "Acute_Interaction"
)

panel_responders <- function(fg, dep, pw) {
  diff <- RING_PALETTES$differential
  intr <- RING_PALETTES$interaction
  specs <- list(
    list(contrast = "Baseline_HRvLR", palette = diff, tag = "A"),
    list(contrast = "Trained_HRvLR", palette = diff, tag = "B"),
    list(contrast = "Acute_HRvLR", palette = diff, tag = "C"),
    list(contrast = "Training_Interaction", palette = intr, tag = "D"),
    list(contrast = "Acute_Interaction", palette = intr, tag = "E")
  )
  build_ring_grid(fg, dep, pw, specs, ncol = 3, byrow = TRUE, guide_area = TRUE)
}

responders <- panel_responders(fg, dep, pw)
F03_PANELS[["responders"]] <- responders$plot
F03_REPORTS[["responders"]] <- responders$reports

responders_top30 <- panel_top30(fg, RESPONDER_CONTRASTS)
F03_SUPP[["responders_top30"]] <- responders_top30$plots
F03_TABLES[["responders_top30"]] <- responders_top30$table
