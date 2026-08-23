# Panel A: the two pooled contrasts, Training and Acute, all 16 subjects,
# no arm split anywhere. Mirrors
# categorical/F03_pathway/a_script/panels/panel_b_within_group.R's neutral
# "responses" palette (the same one categorical uses for its own
# no-group-contrast rings) rather than the differential/interaction
# palettes, which encode a between-arm comparison this tree doesn't make.
#
# This panel owns its supplement: the un-deduplicated top-30 audit for
# both contrasts, so the rings' parsimony can be checked against the full
# ranked enrichment.
if (!exists("fg")) {
  source(here::here(
    "03_Analysis", "continuous", "F03_pathway", "a_script", "setup.R"
  ))
}

panel_pooled <- function(fg, dep, pw) {
  resp <- RING_PALETTES$responses
  specs <- list(
    list(contrast = "Training", palette = resp, tag = "A"),
    list(contrast = "Acute", palette = resp, tag = "B")
  )
  build_ring_grid(fg, dep, pw, specs, ncol = 2, byrow = TRUE, legend = "right")
}

pooled <- panel_pooled(fg, dep, pw)
F03_PANELS[["pooled"]] <- pooled$plot
F03_REPORTS[["pooled"]] <- pooled$reports

pooled_top30 <- panel_top30(fg, CONTRASTS)
F03_SUPP[["pooled_top30"]] <- pooled_top30$plots
F03_TABLES[["pooled_top30"]] <- pooled_top30$table
