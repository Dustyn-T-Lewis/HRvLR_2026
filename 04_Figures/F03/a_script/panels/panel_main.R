# Main composite: the four within-group responses (red/blue) beside the two
# interaction contrasts (orange/purple), six rings in three columns -- training,
# acute, interaction. Tags run down each column (A/B, C/D, E/F) so a column reads
# as one comparison type. Interactions are the lowest-powered cells (n = 8 vs 8),
# read as hypothesis-generating.
panel_main <- function(fg, dep, pw) {
  resp <- RING_PALETTES$responses
  intr <- RING_PALETTES$interaction
  specs <- list(
    list(contrast = "Training_HR", palette = resp, tag = "A"),
    list(contrast = "Training_LR", palette = resp, tag = "B"),
    list(contrast = "Acute_HR", palette = resp, tag = "C"),
    list(contrast = "Acute_LR", palette = resp, tag = "D"),
    list(contrast = "Training_Interaction", palette = intr, tag = "E"),
    list(contrast = "Acute_Interaction", palette = intr, tag = "F")
  )
  build_ring_grid(fg, dep, pw, specs, ncol = 3, byrow = FALSE)
}
