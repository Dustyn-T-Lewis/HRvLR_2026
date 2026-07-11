# Interaction: where the two groups' responses diverge, over training and over
# the acute bout. These carry the central divergence claim but sit at the lowest-
# powered cells of an unbalanced 16-subject design, so the narrative reads them as
# exploratory. Blue (down) / red (up) NES ramp.
panel_interaction <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = c("Training_Interaction", "Acute_Interaction"),
    palette = RING_PALETTES$interaction, ncol = 2
  )
}
