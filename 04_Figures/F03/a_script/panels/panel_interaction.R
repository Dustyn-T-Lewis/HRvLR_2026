# Interaction (training over acute): where the two groups' responses diverge.
# Its own stacked column beside the responses. Exploratory at n = 8 vs 8, two
# partial subjects. Orange/purple ramp.
panel_interaction <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = c("Training_Interaction", "Acute_Interaction"),
    palette = RING_PALETTES$interaction, ncol = 1
  )
}
