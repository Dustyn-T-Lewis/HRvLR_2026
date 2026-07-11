# Responses: within-group HR and LR trajectories over training (T2-T1) and the
# acute bout (T3-T2). Green (up) / purple (down) NES ramp.
panel_responses <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = c("Training_HR", "Training_LR", "Acute_HR", "Acute_LR"),
    palette = RING_PALETTES$responses, ncol = 4
  )
}
