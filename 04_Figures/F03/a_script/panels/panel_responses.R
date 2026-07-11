# Responses: within-group HR/LR trajectories. Row-major fill with ncol = 2 stacks
# training in the left column and acute in the right, HR over LR. Red/blue ramp.
panel_responses <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = c("Training_HR", "Acute_HR", "Training_LR", "Acute_LR"),
    palette = RING_PALETTES$responses, ncol = 2
  )
}
