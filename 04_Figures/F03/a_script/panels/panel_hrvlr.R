# Between-responder differences across states, its own composite: where HR and LR
# differ at baseline (T1), trained (T2), and acute (T3). Green/purple ramp.
panel_hrvlr <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = c("Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR"),
    palette = RING_PALETTES$differential, ncol = 3
  )
}
