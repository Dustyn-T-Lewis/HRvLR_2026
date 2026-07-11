# Supplement: the trained- and acute-state between-responder differences
# (Trained_HRvLR, Acute_HRvLR), demoted from the main composite because the
# interaction contrasts show the same divergence with the timepoint controlled.
# Same orange/purple ramp as the baseline differential panel.
panel_between_state <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = c("Trained_HRvLR", "Acute_HRvLR"),
    palette = RING_PALETTES$differential, ncol = 2
  )
}
