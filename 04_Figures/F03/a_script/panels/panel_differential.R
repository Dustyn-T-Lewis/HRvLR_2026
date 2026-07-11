# Differential: the untrained baseline difference between responders (HR vs LR at
# T1). The trained- and acute-state differences move to the supplement, since the
# interaction contrasts carry the same divergence with the timepoint controlled.
# Orange (up) / purple (down) NES ramp.
panel_differential <- function(fg, dep, pw) {
  render_family(
    fg, dep, pw,
    contrasts = "Baseline_HRvLR",
    palette = RING_PALETTES$differential, ncol = 1
  )
}
