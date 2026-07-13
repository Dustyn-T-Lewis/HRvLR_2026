# Between-responder composite, its own figure: where HR and LR differ at baseline
# (T1), trained (T2), and acute (T3). Three rings in one row, green/purple.
panel_hrvlr <- function(fg, dep, pw) {
  diff <- RING_PALETTES$differential
  specs <- list(
    list(contrast = "Baseline_HRvLR", palette = diff, tag = "A"),
    list(contrast = "Trained_HRvLR", palette = diff, tag = "B"),
    list(contrast = "Acute_HRvLR", palette = diff, tag = "C")
  )
  build_ring_grid(fg, dep, pw, specs, ncol = 3, byrow = TRUE)
}
