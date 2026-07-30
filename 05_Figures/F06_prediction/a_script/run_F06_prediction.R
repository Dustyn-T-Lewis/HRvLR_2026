suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_pred_leaf.R"))
}))

# args: [1] comma-separated levels (default all), [2] max B (default 1000).
# Fast levels sweep 0/200/1000; proteins run in the background at 0/200. Each
# leaf sweeps the six adaptation deltas internally.
args <- commandArgs(trailingOnly = TRUE)
levels_run <- if (length(args) >= 1 && nzchar(args[1])) {
  strsplit(args[1], ",")[[1]]
} else {
  SWEEP_LEVELS
}
bmax <- if (length(args) >= 2) as.integer(args[2]) else 1000L
b_grid <- if (bmax >= 1000L) c(0L, 200L, 1000L) else c(0L, bmax)

bundle <- pred_load()
fingerprint <- sweep_fingerprint(bundle, b_grid)
grid <- do.call(rbind, lapply(levels_run, function(lv) {
  expand.grid(
    level = lv, config = SWEEP_CONFIGS,
    method = methods_for_level(CONT_METHODS, lv),
    stringsAsFactors = FALSE
  )
}))
# plain needs p < n; the trajectory encoding triples the module count past n, so
# the plain model is undefined there and is dropped.
grid <- grid[!(grid$method == "plain" & grid$config == "trajectory"), ]

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  if (leaf_done("F06_prediction", g$level, g$config, g$method, fingerprint)) {
    message(sprintf(
      "[cont %d/%d] %s | %s | %s  skip (done)",
      i, nrow(grid), g$level, g$config, g$method
    ))
    next
  }
  t0 <- Sys.time()
  s <- suppressWarnings(
    run_cont_sweep_leaf(bundle, g$level, g$config, g$method, b_grid = b_grid)
  )
  top <- s[s$B == max(s$B), ]
  best <- top[which.max(top$q2), ]
  message(sprintf(
    "[cont %d/%d] %s | %s | %s  best Q2=%.2f (%s) p=%.3f  (%.0fs)",
    i, nrow(grid), g$level, g$config, g$method,
    best$q2, best$outcome, best$perm_p_q2,
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
}
message("continuous run done for: ", paste(levels_run, collapse = ", "))
