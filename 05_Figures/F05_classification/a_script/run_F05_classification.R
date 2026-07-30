suppressWarnings(suppressMessages({
  library(here)
  source(here("functions", "sweep_pred_leaf.R"))
}))

# args: [1] comma-separated levels (default all), [2] max B (default 1000).
# The B grid is 0/200/1000 for the fast levels and 0/200 for proteins, set by
# the caller so the long protein pole runs in the background at B = 200.
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
    method = methods_for_level(CLASS_METHODS, lv),
    stringsAsFactors = FALSE
  )
}))
# plain needs p < n; the trajectory encoding triples the module count past n, so
# the plain model is undefined there and is dropped.
grid <- grid[!(grid$method == "plain" & grid$config == "trajectory"), ]

for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  done <- leaf_done(
    "F05_classification", g$level, g$config, g$method, fingerprint
  )
  if (done) {
    message(sprintf(
      "[class %d/%d] %s | %s | %s  skip (done)",
      i, nrow(grid), g$level, g$config, g$method
    ))
    next
  }
  t0 <- Sys.time()
  s <- suppressWarnings(
    run_class_sweep_leaf(bundle, g$level, g$config, g$method, b_grid = b_grid)
  )
  row <- s[s$B == max(s$B), ]
  message(sprintf(
    "[class %d/%d] %s | %s | %s  AUC=%.2f p=%.3f  (%.0fs)",
    i, nrow(grid), g$level, g$config, g$method,
    row$estimate, row$perm_p,
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
}
message("classification run done for: ", paste(levels_run, collapse = ", "))
