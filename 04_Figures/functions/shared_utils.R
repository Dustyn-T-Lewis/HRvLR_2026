# Shared pipeline utilities, sourced across stages and figures (alongside pca.R).

# Empty a directory of everything but its .gitkeep, creating it if absent. Used
# before a stage or figure writes, so a rerun never leaves stale outputs behind.
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}
