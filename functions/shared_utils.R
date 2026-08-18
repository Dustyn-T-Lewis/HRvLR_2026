# Shared pipeline utilities, sourced across stages and figures (alongside pca.R).

pacman::p_load(digest, openxlsx)

# Empty a directory of everything but its .gitkeep, creating it if absent. Used
# before a stage or figure writes, so a rerun never leaves stale outputs behind.
clear_dir <- function(d) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  unlink(setdiff(list.files(d, full.names = TRUE), file.path(d, ".gitkeep")), recursive = TRUE)
}

# Digest of whatever inputs a stage was built from, stamped into its workbook so
# a rerun can tell whether the output is stale.
input_fingerprint <- function(...) {
  substr(digest::digest(list(...)), 1, 12)
}

write_workbook <- function(path, sheets, fingerprint = NULL) {
  if (!is.null(fingerprint)) {
    sheets$provenance <- data.frame(fingerprint = fingerprint)
  }
  wb <- createWorkbook()
  for (nm in names(sheets)) {
    addWorksheet(wb, nm)
    writeData(wb, nm, sheets[[nm]])
  }
  saveWorkbook(wb, path, overwrite = TRUE)
}
