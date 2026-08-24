# HRvLR 2026 - restore the package library and record what was installed.
# Run once after cloning: Rscript setup.R
#
# renv.lock is the authority: 392 packages, 43 from Bioconductor and 3 from
# GitHub (proteoDA, enrichVolcano, RRHO2), pinned to the versions the committed
# figures were rendered with. Nothing is listed here a second time, because a
# hand-kept list drifts from the lockfile the moment either changes.
#
# ragg matters more than its absence from any script suggests. ggplot2::ggsave()
# picks ragg::agg_png() when ragg is installed and grDevices::png() when it is
# not, and the two shape text differently, so a library without it re-renders
# every committed figure to different bytes with no error. _dependencies.R
# exists to keep ragg, lintr and styler visible to renv; do not delete it.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

renv::restore(prompt = FALSE)

lock <- jsonlite::fromJSON("renv.lock")
pkgs <- sort(names(lock$Packages))

# Read DESCRIPTION rather than load the package: rgl needs OpenGL, which this
# machine lacks, and requireNamespace() would report an installed package as
# missing purely because its shared library will not open.
describe <- function(pkg) {
  desc <- utils::packageDescription(pkg)
  if (!is.list(desc)) {
    return("not installed")
  }
  if (is.null(desc$RemoteSha)) {
    desc$Version
  } else {
    paste0(desc$Version, " (", substr(desc$RemoteSha, 1, 7), ")")
  }
}

writeLines(
  c(
    R.version.string,
    paste("Platform:", R.version$platform),
    paste("Recorded:", format(Sys.Date())),
    paste("Lockfile:", length(pkgs), "packages"),
    "",
    sprintf("%-28s %s", pkgs, vapply(pkgs, describe, character(1)))
  ),
  "package_versions.txt"
)

cat(sprintf("restored %d packages; wrote package_versions.txt\n", length(pkgs)))
