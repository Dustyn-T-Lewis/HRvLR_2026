# F01 phenotype atlas — render the standalone panels (a-i) into b_reports.
# No composite; panels are assembled into the figure manually.
setwd(rprojroot::find_rstudio_root_file())

RPT <- "04_Figures/F01/b_reports"
unlink(setdiff(list.files(RPT, full.names = TRUE), file.path(RPT, ".gitkeep")),
  recursive = TRUE
)
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)

for (p in LETTERS[1:9]) {
  source(sprintf("04_Figures/F01/a_script/panel_%s.R", p))
}

source("04_Figures/F01/a_script/HRvLR_F01_composite.R")

cat("F01 phenotype panels + composite rendered to", RPT, "\n")
