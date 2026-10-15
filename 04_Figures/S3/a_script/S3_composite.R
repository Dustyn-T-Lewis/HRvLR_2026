# The panels together run past 240 mm, so the composite is one PDF with a page per panel, each page
# as tall as its panel.
library(here)
library(dplyr)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S3")
stopifnot(names(panels) == c("a_tree_matrix", "b_cut_075", "c_cut_09", "d_go_slim"))

pages <- purrr::imap_chr(panels, \(env, name) {
  file <- file.path(tempdir(), paste0(name, ".pdf"))
  ggsave(
    file, env$panel + labs(tag = substr(name, 1, 1)) + figure_tags,
    width = 180, height = env$panel_height, units = "mm", device = cairo_pdf
  )
  file
})
invisible(qpdf::pdf_combine(pages, here("04_Figures", "S3", "b_reports", "S3_composite.pdf")))

write_figure_workbook(panels, "S3")
