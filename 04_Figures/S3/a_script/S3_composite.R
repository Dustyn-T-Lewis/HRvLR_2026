# S3 runs far past one page, so the composite is one PDF with a page per panel (panel a: a page per
# window), each 250 mm wide and as tall as its panel; e is the Hallmark score heatmap.
library(here)
library(ggplot2)
source(here("04_Figures", "style", "figure_style.R"))
source(here("04_Figures", "style", "composite.R"))

panels <- run_panels("S3")
stopifnot(names(panels) == c(
  "a_tree_matrix", "b_cut_075", "c_cut_09", "d_go_slim", "e_hallmark_heatmap"
))

page <- \(plot, tag, height, name) {
  file <- file.path(tempdir(), paste0(name, ".pdf"))
  ggsave(
    file, plot + labs(tag = tag) + figure_tags,
    width = 250, height = height, units = "mm", device = grDevices::cairo_pdf
  )
  file
}
tree <- panels$a_tree_matrix
pages <- c(
  purrr::imap_chr(tree$panel, \(plot, w) page(plot, "a", tree$panel_height[[w]], paste0("a_", w))),
  purrr::imap_chr(panels[-1], \(env, name) {
    page(env$panel, substr(name, 1, 1), env$panel_height, name)
  })
)
invisible(qpdf::pdf_combine(pages, here("04_Figures", "S3", "b_reports", "S3_composite.pdf")))

write_figure_workbook(panels, "S3")
