# What every composite shares: run the figure's panel scripts, tag a panel, write the workbook.
# Each panel script leaves panel (the plot), panel_data (the numbers drawn), holds (one line on what
# panel_data holds) and inputs (the files it read, relative to the project root). Expects
# figure_style.R sourced.
source(here::here("00_Setup", "Shared_functions", "workbooks.R"))

# Every panel script of a figure, in letter order, each run in its own environment.
run_panels <- function(figure) {
  files <- list.files(
    here::here("04_Figures", figure, "a_script", "panels"), "\\.R$",
    full.names = TRUE
  )
  envs <- lapply(files, \(file) {
    env <- new.env()
    sys.source(file, envir = env)
    env
  })
  rlang::set_names(envs, sub("\\.R$", "", basename(files)))
}

# A panel that is itself a patchwork carries its letter on its first plot.
tag_panel <- function(plot, letter) {
  if (inherits(plot, "patchwork")) {
    plot[[1]] <- plot[[1]] + ggplot2::labs(tag = letter) + figure_tags
    return(plot)
  }
  plot + ggplot2::labs(tag = letter) + figure_tags
}

# One sheet per panel, named by its script, then the files read and the packages loaded.
write_figure_workbook <- function(envs, figure) {
  inputs <- unlist(unname(lapply(envs, `[[`, "inputs")))
  inputs <- inputs[!duplicated(inputs)]
  write_workbook(
    lapply(envs, `[[`, "panel_data"), unname(vapply(envs, `[[`, "", "holds")),
    here::here("04_Figures", figure, "c_data", paste0(figure, ".xlsx")), inputs,
    here::here(inputs)
  )
}

# A figure's composite PDF at 04_Figures/<figure>/b_reports/<figure>_composite.pdf, landscape by
# default (250 mm wide).
save_figure <- function(plot, figure, height, width = 250) {
  ggplot2::ggsave(
    here::here("04_Figures", figure, "b_reports", paste0(figure, "_composite.pdf")), plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf
  )
}
