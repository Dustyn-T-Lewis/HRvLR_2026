# What every composite shares: run the figure's panel scripts, tag a panel, save the figure and its
# workbook. Each panel script (a_*.R, b_*.R, ... beside the composite; in supp/, prefixed with the
# supplement's code, S4_a_*.R) leaves panel (the plot) and panel_data (the numbers drawn). Expects
# figure_style.R sourced.

# Every panel script of a figure, in letter order, each run in its own environment.
run_panels <- function(figure) {
  pattern <- paste0("^", panel_prefix(figure), "[a-z]_.+\\.R$")
  files <- list.files(figure_dir(figure, "a_script"), pattern, full.names = TRUE)
  envs <- lapply(files, \(file) {
    env <- new.env()
    sys.source(file, envir = env)
    env
  })
  names <- sub(paste0("^", panel_prefix(figure), "(.+)\\.R$"), "\\1", basename(files))
  rlang::set_names(envs, names)
}

# A panel that is itself a patchwork carries its letter on its first plot.
tag_panel <- function(plot, letter) {
  if (inherits(plot, "patchwork")) {
    plot[[1]] <- plot[[1]] + ggplot2::labs(tag = letter) + figure_tags
    return(plot)
  }
  plot + ggplot2::labs(tag = letter) + figure_tags
}

# The figure's workbook: one sheet per panel, named by its script, holding the numbers drawn.
write_figure_workbook <- function(envs, figure) {
  writexl::write_xlsx(
    lapply(envs, `[[`, "panel_data"),
    file.path(figure_dir(figure, "c_data"), paste0(figure, ".xlsx"))
  )
}

# The composite PDF beside the panels, landscape by default (250 mm wide).
save_figure <- function(plot, figure, height, width = 250) {
  ggplot2::ggsave(
    file.path(figure_dir(figure, "b_reports"), paste0(figure, "_composite.pdf")), plot,
    width = width, height = height, units = "mm", device = grDevices::cairo_pdf, create.dir = TRUE
  )
}
