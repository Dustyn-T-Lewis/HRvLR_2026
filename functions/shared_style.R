# Palettes, theme and the packet writer every stage's last step uses.

pacman::p_load(here, ggplot2, scales, stringr, qpdf)

source(here("functions", "shared_utils.R"))

# Two blue/red mappings coexist and must not be conflated: GROUP_COLORS encode
# responder (HR dark blue, LR dark red); DIR_COLORS encode direction (up light
# red, down light blue). The group hues are the darker shades so a responder
# legend never reads as a direction legend.
GROUP_COLORS <- c(HR = "#2166AC", LR = "#B2182B")
DIR_COLORS <- c(Up = "#D6604D", Down = "#4393C3", NS = "grey70")

# Okabe-Ito, colourblind-safe: baseline, trained, acute.
TIME_COLORS <- c(T1 = "#E69F00", T2 = "#0072B2", T3 = "#009E73")

# Responder family sets the hue (HR blue, LR red, between-arm green,
# interaction purple); timepoint sets the shade.
CONTRAST_COLORS <- c(
  Training_HR          = "#6BAED6",
  Acute_HR             = "#2166AC",
  Training_LR          = "#FC9272",
  Acute_LR             = "#B2182B",
  Baseline_HRvLR       = "#A1D99B",
  Trained_HRvLR        = "#41AB5D",
  Acute_HRvLR          = "#238B45",
  Training_Interaction = "#9E9AC8",
  Acute_Interaction    = "#6A51A3"
)

DB_COLORS <- c(
  Hallmark = "#AA336A", Reactome = "#1565C0", "GO:BP" = "#00796B",
  "GO Slim" = "#5D4037"
)

FIG_THEME <- theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, color = "grey30"),
    plot.caption = element_text(hjust = 0, size = 8, color = "grey25"),
    plot.caption.position = "plot",
    plot.title.position = "plot",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )

# Captions are long by design; wrap them to the page rather than by hand.
caption <- function(...) str_wrap(paste0(...), width = 150)

# MSigDB names are shouting snake case with a collection prefix; a label needs
# neither.
clean_set_name <- function(set, width = 50) {
  set |>
    str_remove("^(HALLMARK|REACTOME|GOBP|GOSLIM)_") |>
    str_replace_all("_", " ") |>
    str_to_sentence() |>
    str_trunc(width)
}

# A packet is a contents page followed by one page per plot, merged into one
# PDF. Each page carries its own title, subtitle (method and counts) and
# caption (what each channel encodes, which table holds the data), so a page
# pulled out of the packet still explains itself.
write_packet <- function(pages, path, title, width = 280, height = 200) {
  page_dir <- file.path(dirname(path), "pages")
  clear_dir(page_dir)
  contents <- ggplot() +
    annotate(
      "text",
      x = 0, y = -seq_along(pages), hjust = 0, size = 4,
      label = sprintf("%2d   %s", seq_along(pages) + 1L, names(pages))
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(-length(pages) - 1, 0)) +
    labs(title = title, subtitle = paste("Generated", Sys.Date())) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.margin = margin(20, 20, 20, 20)
    )
  pages <- c(list(contents), pages)
  files <- file.path(page_dir, sprintf("p%02d.pdf", seq_along(pages)))
  for (i in seq_along(pages)) {
    ggsave(files[i], pages[[i]],
      width = width, height = height, units = "mm", device = "pdf",
      bg = "white"
    )
  }
  qpdf::pdf_combine(files, path)
  unlink(page_dir, recursive = TRUE)
  invisible(path)
}
