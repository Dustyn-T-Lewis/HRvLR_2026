# Colours and theme shared by every figure in 06_Figures. The contrast palette is sampled from the
# design panel (06_Figures/F01_Overview/assets/design.png), so a contrast family looks the same in
# the schematic and in every plot.
#
# Each family has three tones: `line` is the panel's own colour, for marks; `fill` is the panel's
# box tint, behind strip labels; `ink` is a darker shade of the same hue for text, because the
# panel's orange and yellow are too light to read at 6 pt on white.

contrast_palette <- tibble::tribble(
  ~family, ~line, ~fill, ~ink,
  "Training", "#317EC2", "#B9DBF4", "#1F5A8F",
  "Acute", "#E7872B", "#F6BF93", "#A85A0E",
  "HR - LR", "#4F609C", "#E5EBFB", "#36447A",
  "Interaction", "#E9C61D", "#F4E6AA", "#8A7208"
)
contrast_tone <- \(tone) tibble::deframe(contrast_palette[c("family", tone)])

# Within Training and Acute, HR takes the ink and LR the line colour, so HR reads darker, as in the
# design panel, and both stay legible.
arm_shade <- function(family, arm) {
  purrr::map2_chr(family, arm, \(f, a) {
    row <- contrast_palette[contrast_palette$family == f, ]
    if (!nrow(row) || !a %in% c("HR", "LR")) {
      return("grey30")
    }
    if (a == "HR") row$ink else row$line
  })
}

# Time points take the colour of the window that produced them: T1 neutral, T2 Training, T3 Acute.
time_colours <- c(
  T1 = "#8C8C8C",
  T2 = contrast_palette$line[contrast_palette$family == "Training"],
  T3 = contrast_palette$line[contrast_palette$family == "Acute"]
)
arm_shapes <- c(HR = 16, LR = 17)
direction_colours <- c(Up = "#4D9221", Down = "#C51B7D")

figure_theme <- ggplot2::theme_minimal(base_size = 7) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 7.5, margin = ggplot2::margin(b = 1)),
    plot.subtitle = ggplot2::element_text(
      size = 6, colour = "grey30", margin = ggplot2::margin(b = 3)
    ),
    plot.title.position = "plot",
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major = ggplot2::element_line(linewidth = 0.2),
    legend.key.size = ggplot2::unit(3, "mm")
  )

# Panel letters sit in a left margin column, level with the title.
figure_tags <- ggplot2::theme(
  plot.tag = ggplot2::element_text(
    face = "bold", size = 9, vjust = 1, margin = ggplot2::margin(r = 2)
  ),
  plot.tag.location = "margin", plot.tag.position = "left"
)

# Facet strips as the design panel's labelled boxes: family tint inside a border in the family
# colour, with the name in bold black. A level outside the palette (the database column in F01)
# gets no box.
contrast_strips <- function(x = NULL, y = NULL) {
  line <- contrast_tone("line")
  fill <- contrast_tone("fill")
  text <- \(levels) ggh4x::elem_list_text(colour = "black", face = "bold")
  box <- \(levels) {
    ggh4x::elem_list_rect(
      fill = unname(dplyr::coalesce(fill[levels], "transparent")),
      colour = unname(dplyr::coalesce(line[levels], "transparent")),
      linewidth = 0.6
    )
  }
  ggh4x::strip_themed(
    text_x = if (!is.null(x)) text(x), background_x = if (!is.null(x)) box(x),
    text_y = if (!is.null(y)) text(y), background_y = if (!is.null(y)) box(y)
  )
}

# Axis labels coloured one by one, rendered with ggtext::element_markdown().
in_colour <- \(text, colour) sprintf("<span style='color:%s'>%s</span>", colour, text)
