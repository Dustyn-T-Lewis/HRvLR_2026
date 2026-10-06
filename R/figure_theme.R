# Colours and theme shared by every figure in every arc. The contrast palette is sampled from the
# design panel (02_HRvLR_groups/Figures/F01_Overview/assets/design.png), so a contrast family looks
# the same in the schematic and in every plot.
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
# HR dark and LR light grey, as the two figures in the design panel.
arm_colours <- c(HR = "#4D4D4D", LR = "#A6A6A6")
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

# The two thresholds every figure from F02 on uses, and one way of printing p and FDR.
alpha_nominal <- 0.05
alpha_fdr <- 0.05
label_p <- scales::label_pvalue(accuracy = 0.001, prefix = c("p < ", "p = ", "p > "))
label_fdr <- scales::label_pvalue(accuracy = 0.001, prefix = c("FDR < ", "FDR = ", "FDR > "))

# One evidence code for every marker: filled at FDR, open when only nominal, grey otherwise.
evidence_levels <- c(
  sprintf("FDR < %s", alpha_fdr), sprintf("p < %s only", alpha_nominal), "not significant"
)
evidence_level <- function(p, fdr) {
  level <- dplyr::case_when(
    fdr < alpha_fdr ~ evidence_levels[1],
    p < alpha_nominal ~ evidence_levels[2],
    .default = evidence_levels[3]
  )
  factor(level, evidence_levels)
}
evidence_shapes <- rlang::set_names(c(16, 1, 16), evidence_levels)
evidence_alpha <- rlang::set_names(c(1, 1, 0.35), evidence_levels)

# Contrast names from the design, as readable labels and with the family that sets their colour.
contrast_label <- function(contrast) {
  contrast |>
    sub(pattern = "HRvLR", replacement = "HR - LR") |>
    gsub(pattern = "_", replacement = " ")
}
contrast_family <- function(contrast) {
  dplyr::case_when(
    grepl("Interaction", contrast) ~ "Interaction",
    grepl("HRvLR", contrast) ~ "HR - LR",
    startsWith(contrast, "Training") ~ "Training",
    startsWith(contrast, "Acute") ~ "Acute"
  )
}
contrast_colour <- function(contrast) {
  arm <- dplyr::case_when(
    endsWith(contrast, "_HR") ~ "HR", endsWith(contrast, "_LR") ~ "LR",
    .default = NA_character_
  )
  family <- contrast_family(contrast)
  ifelse(is.na(arm), contrast_tone("line")[family], arm_shade(family, arm))
}

# The YvO figure theme (theme_bw, Helvetica, 6 pt base), used from F02 on.
paper_theme <- ggplot2::theme_bw(base_size = 6, base_family = "Helvetica") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 7, margin = ggplot2::margin(b = 1)),
    plot.subtitle = ggplot2::element_text(
      face = "bold.italic", size = 6, colour = "grey30", margin = ggplot2::margin(b = 2)
    ),
    plot.title.position = "plot",
    plot.tag = ggplot2::element_text(face = "bold", size = 8),
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold", size = 5),
    axis.title = ggplot2::element_text(face = "bold", size = 5),
    axis.text = ggplot2::element_text(size = 5, colour = "grey15"),
    legend.title = ggplot2::element_text(face = "bold", size = 5, colour = "grey20"),
    legend.text = ggplot2::element_text(face = "bold", size = 6, colour = "grey15"),
    legend.key.size = ggplot2::unit(2.5, "mm"),
    panel.grid.minor = ggplot2::element_blank()
  )
