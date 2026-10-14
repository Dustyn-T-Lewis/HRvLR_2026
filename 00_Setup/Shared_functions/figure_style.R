# Colours and theme shared by every figure in every folder. The contrast palette is sampled from the
# design panel (05_Figures/HRvLR/F1_Overview/assets/design.png), so a contrast family looks
# the same in the schematic and in every plot.
#
# Each family has three tones: line is the panel's own colour, for marks; fill is the panel's
# box tint, behind strip labels; ink is a darker shade of the same hue for text, because the
# panel's orange and yellow are too light to read at 6 pt on white.

contrast_palette <- tibble::tribble(
  ~family, ~line, ~fill, ~ink,
  "Training", "#317EC2", "#B9DBF4", "#1F5A8F",
  "Acute", "#E7872B", "#F6BF93", "#A85A0E",
  "HR - LR", "#4F609C", "#E5EBFB", "#36447A",
  "Hi - Lo", "#4F609C", "#E5EBFB", "#36447A",
  "C1 - C2", "#4F609C", "#E5EBFB", "#36447A",
  "Between groups", "#4F609C", "#E5EBFB", "#36447A",
  "Interaction", "#E9C61D", "#F4E6AA", "#8A7208"
)
contrast_tone <- \(tone) tibble::deframe(contrast_palette[c("family", tone)])

# Within Training and Acute, the first group (HR) takes the ink and the second (LR) the line colour,
# so the first reads darker, as in the design panel, and both stay legible.
arm_shade <- function(family, arm, groups = c("HR", "LR")) {
  purrr::map2_chr(family, arm, \(f, a) {
    row <- contrast_palette[contrast_palette$family == f, ]
    if (!nrow(row) || !a %in% groups) {
      return("grey30")
    }
    if (a == groups[1]) row$ink else row$line
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
# colour, with the name in bold black. A level outside the palette (the database column in Figure 1)
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

# The two thresholds the figures use, and one way of printing p and FDR.
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

# Contrast names from a design, for any pair of groups (HR/LR, Hi/Lo, C1/C2): a readable label, the
# family that sets the colour, and the window that groups the rows. A between-group suffix reads
# "HRvLR"; a within-group suffix is one group's name. The first group takes the darker shade.
between <- "[A-Za-z0-9]+v[A-Za-z0-9]+$"
contrast_label <- function(contrast) {
  contrast |>
    sub(pattern = "([A-Za-z0-9]+)v([A-Za-z0-9]+)$", replacement = "\\1 - \\2") |>
    gsub(pattern = "_", replacement = " ")
}
contrast_family <- function(contrast) {
  dplyr::case_when(
    grepl("Interaction", contrast) ~ "Interaction",
    grepl(between, contrast) ~ "HR - LR",
    startsWith(contrast, "Training") ~ "Training",
    startsWith(contrast, "Acute") ~ "Acute"
  )
}
contrast_window <- function(contrast) {
  dplyr::case_when(
    grepl(between, contrast) ~ "Between groups",
    startsWith(contrast, "Training") ~ "Training",
    startsWith(contrast, "Acute") ~ "Acute"
  )
}

# The strip a contrast sits under from F2 on, which names no group, and its short column label:
# the group for a within-group window, the time point for a between-group level, the window for an
# interaction.
panel_families <- c("Training", "Acute", "Between groups", "Interaction")
contrast_panel <- function(contrast) {
  family <- dplyr::if_else(grepl("Interaction", contrast), "Interaction", contrast_window(contrast))
  factor(family, panel_families)
}
contrast_short <- function(contrast) {
  prefix <- sub("_.*", "", contrast)
  dplyr::case_when(
    grepl("Interaction", contrast) ~ prefix,
    grepl(between, contrast) ~ c(Baseline = "T1", Trained = "T2", Acute = "T3")[prefix],
    .default = sub(".*_", "", contrast)
  )
}
# The two groups of any grouping in the design panel's HR (dark) and LR (light) greys.
group_colours <- \(groups) rlang::set_names(unname(arm_colours), groups)

# Trait column names as axis text: "pct_leg_press_1rm" reads "leg press 1RM (% change)".
trait_label <- function(trait) {
  base <- trait |>
    sub(pattern = "^pct_", replacement = "") |>
    gsub(pattern = "_", replacement = " ") |>
    sub(pattern = "1rm", replacement = "1RM") |>
    sub(pattern = "fcsa", replacement = "fCSA") |>
    sub(pattern = "mcsa", replacement = "mCSA")
  ifelse(startsWith(trait, "pct_"), paste(base, "(% change)"), base)
}
