# Draws pipeline.png at the repository root: every step in run_all.R, one row per folder, and the
# files each step read, taken from the input_manifest sheet of its workbooks. Run it after
# run_all.R so the manifests are current.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(ggplot2)
})

run_all <- parse(here("run_all.R"))
stages_at <- which(map_lgl(as.list(run_all), \(e) identical(e[[2]], quote(stages))))
stages <- eval(run_all[[stages_at]][[3]])

# The three cluster groupings run the same steps, so the row shows Muscle_size's, plus the
# one step only Proteome has.
rows <- tibble(stage = stages) |>
  mutate(
    folder = sub("/.*", "", stage),
    shown = !grepl("^03_Cluster_groups/(Muscle_strength|Proteome)/", stage) |
      grepl("08_Cluster_Phenotype", stage),
    step = sub("^(Muscle_size|Proteome)/", "", sub("^[^/]+/", "", stage))
  )

manifest <- map(set_names(stages), \(stage) {
  files <- list.files(here(stage), "[.]xlsx$", recursive = TRUE, full.names = TRUE)
  files <- files[!grepl("_all_inputs[.]xlsx$", files)]
  map(files, \(f) readxl::read_excel(f, "input_manifest")$path) |>
    unlist() |>
    unique()
})
producer <- function(path) {
  if (startsWith(path, "00_Input")) {
    return("00_Input")
  }
  hit <- stages[startsWith(path, paste0(stages, "/"))]
  stopifnot(length(hit) > 0)
  hit[which.max(nchar(hit))]
}
edges <- imap(manifest, \(paths, stage) tibble(from = map_chr(paths, producer), to = stage)) |>
  list_rbind() |>
  distinct() |>
  filter(from != to)

lanes <- c(
  "01_Processing" = "01_Processing: shared steps",
  "02_HRvLR_groups" = "02_HRvLR_groups: the study's HR and LR",
  "03_Cluster_groups" = "03_Cluster_groups: Muscle_size, Muscle_strength, Proteome (same steps)",
  "04_Continuous_no_group" = "04_Continuous_no_group: traits as numbers, no groups"
)
nodes <- rows |>
  filter(shown) |>
  mutate(
    lane = match(folder, names(lanes)),
    y = length(lanes) - lane + 1,
    x = row_number() + 1, .by = folder
  ) |>
  mutate(
    figure = grepl("Figures/", stage),
    label = gsub("[_/]", " ", sub("^Figures/", "", step)) |> stringr::str_wrap(13)
  )
# Steps of the other two groupings point at their Muscle_size twin.
twin <- function(stage) {
  sub("^03_Cluster_groups/(Muscle_strength|Proteome)/", "03_Cluster_groups/Muscle_size/", stage) |>
    (\(s) if_else(s %in% nodes$stage, s, stage))()
}
placed <- edges |>
  mutate(from = if_else(from == "00_Input", from, twin(from)), to = twin(to)) |>
  distinct() |>
  filter(from == "00_Input" | from %in% nodes$stage, to %in% nodes$stage) |>
  left_join(select(nodes, from = stage, from_folder = folder, x0 = x, y0 = y), by = "from") |>
  left_join(select(nodes, to = stage, to_folder = folder, x1 = x, y1 = y), by = "to")

within <- filter(placed, from_folder == to_folder)
across <- placed |>
  filter(from == "00_Input" | from_folder != to_folder) |>
  mutate(source = if_else(from == "00_Input", "00_Input", from_folder)) |>
  distinct(source, to_folder, y1) |>
  filter(source != to_folder)
input_rows <- filter(across, source == "00_Input")
processing_rows <- filter(across, source == "01_Processing")
between_folders <- placed |>
  filter(from != "00_Input", from_folder != to_folder, from_folder != "01_Processing")

box_w <- 0.86
box_h <- 0.42
# Each arc rises the same small height whatever its span, so arcs stay between the rows.
arc_layers <- within |>
  mutate(span = x1 - x0) |>
  split(~span) |>
  map(\(e) {
    geom_curve(
      data = e, aes(x = x0, y = y0 + box_h / 2, xend = x1, yend = y1 + box_h / 2),
      curvature = -(0.5 + 0.03 * e$span[1]) / e$span[1], ncp = 10,
      arrow = arrow(length = unit(1.2, "mm"), type = "closed"),
      linewidth = 0.25, colour = "#2166AC", alpha = 0.5
    )
  })
plot <- ggplot() +
  annotate("rect",
    xmin = 0.45, xmax = 1.2, ymin = 0.55, ymax = length(lanes) + 0.3,
    fill = "#F2F2F2", colour = "grey40", linewidth = 0.3
  ) +
  annotate("text",
    x = 0.825, y = (length(lanes) + 0.85) / 2, label = "00_Input\nstudy data\nand reference\nfiles",
    size = 2.6, lineheight = 0.95
  ) +
  geom_segment(
    data = input_rows, aes(x = 1.2, xend = 2 - box_w / 2, y = y1, yend = y1),
    arrow = arrow(length = unit(1.4, "mm"), type = "closed"), linewidth = 0.3, colour = "grey35"
  ) +
  geom_segment(
    data = processing_rows,
    aes(x = 1.45, xend = 1.45, y = length(lanes) - box_h / 2, yend = y1 + 0.12),
    linewidth = 0.4, colour = "#B2182B"
  ) +
  geom_segment(
    data = processing_rows, aes(x = 1.45, xend = 2 - box_w / 2, y = y1 + 0.12, yend = y1 + 0.12),
    arrow = arrow(length = unit(1.4, "mm"), type = "closed"), linewidth = 0.4, colour = "#B2182B"
  ) +
  arc_layers +
  geom_curve(
    data = between_folders, aes(x = x0, y = y0 + box_h / 2, xend = x1, yend = y1 - box_h / 2),
    curvature = 0.15, linetype = "dashed", arrow = arrow(length = unit(1.4, "mm"), type = "closed"),
    linewidth = 0.35, colour = "grey25"
  ) +
  geom_tile(
    data = nodes, aes(x, y, fill = figure), width = box_w, height = box_h,
    colour = "grey30", linewidth = 0.25
  ) +
  geom_text(data = nodes, aes(x, y, label = label), size = 2, lineheight = 0.9) +
  annotate("text",
    x = 2 - box_w / 2, y = rev(seq_along(lanes)) - box_h / 2 - 0.1, label = lanes, hjust = 0,
    vjust = 1, size = 2.8, fontface = "bold"
  ) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#FDF1D8"), guide = "none") +
  scale_x_continuous(expand = expansion(add = 0.1)) +
  scale_y_continuous(limits = c(0.4, length(lanes) + 0.6), expand = expansion(0)) +
  labs(
    title = "HRvLR proteomics pipeline",
    subtitle = paste(
      "Blue: a step reads an earlier step's output. Red: the analysis folders read 01_Processing.",
      "Dashed: a step reads another folder. Yellow: figures. Edges come from each workbook's",
      "input_manifest."
    )
  ) +
  theme_void(base_size = 9) +
  theme(
    plot.background = element_rect(fill = "white", colour = NA),
    plot.title = element_text(face = "bold", margin = margin(b = 2)),
    plot.subtitle = element_text(size = 7.5, colour = "grey30", margin = margin(b = 6)),
    plot.margin = margin(8, 8, 8, 8)
  )
ggsave(here("pipeline.png"), plot, width = 13, height = 6.5, dpi = 200)
