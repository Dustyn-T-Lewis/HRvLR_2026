# F04 Panel B: module-level NES scatters. Each WGCNA module is a gene set; fgsea
# runs its members' moderated-t ranks per contrast. Top scatter reads training
# concordance (do HR and LR move the same modules the same way); bottom reads
# whether the acute response reinforces or reverses the chronic (training) one.
pacman::p_load(
  here, readr, dplyr, tidyr, tibble, stringr, purrr, ggplot2, patchwork, cowplot
)
source(here("04_Figures", "functions", "style.R"))
source(here("04_Figures", "functions", "pathway_utils.R"))

if (!exists("RPT_DIR")) RPT_DIR <- here("04_Figures", "F04", "b_reports")
if (!exists("DAT_DIR")) DAT_DIR <- here("04_Figures", "F04", "c_data")
PANEL_DIR <- file.path(RPT_DIR, "panels")
dir.create(PANEL_DIR, recursive = TRUE, showWarnings = FALSE)

is_light_color <- function(col) {
  rgb <- grDevices::col2rgb(col) / 255
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  lum > 0.6
}

mem <- read_csv(file.path(DAT_DIR, "wgcna_membership.csv"), show_col_types = FALSE)
combined <- read_csv(
  here("03_DEP", "a_non_imputed", "c_data", "03_combined_results.csv"),
  show_col_types = FALSE
)
ora <- read_csv(file.path(DAT_DIR, "module_ora.csv"), show_col_types = FALSE)

mem_real <- mem |> filter(group_id != "grey")
module_sets <- split(mem_real$protein_id, mem_real$group_id)
mod_sizes <- vapply(module_sets, length, integer(1))

# M1..Mk by size (largest first), matching Panel A's row order.
mod_order <- names(sort(mod_sizes, decreasing = TRUE))
mod_id <- setNames(sprintf("M%d", seq_along(mod_order)), mod_order)

bio_labels <- ora |>
  group_by(group) |>
  slice_min(padj, n = 1, with_ties = FALSE) |>
  ungroup() |>
  transmute(module_color = group, bio_label = clean_pathway_name(pathway))

build_ranks <- function(df, col) {
  vals <- df[[col]]
  names(vals) <- df$uniprot_id
  vals <- vals[!is.na(vals)]
  sort(vals, decreasing = TRUE)
}

run_module_fgsea <- function(ranks, module_sets) {
  set.seed(42)
  res <- fgsea::fgseaMultilevel(
    pathways = module_sets, stats = ranks,
    minSize = 15, maxSize = 500, nPermSimple = 10000, eps = 0
  )
  as.data.frame(res)
}
fgsea_tr_hr <- run_module_fgsea(build_ranks(combined, "t_Training_HR"), module_sets)
fgsea_tr_lr <- run_module_fgsea(build_ranks(combined, "t_Training_LR"), module_sets)
fgsea_ac_hr <- run_module_fgsea(build_ranks(combined, "t_Acute_HR"), module_sets)

merge_fgsea <- function(res, suffix) {
  res |>
    as_tibble() |>
    dplyr::select(pathway, NES, padj, size) |>
    dplyr::rename_with(~ paste0(., "_", suffix), c(NES, padj, size))
}
fgsea_wide <- merge_fgsea(fgsea_tr_hr, "Training_HR") |>
  left_join(merge_fgsea(fgsea_tr_lr, "Training_LR"), by = "pathway") |>
  left_join(merge_fgsea(fgsea_ac_hr, "Acute_HR"), by = "pathway") |>
  mutate(
    module_color = pathway,
    n_proteins = mod_sizes[pathway],
    module_id = mod_id[pathway]
  ) |>
  left_join(bio_labels, by = "module_color") |>
  mutate(
    bio_label = ifelse(is.na(bio_label), str_to_title(module_color), bio_label),
    bio_label = str_wrap(str_trunc(bio_label, 26, ellipsis = "..."), width = 14)
  )

# green reads visually light despite is_light_color returning FALSE
label_text_color <- function(cols) {
  out <- ifelse(vapply(cols, is_light_color, logical(1)), "black", "white")
  out[cols %in% c("green")] <- "black"
  out
}

build_scatter <- function(df, x_col, y_col, x_lab, y_lab, quad) {
  x_vals <- df[[x_col]]
  y_vals <- df[[y_col]]
  nes_lim <- max(abs(c(x_vals, y_vals)), na.rm = TRUE) * 1.35

  q_tr <- sum(x_vals > 0 & y_vals > 0, na.rm = TRUE)
  q_bl <- sum(x_vals < 0 & y_vals < 0, na.rm = TRUE)
  q_br <- sum(x_vals > 0 & y_vals < 0, na.rm = TRUE)
  q_tl <- sum(x_vals < 0 & y_vals > 0, na.rm = TRUE)

  df$label_col <- label_text_color(df$module_color)
  txt <- 5.0

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    annotate("rect",
      xmin = 0, xmax = Inf, ymin = 0, ymax = Inf,
      fill = quad$fill[1], alpha = 0.20
    ) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0,
      fill = quad$fill[2], alpha = 0.20
    ) +
    annotate("rect",
      xmin = 0, xmax = Inf, ymin = -Inf, ymax = 0,
      fill = quad$fill[3], alpha = 0.20
    ) +
    annotate("rect",
      xmin = -Inf, xmax = 0, ymin = 0, ymax = Inf,
      fill = quad$fill[4], alpha = 0.20
    ) +
    geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed",
      color = "black", linewidth = 0.3
    ) +
    geom_point(aes(size = n_proteins),
      fill = df$module_color, color = "black",
      shape = 21, alpha = 0.85, stroke = 0.5
    ) +
    geom_text(aes(label = module_id),
      color = df$label_col,
      size = txt * 0.65, fontface = "bold", show.legend = FALSE
    ) +
    annotate("label",
      x = Inf, y = Inf,
      label = sprintf("%s  n=%d", quad$label[1], q_tr),
      hjust = 1, vjust = 1, size = txt, fontface = "bold",
      color = quad$text_col[1], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = -Inf, y = -Inf,
      label = sprintf("%s  n=%d", quad$label[2], q_bl),
      hjust = 0, vjust = 0, size = txt, fontface = "bold",
      color = quad$text_col[2], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = Inf, y = -Inf,
      label = sprintf("%s  n=%d", quad$label[3], q_br),
      hjust = 1, vjust = 0, size = txt, fontface = "bold",
      color = quad$text_col[3], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    annotate("label",
      x = -Inf, y = Inf,
      label = sprintf("%s  n=%d", quad$label[4], q_tl),
      hjust = 0, vjust = 1, size = txt, fontface = "bold",
      color = quad$text_col[4], fill = scales::alpha("white", 0.92),
      label.padding = unit(2, "pt")
    ) +
    scale_size_continuous(
      range = c(3, 9), name = "Proteins",
      breaks = c(50, 100, 200, 300)
    ) +
    scale_x_continuous(expand = expansion(mult = 0.02)) +
    scale_y_continuous(expand = expansion(mult = 0.02)) +
    coord_fixed(
      ratio = 1, xlim = c(-nes_lim, nes_lim),
      ylim = c(-nes_lim, nes_lim), clip = "off"
    ) +
    labs(x = x_lab, y = y_lab) +
    FIG_THEME +
    theme(
      axis.text = element_text(size = 11, face = "bold", color = "grey30"),
      axis.title.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = -4)),
      legend.position = "none",
      plot.margin = margin(2, 4, 2, -8)
    )
}

warm <- unname(DIR_COLORS["Up"])
cool <- unname(DIR_COLORS["Down"])
green <- "#2E7D32"

quad_conc <- list(
  label = c("Concordant Up", "Concordant Down", "Discordant", "Discordant"),
  fill = c(warm, warm, cool, cool),
  text_col = unname(GROUP_COLORS[c("LR", "LR", "HR", "HR")])
)
p_top <- build_scatter(
  fgsea_wide, "NES_Training_LR", "NES_Training_HR",
  "NES (Training LR)", "NES (Training HR)", quad_conc
)

quad_rev <- list(
  label = c("Reinforced", "Reinforced", "Reversed", "Reversed"),
  fill = c(warm, warm, green, green),
  text_col = c(unname(GROUP_COLORS["LR"]), unname(GROUP_COLORS["LR"]), green, green)
)
p_bottom <- build_scatter(
  fgsea_wide, "NES_Acute_HR", "NES_Training_HR",
  "NES (Acute HR)", "NES (Training HR)", quad_rev
)

mod_key_df <- fgsea_wide |>
  dplyr::select(module_color, module_id, bio_label) |>
  arrange(match(module_id, mod_id))

build_key_cell <- function(mod, mid, lab) {
  ggplot() +
    annotate("rect",
      xmin = 0.02, xmax = 0.10, ymin = 0.25, ymax = 0.75,
      fill = mod, color = "grey20", linewidth = 0.4
    ) +
    annotate("text",
      x = 0.06, y = 0.50, label = mid,
      color = ifelse(is_light_color(mod) | mod %in% c("green"), "black", "white"),
      fontface = "bold", size = 3.2
    ) +
    annotate("text",
      x = 0.13, y = 0.50, label = lab,
      hjust = 0, vjust = 0.5, size = 3.0, color = "grey10"
    ) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(plot.margin = margin(4, 1, 4, 1))
}
key_cells <- purrr::pmap(
  list(mod_key_df$module_color, mod_key_df$module_id, mod_key_df$bio_label),
  build_key_cell
)
key_ncol <- 3L
key_nrow <- ceiling(length(key_cells) / key_ncol)
module_key_strip <- wrap_plots(key_cells, nrow = key_nrow, ncol = key_ncol)

scatters_panel <- (p_top / p_bottom / module_key_strip) +
  plot_layout(heights = c(1, 1, 0.42))

PB_W <- 220
PB_H <- 320
save_png(scatters_panel, file.path(PANEL_DIR, "panel_b_nes_scatters"), PB_W, PB_H)
ggsave(file.path(PANEL_DIR, "panel_b_nes_scatters.pdf"), scatters_panel,
  width = PB_W, height = PB_H, units = "mm", device = PDF_DEVICE
)

p_legend_src <- p_top +
  scale_size_continuous(
    range = c(3, 9), name = "Proteins",
    breaks = c(50, 100, 200, 300),
    guide = guide_legend(
      nrow = 1,
      override.aes = list(alpha = 0.7, fill = "grey60", stroke = 0)
    )
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_rect(fill = NA, color = NA),
    legend.key.size = unit(5, "mm"),
    legend.background = element_rect(fill = NA, color = NA)
  )
legend_grob <- cowplot::get_plot_component(p_legend_src, "guide-box-bottom",
  return_all = FALSE
)
p_legend <- cowplot::ggdraw(legend_grob)
save_png(p_legend, file.path(PANEL_DIR, "panel_b_nes_legend"), 90, 16)

write_csv(fgsea_wide, file.path(DAT_DIR, "panel_b_module_fgsea.csv"))
message("F04 Panel B (NES scatters) done.")
