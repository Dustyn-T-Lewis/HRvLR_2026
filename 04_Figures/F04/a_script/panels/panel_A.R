# F04 Panel A: WGCNA module atlas as one integrated heatmap. Rows = modules
# (largest at top); far left = protein-count bars, then each module's colour block
# with its top ORA term (from panel_B), then two bracketed trait sections:
#   baseline (T1) eigengene -> adaptation (LOO-CV) and responder (HR - LR) per timepoint.
pacman::p_load(here, dplyr, tidyr, ggplot2, ggnewscale, patchwork, stringr, scales, grid)
if (!exists("wgcna_mem")) source(here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))

# panel_B attaches AnnotationDbi/clusterProfiler, which mask dplyr verbs (select);
# re-attach dplyr/tidyr last so the unqualified verbs below resolve correctly.
suppressPackageStartupMessages({
  for (p in c("tidyr", "dplyr")) {
    pk <- paste0("package:", p)
    if (pk %in% search()) suppressWarnings(detach(pk, character.only = TRUE))
    library(p, character.only = TRUE)
  }
})

is_light_color <- function(col) {
  rgb <- grDevices::col2rgb(col) / 255
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]
  lum > 0.6
}

mod_sizes <- wgcna_mem |>
  filter(group_id != "grey") |>
  count(group_id, name = "n_proteins") |>
  arrange(desc(n_proteins)) |>
  mutate(mod_col = group_id)
y_levels <- rev(mod_sizes$group_id)
mod_sizes <- mod_sizes |>
  mutate(
    module = factor(group_id, levels = y_levels),
    block_col = ifelse(is_light_color(mod_col), "grey10", "white")
  )

ora_path <- file.path(DAT_DIR, "module_ora.csv")
bio_labels <- if (file.exists(ora_path)) {
  read.csv(ora_path, stringsAsFactors = FALSE) |>
    group_by(group) |>
    slice_min(padj, n = 1, with_ties = FALSE) |>
    ungroup() |>
    transmute(group_id = group, bio = str_wrap(clean_pathway_name(pathway), 14))
} else {
  tibble(group_id = character(), bio = character())
}
mod_sizes <- mod_sizes |>
  left_join(bio_labels, by = "group_id") |>
  mutate(bio = ifelse(is.na(bio), "", bio))

trait_labels <- c(
  comp_hypertrophy = "Composite hypertrophy", d_fcsa_I = "Δ fCSA I",
  d_fcsa_II = "Δ fCSA II",
  d_mcsa = "Δ mCSA", d_1rm_legpress = "Δ 1RM leg-press", d_1rm_ext = "Δ 1RM leg-ext"
)
pred_traits <- names(trait_labels)
resp_times <- c("T1", "T2", "T3")

gap <- 0.6
pred_x <- setNames(seq_along(pred_traits), pred_traits)
resp_x <- setNames(max(pred_x) + gap + seq_along(resp_times), resp_times)
xmin_all <- min(pred_x) - 0.55
xmax_all <- max(resp_x) + 0.55
col_labels <- c(unname(trait_labels), resp_times)
col_breaks <- c(pred_x, resp_x)

L <- max(abs(c(pred$r, resp$r)), na.rm = TRUE)
Lr <- round(L, 1)

cell_txt <- function(r) ifelse(!is.na(r) & abs(r) >= 0.55 * L, "white", "grey10")

pred_long <- pred |>
  mutate(
    module = factor(module, levels = y_levels),
    xpos = pred_x[trait],
    r_lab = sprintf("%.2f", r),
    stars = sig_stars(p_bh),
    border = !is.na(p) & p < 0.05,
    dot = !is.na(q2) & q2 > 0,
    txt_col = cell_txt(r)
  )
resp_long <- resp |>
  mutate(
    module = factor(module, levels = y_levels),
    xpos = resp_x[as.character(timepoint)],
    r_lab = sprintf("%.2f", r),
    stars = sig_stars(p_bh),
    border = !is.na(p) & p < 0.05,
    dot = FALSE,
    txt_col = cell_txt(r)
  )
all_cells <- bind_rows(
  pred_long |> select(module, xpos, r, r_lab, stars, border, dot, txt_col),
  resp_long |> select(module, xpos, r, r_lab, stars, border, dot, txt_col)
)

ttl <- function() element_text(hjust = 0.5, size = FIG_SUBTITLE_SIZE, face = "bold")

p_count <- ggplot(mod_sizes, aes(x = n_proteins, y = module)) +
  geom_col(aes(fill = mod_col), width = 0.85, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_proteins), hjust = -0.15, size = 2.4, fontface = "bold") +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, max(mod_sizes$n_proteins) * 1.28),
    expand = expansion(mult = c(0, 0)), breaks = scales::breaks_pretty(2)
  ) +
  scale_y_discrete(labels = function(x) stringr::str_to_title(x)) +
  labs(title = "Protein count", x = "proteins", y = NULL) +
  FIG_THEME +
  theme(
    axis.text.y = element_text(face = "bold", size = FIG_AXIS_TEXT),
    axis.ticks.y = element_blank(), panel.grid = element_blank(),
    plot.title = ttl(), plot.margin = margin(2, 2, 2, 2)
  )

p_block <- ggplot(mod_sizes, aes(x = 1, y = module)) +
  geom_tile(aes(fill = mod_col), width = 1, height = 0.9, color = "black", linewidth = 0.3) +
  geom_text(aes(label = bio, color = block_col),
    size = 2.0, fontface = "bold", lineheight = 0.82
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  scale_x_continuous(expand = c(0, 0)) +
  labs(title = "Module biology", x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    panel.grid = element_blank(), panel.border = element_blank(),
    plot.title = ttl(), plot.margin = margin(2, 2, 2, 0)
  )

p_heat <- ggplot() +
  geom_tile(
    data = pred_long, aes(xpos, module, fill = r),
    width = 1, height = 1, color = "grey85", linewidth = 0.3
  ) +
  scale_fill_gradient2(
    low = GROUP_COLORS[["HR"]], mid = "white", high = GROUP_COLORS[["LR"]],
    midpoint = 0, limits = c(-L, L), breaks = c(-Lr, 0, Lr),
    name = "r  baseline → adaptation  "
  ) +
  new_scale_fill() +
  geom_tile(
    data = resp_long, aes(xpos, module, fill = r),
    width = 1, height = 1, color = "grey85", linewidth = 0.3
  ) +
  scale_fill_gradient2(
    low = DIR_COLORS[["Down"]], mid = "white", high = DIR_COLORS[["Up"]],
    midpoint = 0, limits = c(-L, L), breaks = c(-Lr, 0, Lr),
    name = "r  responder (HR − LR)  "
  ) +
  geom_tile(
    data = filter(all_cells, border), aes(xpos, module),
    width = 1, height = 1, fill = NA, color = "black", linewidth = 0.9
  ) +
  geom_point(
    data = filter(pred_long, dot), aes(xpos, module),
    shape = 21, fill = "white", color = "black", size = 1.6, stroke = 0.5,
    position = position_nudge(y = -0.30)
  ) +
  geom_text(
    data = all_cells, aes(xpos, module, label = r_lab, color = txt_col),
    size = 2.3, fontface = "bold"
  ) +
  geom_text(
    data = filter(all_cells, stars != ""),
    aes(xpos, module, label = stars, color = txt_col),
    size = 2.7, fontface = "bold", position = position_nudge(y = 0.30)
  ) +
  scale_color_identity() +
  scale_x_continuous(
    breaks = col_breaks, labels = col_labels,
    limits = c(xmin_all, xmax_all), expand = c(0, 0)
  ) +
  scale_y_discrete(labels = NULL) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT - 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), legend.position = "bottom",
    plot.margin = margin(2, 4, 2, 2)
  ) +
  guides(fill = guide_colorbar(
    barwidth = unit(34, "mm"), barheight = unit(3.5, "mm"),
    title.position = "top", title.hjust = 0.5
  ))

brackets <- tibble(
  label = c("Baseline (T1) → adaptation", "Responder (HR − LR)"),
  start = c(pred_x[[1]], resp_x[[1]]),
  end   = c(pred_x[[length(pred_x)]], resp_x[[length(resp_x)]])
) |>
  mutate(mid = (start + end) / 2)

p_brackets <- ggplot(brackets) +
  geom_segment(aes(x = start - 0.4, xend = end + 0.4, y = 0.2, yend = 0.2),
    linewidth = 0.4, color = "grey30"
  ) +
  geom_segment(aes(x = start - 0.4, xend = start - 0.4, y = 0.2, yend = 0.1),
    linewidth = 0.4, color = "grey30"
  ) +
  geom_segment(aes(x = end + 0.4, xend = end + 0.4, y = 0.2, yend = 0.1),
    linewidth = 0.4, color = "grey30"
  ) +
  geom_text(aes(x = mid, y = 0.55, label = label),
    size = 2.9, fontface = "bold", color = "grey20"
  ) +
  scale_x_continuous(limits = c(xmin_all, xmax_all), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  theme_void() +
  theme(plot.margin = margin(2, 4, 0, 2))

design <- c(
  area(2, 1, 11, 2),
  area(2, 3, 11, 4),
  area(2, 5, 11, 12),
  area(1, 5, 1, 12)
)
# Title, subtitle, and tag are supplied by the composite (HRvLR_F04_composite.R),
# so the panel renders header-free to sit cleanly under the composite's draw_labels.
panel_a <- p_count + p_block + p_heat + p_brackets +
  plot_layout(design = design, guides = "collect") &
  theme(legend.position = "bottom")

save_panel(panel_a, file.path(RPT_DIR, "panels", "panel_a_wgcna_map"), 400, 280)
cat("F04 Panel A done.\n")
