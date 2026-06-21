# Figure 3 — Panel B: DEP Rank Location (Barcode Plot) ----
# Shows where Pi < 0.05 DEPs sit in the t-statistic-ranked proteome.
# Dual density traces (Up/Down) + direction-colored barcode ticks.
# 5 main contrasts from 03_combined_results.csv
# Outputs: pB (ggplot object), panel_B_barcode_MAIN.pdf/.png

if (!exists("meta")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")

RPT <- RPT_DIR
DAT <- DAT_DIR

CONTRASTS_B <- c("Training_HR", "Training_LR", "Acute_HR",
                 "Acute_LR", "Baseline_HRvLR")

PB_W <- 170
PB_H <- 200

# --- Build long-form data: rank position + DEP status per contrast ----------

rank_list <- lapply(CONTRASTS_B, function(ctr) {
  t_col   <- paste0("t_", ctr)
  pi_col  <- paste0("pi_score_", ctr)
  lfc_col <- paste0("logFC_", ctr)

  dep_df |>
    filter(!is.na(.data[[t_col]])) |>
    arrange(.data[[t_col]]) |>
    mutate(
      rank_frac = seq_len(n()) / n(),
      is_dep    = !is.na(.data[[pi_col]]) & .data[[pi_col]] < 0.05,
      direction = case_when(
        !is_dep              ~ NA_character_,
        .data[[lfc_col]] > 0 ~ "Up",
        TRUE                 ~ "Down"
      ),
      contrast = ctr
    ) |>
    select(gene, contrast, rank_frac, is_dep, direction)
})
rank_df <- bind_rows(rank_list)
rank_df$contrast <- factor(rank_df$contrast, levels = CONTRASTS_B)

# DEP subset for ticks and density
dep_only <- rank_df |> filter(is_dep)
dep_only$direction <- factor(dep_only$direction, levels = c("Up", "Down"))

# DEP counts per contrast for annotation
dep_counts <- dep_only |>
  group_by(contrast) |>
  summarise(
    n_up    = sum(direction == "Up"),
    n_down  = sum(direction == "Down"),
    n_total = n(),
    .groups = "drop"
  ) |>
  mutate(label = sprintf("n = %d  (%d\u2191  %d\u2193)", n_total, n_up, n_down))

write.csv(dep_counts, file.path(DAT, "panel_B_barcode_enrichment.csv"),
          row.names = FALSE)

# --- Pre-compute density curves (normalize peak = 1.0) ----------------------

DENS_PAD <- 0.06
dens_list <- lapply(split(dep_only, dep_only$contrast), function(ctr_df) {
  lapply(split(ctr_df, ctr_df$direction, drop = TRUE), function(dir_df) {
    if (nrow(dir_df) < 2) return(NULL)
    d <- density(dir_df$rank_frac, adjust = 1.8,
                 from = -DENS_PAD, to = 1 + DENS_PAD, n = 512)
    tibble(x = d$x, y = d$y, direction = dir_df$direction[1],
           contrast = dir_df$contrast[1])
  }) |> bind_rows()
}) |> bind_rows()

dens_list <- dens_list |>
  group_by(contrast) |>
  mutate(y_norm = y / max(y)) |>
  ungroup()
dens_list$direction <- factor(dens_list$direction, levels = c("Up", "Down"))
dens_list$contrast  <- factor(dens_list$contrast, levels = CONTRASTS_B)

TICK_DEPTH <- -0.40
ANNOT_SZ   <- scale_text(BASE_STAT - 0.8, PB_W)

# --- Peak positions for label placement -------------------------------------

peak_pos <- dens_list |>
  group_by(contrast, direction) |>
  slice_max(y_norm, n = 1, with_ties = FALSE) |>
  ungroup() |>
  select(contrast, direction, peak_x = x, peak_y = y_norm)

LABEL_NUDGE <- 0.03

n_down <- dep_only |> filter(direction == "Down") |> count(contrast) |>
  tibble::deframe()
n_up   <- dep_only |> filter(direction == "Up")   |> count(contrast) |>
  tibble::deframe()

DESC_DOWN <- c(Training_HR    = "dec. with training (HR)",
               Training_LR    = "dec. with training (LR)",
               Acute_HR       = "dec. acutely (HR)",
               Acute_LR       = "dec. acutely (LR)",
               Baseline_HRvLR = "lower in HR vs LR")
DESC_UP   <- c(Training_HR    = "inc. with training (HR)",
               Training_LR    = "inc. with training (LR)",
               Acute_HR       = "inc. acutely (HR)",
               Acute_LR       = "inc. acutely (LR)",
               Baseline_HRvLR = "higher in HR vs LR")

annot_down <- peak_pos |>
  filter(direction == "Down") |>
  mutate(
    label_x = peak_x + LABEL_NUDGE,
    label_y = peak_y,
    label = paste(n_down[as.character(contrast)],
                  DESC_DOWN[as.character(contrast)])
  )

annot_up <- peak_pos |>
  filter(direction == "Up") |>
  mutate(
    label_x = peak_x - LABEL_NUDGE,
    label_y = peak_y,
    label = paste(n_up[as.character(contrast)],
                  DESC_UP[as.character(contrast)])
  )

conn_down <- annot_down |>
  mutate(x_start = peak_x, y_start = peak_y,
         x_end = label_x, y_end = label_y)
conn_up <- annot_up |>
  mutate(x_start = peak_x, y_start = peak_y,
         x_end = label_x, y_end = label_y)

# --- Plot --------------------------------------------------------------------

pB <- ggplot() +
  # Background contrast wash
  geom_rect(data = tibble(contrast = factor(CONTRASTS_B, levels = CONTRASTS_B),
                           fill = unname(CONTRAST_COLORS[CONTRASTS_B])),
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            fill = rep(unname(CONTRAST_COLORS[CONTRASTS_B]), 1),
            alpha = 0.10) +
  # Density traces
  geom_ribbon(data = dens_list,
              aes(x = x, ymin = 0, ymax = y_norm, fill = direction),
              alpha = 0.30, outline.type = "upper") +
  geom_line(data = dens_list,
            aes(x = x, y = y_norm, color = direction),
            linewidth = 0.5) +
  # Barcode ticks
  geom_segment(data = dep_only,
               aes(x = rank_frac, xend = rank_frac,
                   y = 0, yend = TICK_DEPTH,
                   color = direction),
               linewidth = 0.35, alpha = 0.8) +
  # Zero line
  geom_hline(yintercept = 0, linewidth = 0.25, color = "grey50") +
  # Connector lines
  geom_segment(data = conn_down,
               aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
               linewidth = 0.3, color = unname(DIR_COLORS["Down"]),
               alpha = 0.4, linetype = "solid") +
  geom_segment(data = conn_up,
               aes(x = x_start, xend = x_end, y = y_start, yend = y_end),
               linewidth = 0.3, color = unname(DIR_COLORS["Up"]),
               alpha = 0.4, linetype = "solid") +
  # Down labels (blue box)
  geom_label(data = annot_down,
             aes(x = label_x, y = label_y, label = label),
             hjust = 0, vjust = 0.5,
             size = ANNOT_SZ,
             color = "white", fontface = "bold",
             fill = scales::alpha(unname(DIR_COLORS["Down"]), 0.85),
             linewidth = 0, label.padding = unit(1.2, "mm"),
             label.r = unit(0.8, "mm")) +
  # Up labels (red box)
  geom_label(data = annot_up,
             aes(x = label_x, y = label_y, label = label),
             hjust = 1, vjust = 0.5,
             size = ANNOT_SZ,
             color = "white", fontface = "bold",
             fill = scales::alpha(unname(DIR_COLORS["Up"]), 0.85),
             linewidth = 0, label.padding = unit(1.2, "mm"),
             label.r = unit(0.8, "mm")) +
  # Scales
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_color_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                 Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.08))) +
  coord_cartesian(xlim = c(-DENS_PAD, 1 + DENS_PAD),
                  ylim = c(TICK_DEPTH, 1.15)) +
  # Facet: one row per contrast
  facet_wrap(~ contrast, ncol = 1, strip.position = "left",
             labeller = as_labeller(CTR_FACET)) +
  labs(
    title    = "DEP Rank Location",
    subtitle = "\u03A0 \u2264 0.05 DEPs in t-stat-ranked proteome",
    x        = "Rank position (by t-statistic)",
    y        = NULL,
    tag      = "B"
  ) +
  FIG_THEME +
  theme(
    legend.position    = "none",
    axis.text.y        = element_blank(),
    axis.ticks.y       = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    strip.text.y.left  = element_text(face = "bold", size = FIG_STRIP_SIZE,
                                       angle = 0, hjust = 1, vjust = 0.5),
    strip.placement    = "outside",
    panel.spacing.y    = unit(1, "mm")
  )

save_panel(pB, file.path(RPT, "panel_B_barcode_MAIN"), PB_W, PB_H)
message("F03 Panel B (barcode) done.")
