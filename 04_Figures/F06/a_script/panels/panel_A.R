# F06 Panel A: WGCNA module atlas - one row per module, three aligned blocks:
#   (1) module size, (2) baseline (T1) prediction of each training adaptation
#   (leave-one-out cross-validated), (3) responder (HR vs LR) signal per timepoint.
pacman::p_load(here, dplyr, tidyr, ggplot2, patchwork, scales)
if (!exists("wgcna_mem")) source(here("04_Figures", "F06", "a_script", "HRvLR_F06_setup.R"))

# Modules ordered by their best out-of-sample prediction (predictive at top)
mod_rank <- pred |>
  group_by(group_id = module) |>
  summarise(best_q2 = max(q2, na.rm = TRUE), .groups = "drop")
mod_sizes <- wgcna_mem |>
  filter(group_id != "grey") |>
  count(group_id, name = "n_proteins") |>
  left_join(mod_rank, by = "group_id") |>
  arrange(desc(best_q2)) |>
  mutate(mod_col = group_id, module = factor(group_id, levels = group_id))
mod_order <- levels(mod_sizes$module)

trait_labels <- c(
  comp_hypertrophy = "Composite hypertrophy", d_fcsa_I = "Δ fCSA I",
  d_fcsa_II = "Δ fCSA II", d_myovision_fcsa_I = "Δ MyoVision fCSA I",
  d_mcsa = "Δ mCSA", d_1rm_legpress = "Δ 1RM leg-press", d_1rm_ext = "Δ 1RM leg-ext"
)
pred_df <- pred |>
  mutate(
    module = factor(module, levels = mod_order),
    trait_label = factor(trait_labels[trait], levels = unname(trait_labels))
  )
resp_df <- resp |>
  mutate(
    module = factor(module, levels = mod_order),
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3"))
  )

L <- max(abs(c(pred_df$r, resp_df$r)), na.rm = TRUE)
predicts <- filter(pred_df, !is.na(q2) & q2 > 0)
pred_sig <- filter(pred_df, p < 0.05)
resp_sig <- filter(resp_df, p < 0.05)
strip_y <- theme(strip.text.y = element_blank(), strip.background.y = element_blank())
ttl <- function() element_text(hjust = 0.5, size = FIG_SUBTITLE_SIZE, face = "bold")

# (1) module size
p_count <- ggplot(mod_sizes, aes(x = n_proteins, y = "")) +
  geom_col(aes(fill = mod_col), width = 1, color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_proteins), hjust = -0.18, size = 2.5, fontface = "bold") +
  facet_grid(module ~ ., switch = "y") +
  scale_fill_identity() +
  scale_x_continuous(
    limits = c(0, max(mod_sizes$n_proteins) * 1.25),
    expand = expansion(mult = c(0, 0)), breaks = scales::breaks_pretty(2)
  ) +
  labs(title = "Protein count", x = "proteins", y = NULL) +
  FIG_THEME +
  theme(
    strip.text.y.left = element_text(angle = 0, face = "bold", size = FIG_AXIS_TEXT),
    strip.placement = "outside", strip.background = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(1.5, "pt"), plot.title = ttl()
  )

# (2) baseline prediction of adaptation: fill = r, dot = predicts OOS, bold = p<0.05
p_pred <- ggplot(pred_df, aes(trait_label, y = "", fill = r)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  geom_tile(data = pred_sig, fill = NA, color = "black", linewidth = 0.9) +
  geom_point(
    data = predicts, aes(trait_label, ""), inherit.aes = FALSE,
    shape = 21, fill = "white", color = "black", size = 1.7, stroke = 0.5
  ) +
  facet_grid(module ~ .) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-L, L), breaks = c(-round(L, 1), 0, round(L, 1)), name = "Correlation r  "
  ) +
  labs(title = "Baseline (T1) module → adaptation", x = NULL, y = NULL) +
  FIG_THEME +
  strip_y +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = FIG_AXIS_TEXT - 1),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(1.5, "pt"),
    legend.position = "bottom", plot.title = ttl()
  )

# (3) responder (HR - LR) signal per timepoint
p_resp <- ggplot(resp_df, aes(timepoint, y = "", fill = r)) +
  geom_tile(color = "grey85", linewidth = 0.3) +
  geom_tile(data = resp_sig, fill = NA, color = "black", linewidth = 0.9) +
  facet_grid(module ~ .) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-L, L), guide = "none"
  ) +
  labs(title = "Responder (HR − LR)", x = NULL, y = NULL) +
  FIG_THEME +
  strip_y +
  theme(
    axis.text.x = element_text(size = FIG_AXIS_TEXT),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    panel.grid = element_blank(), panel.spacing = unit(1.5, "pt"), plot.title = ttl()
  )

panel_a <- p_count + p_pred + p_resp +
  plot_layout(widths = c(0.8, 3, 1.1), guides = "collect") +
  plot_annotation(
    title = "WGCNA module atlas: size, baseline prediction of adaptation, and responder signal",
    subtitle = sprintf(
      "Full-proteome signed network: %d proteins in %d modules (rows ordered by best out-of-sample prediction). Middle: correlation of each module's baseline (T1) eigengene with each training adaptation (Δ T1→T2); white dot = predicts out-of-sample (leave-one-out Q²>0), bold = in-sample p<0.05. Right: module–responder (HR−LR) correlation per timepoint. The pink module's baseline predicts Δ mCSA out-of-sample.",
      sum(mod_sizes$n_proteins), length(mod_order)
    ),
    tag_levels = list(c("A", "", "")),
    theme = theme(
      plot.title = element_text(face = "bold", size = FIG_TITLE_SIZE),
      plot.subtitle = element_text(face = "italic", size = FIG_SUBTITLE_SIZE - 1, color = "grey30"),
      plot.tag = element_text(face = "bold", size = FIG_TAG_SIZE)
    )
  ) &
  theme(legend.position = "bottom")

save_panel(panel_a, file.path(RPT_DIR, "panels", "panel_a_wgcna_map"), 270, 200)
cat("F06 Panel A done.\n")
