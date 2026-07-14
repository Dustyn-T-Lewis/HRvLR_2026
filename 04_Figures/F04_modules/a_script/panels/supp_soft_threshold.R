# F04 supplement: WGCNA soft-threshold power sweep. Scale-free topology R^2 and
# mean connectivity across candidate powers; the chosen beta is highlighted.
pacman::p_load(here, dplyr, tibble, ggplot2, ggrepel, patchwork)
if (!exists("wgcna_sft")) source(here("04_Figures", "F04", "a_script", "HRvLR_F04_setup.R"))

fit <- wgcna_sft$fitIndices
fit_df <- tibble(
  power = fit$Power,
  r2 = -sign(fit$slope) * fit$SFT.R.sq,
  mean_k = fit$mean.k.
) |>
  mutate(selected = power == wgcna_power)

sel_r2 <- fit_df$r2[fit_df$power == wgcna_power]
sft_r2_threshold <- 0.87

p1 <- ggplot(fit_df, aes(power, r2)) +
  geom_hline(yintercept = sft_r2_threshold, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = selected), shape = 21, size = 2.5, stroke = 0.4, color = "black") +
  ggrepel::geom_text_repel(aes(label = power),
    size = 2.6, color = "grey25", fontface = "bold",
    max.overlaps = 20, seed = 42, min.segment.length = 0.3,
    segment.size = 0.25, segment.color = "grey60"
  ) +
  scale_fill_manual(values = c("TRUE" = DIR_COLORS[["Up"]], "FALSE" = "grey70"), guide = "none") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  annotate("text",
    x = 19, y = 0.90,
    label = as.expression(bquote(R^2 == .(sft_r2_threshold) ~ threshold)),
    size = 2.4, color = "grey40", hjust = 1
  ) +
  labs(
    x = "Soft-threshold (power)", y = expression(Scale - free ~ topology ~ R^2),
    title = "Scale independence"
  ) +
  FIG_THEME +
  theme(plot.title = element_text(size = 10, face = "bold"))

p2 <- ggplot(fit_df, aes(power, mean_k)) +
  geom_point(aes(fill = selected), shape = 21, size = 2.5, stroke = 0.4, color = "black") +
  ggrepel::geom_text_repel(aes(label = power),
    size = 2.6, color = "grey25", fontface = "bold",
    max.overlaps = 20, seed = 42, min.segment.length = 0.3,
    segment.size = 0.25, segment.color = "grey60"
  ) +
  scale_fill_manual(values = c("TRUE" = DIR_COLORS[["Up"]], "FALSE" = "grey70"), guide = "none") +
  scale_x_continuous(breaks = seq(2, 20, 2)) +
  labs(x = "Soft-threshold (power)", y = "Mean connectivity", title = "Mean connectivity") +
  FIG_THEME +
  theme(plot.title = element_text(size = 10, face = "bold"))

supp_sft <- (p1 | p2) +
  plot_annotation(
    title = "WGCNA soft-threshold selection",
    subtitle = sprintf(
      "Signed network; selected power = %d (R² = %.3f); %d proteins",
      wgcna_power, sel_r2, nrow(wgcna_mem)
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(face = "bold.italic", size = 10, color = "grey30")
    )
  )

dir.create(file.path(RPT_DIR, "supp"), recursive = TRUE, showWarnings = FALSE)
save_panel(supp_sft, file.path(RPT_DIR, "supp", "supp_soft_threshold"), 240, 110)
message("F04 supplement (soft-threshold) done.")
