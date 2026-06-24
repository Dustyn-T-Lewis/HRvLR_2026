# F03 Panel B: WGCNA module eigengene trajectories ----
setwd(rprojroot::find_rstudio_root_file())
if (!exists("wgcna_eig")) source("04_Figures/F03/a_script/HRvLR_F03_setup.R")
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

# Ordered factors: top-size-first module order, timepoints T1 -> T3
mod_order <- c("turquoise", "blue", "brown", "yellow", "green")

plot_eig <- wgcna_eig |>
  mutate(
    timepoint = factor(timepoint, levels = c("T1", "T2", "T3")),
    group_id  = factor(group_id, levels = mod_order)
  )

p <- ggplot(plot_eig, aes(x = timepoint, y = ME)) +
  geom_line(
    aes(group = subject, colour = group_arm),
    alpha = 0.20, linewidth = 0.3
  ) +
  stat_summary(
    aes(group = group_arm, colour = group_arm),
    fun = mean, geom = "line", linewidth = 1.1
  ) +
  stat_summary(
    aes(group = group_arm, colour = group_arm),
    fun = mean, geom = "point", size = 2
  ) +
  scale_colour_manual(values = GROUP_COLORS, name = NULL) +
  facet_wrap(~group_id, nrow = 1) +
  labs(
    title    = "Module eigengene trajectories across training and acute exercise",
    subtitle = "Mean +/- per-subject lines by responder group; T1 baseline, T2 trained, T3 acute (n approx 16)",
    x        = NULL,
    y        = "Module eigengene (PC1)",
    tag      = "B"
  ) +
  FIG_THEME

save_panel(p, file.path(RPT_DIR, "panel_b_trajectories"), 230, 80)

# Audit: mean ME per module x timepoint x group_arm
audit_means <- plot_eig |>
  group_by(group_id, timepoint, group_arm) |>
  summarise(mean_ME = mean(ME, na.rm = TRUE), .groups = "drop")

write.csv(audit_means, file.path(DAT_DIR, "audit_panel_B_eigengene_means.csv"), row.names = FALSE)

cat("F03 Panel B done.\n")
