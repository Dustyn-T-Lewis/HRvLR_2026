# F02 Panel C: CV% Violins (Inter-Individual Variability) ----
# CV% per protein within each Group_Time, violin + jitter

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

library(dplyr)
library(tidyr)
library(ggplot2)

PC_W <- 140
PC_H <- 120
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

GROUP_ORDER <- c("HR_T1", "HR_T2", "HR_T3", "LR_T1", "LR_T2", "LR_T3")
GROUP_COLS <- GROUP_FILL

# CV on linear scale per protein within each group
lin_mat <- 2^as.matrix(imp_mat)
cv_list <- lapply(GROUP_ORDER, function(g) {
  idx <- meta$Col_ID[meta$Group_Time == g]
  idx <- intersect(idx, colnames(lin_mat))
  if (length(idx) < 2) {
    return(NULL)
  }
  sub <- lin_mat[, idx, drop = FALSE]
  cv <- apply(sub, 1, function(x) sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100)
  tibble(cv = cv, group = g)
})
cv_df <- bind_rows(cv_list) |>
  mutate(group = factor(group, levels = GROUP_ORDER))

cv_med <- cv_df |>
  group_by(group) |>
  summarise(med = median(cv, na.rm = TRUE), .groups = "drop")

pC <- ggplot(cv_df, aes(x = group, y = cv, fill = group)) +
  geom_violin(alpha = 0.4, color = NA, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", alpha = 0.6) +
  geom_text(
    data = cv_med, aes(x = group, y = med, label = sprintf("%.1f%%", med)),
    vjust = -0.5, size = 3, fontface = "bold"
  ) +
  scale_fill_manual(values = GROUP_COLS) +
  labs(
    title = "Inter-Individual Variability",
    subtitle = "CV% per protein (linear scale)",
    x = NULL, y = "CV%", tag = "C"
  ) +
  FIG_THEME +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 30, hjust = 1)
  )

save_panel(pC, file.path(RPT_DIR, "panel_c_cv"), PC_W, PC_H)
write.csv(cv_med, file.path(DAT_DIR, "audit_panel_C_cv_medians.csv"), row.names = FALSE)
cat("F02 Panel C done.\n")
