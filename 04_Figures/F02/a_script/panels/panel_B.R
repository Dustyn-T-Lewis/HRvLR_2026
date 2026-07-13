# F02 Panel B: response-trajectory magnitude (RRPP permutation)
# Per-subject change-vector length ‖Δ‖ (training T2 - T1, acute T3 - T2), HR vs LR,
# tested with RRPP (Collyer, Sekora & Adams 2015). How far the proteome moved; the
# direction geometry is in the supplement. Title drawn on the composite.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, RRPP)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

PB_W <- 120
PB_H <- 110

subject_of <- function(ids) sub("_T[123]$", "", ids)
subj_group <- setNames(meta$Group, subject_of(meta$Col_ID))

delta_matrix <- function(from_tp, to_tp) {
  from_ids <- meta$Col_ID[meta$Timepoint == from_tp & meta$Col_ID %in% colnames(imp_mat)]
  to_ids <- meta$Col_ID[meta$Timepoint == to_tp & meta$Col_ID %in% colnames(imp_mat)]
  paired <- intersect(subject_of(from_ids), subject_of(to_ids))
  d <- t(imp_mat[, to_ids[match(paired, subject_of(to_ids))]]) -
    t(imp_mat[, from_ids[match(paired, subject_of(from_ids))]])
  rownames(d) <- paired
  d[, apply(d, 2, var) > 0]
}

phase_levels <- c("Training (T2 - T1)", "Acute (T3 - T2)")
phase_defs <- list(Training = c("T1", "T2"), Acute = c("T2", "T3"))

mag_fit <- function(phase) {
  dm <- delta_matrix(phase_defs[[phase]][1], phase_defs[[phase]][2])
  grp <- droplevels(subj_group[rownames(dm)])
  mag <- sqrt(rowSums(dm^2))
  set.seed(42)
  tab <- anova(lm.rrpp(mag ~ grp, data = rrpp.data.frame(mag = mag, grp = grp), iter = 999))$table
  list(
    mag = tibble(subject = rownames(dm), phase = phase, Group = grp, magnitude = mag),
    stat = tibble(phase = phase, label = sprintf("RRPP %s", fmt_p(tab$`Pr(>F)`[1])))
  )
}

fits <- lapply(c("Training", "Acute"), mag_fit)
relabel <- function(x) factor(recode(x, Training = phase_levels[1], Acute = phase_levels[2]), levels = phase_levels)
mag_df <- bind_rows(lapply(fits, `[[`, "mag")) |> mutate(phase = relabel(phase))
stat_df <- bind_rows(lapply(fits, `[[`, "stat")) |> mutate(phase = relabel(phase))

pB <- ggplot(mag_df, aes(Group, magnitude, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.18, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  geom_jitter(width = 0.12, height = 0, size = 1.6, alpha = 0.85) +
  geom_text(
    data = stat_df, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
    vjust = 1.4, size = FIG_GEOM_TEXT, fontface = "bold", color = "grey25"
  ) +
  facet_wrap(~phase) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(x = NULL, y = "‖Δ‖ (log2 units)") +
  FIG_THEME +
  theme(legend.position = "none")

save_png(pB, file.path(RPT_DIR, "panels", "panel_b_trajectory"), PB_W, PB_H)
write.csv(mag_df, file.path(DAT_DIR, "audit_panel_B_trajectory.csv"), row.names = FALSE)
cat("F02 Panel B done.\n")
