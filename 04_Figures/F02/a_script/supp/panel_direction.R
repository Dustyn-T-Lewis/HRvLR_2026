# F02 Supp: response-trajectory direction (Δ-score PCA + RRPP vector test)
# Geometry of the per-subject change vectors, training and acute, with the overall
# change-vector difference tested by RRPP and the angle between the two group-mean
# vectors. Companion to Panel B (magnitude); both are null.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, patchwork, RRPP)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

supp_dir <- file.path(RPT_DIR, "supp")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

PD_W <- 200
PD_H <- PANEL_H

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

vector_angle <- function(dm, grp) {
  mu <- vapply(levels(grp), function(g) colMeans(dm[grp == g, , drop = FALSE]), numeric(ncol(dm)))
  cos_ab <- sum(mu[, 1] * mu[, 2]) / (sqrt(sum(mu[, 1]^2)) * sqrt(sum(mu[, 2]^2)))
  acos(pmin(pmax(cos_ab, -1), 1)) * 180 / pi
}

phase_levels <- c("Training (T2 - T1)", "Acute (T3 - T2)")
phase_defs <- list(Training = c("T1", "T2"), Acute = c("T2", "T3"))

dir_fit <- function(phase) {
  dm <- delta_matrix(phase_defs[[phase]][1], phase_defs[[phase]][2])
  grp <- droplevels(subj_group[rownames(dm)])
  set.seed(42)
  vec_tab <- anova(lm.rrpp(dm ~ grp, data = rrpp.data.frame(dm = dm, grp = grp), iter = 999))$table
  pca <- prcomp(dm, center = TRUE, scale. = FALSE)
  list(
    scores = as_tibble(pca$x[, 1:2]) |> mutate(Group = grp, phase = phase),
    stat = tibble(phase = phase, label = sprintf("vector %s\nangle = %.0f°", fmt_p(vec_tab$`Pr(>F)`[1]), vector_angle(dm, grp)))
  )
}

fits <- lapply(c("Training", "Acute"), dir_fit)
relabel <- function(x) factor(recode(x, Training = phase_levels[1], Acute = phase_levels[2]), levels = phase_levels)
score_df <- bind_rows(lapply(fits, `[[`, "scores")) |> mutate(phase = relabel(phase))
stat_df <- bind_rows(lapply(fits, `[[`, "stat")) |> mutate(phase = relabel(phase))

pDir <- ggplot(score_df, aes(PC1, PC2, color = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group), level = 0.80, linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.7, alpha = 0.9) +
  geom_text(
    data = stat_df, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.15, size = FIG_GEOM_TEXT, fontface = "bold", color = "grey25"
  ) +
  facet_wrap(~phase, scales = "free") +
  scale_color_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    title = "Response Direction (Δ-score PCA)",
    subtitle = "geometry of the change vectors | RRPP vector test",
    x = "PC1", y = "PC2", tag = "O"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

save_png(pDir, file.path(supp_dir, "panel_o_direction"), PD_W, PD_H)
write.csv(score_df, file.path(DAT_DIR, "audit_panel_O_direction.csv"), row.names = FALSE)
cat("F02 Supp Panel O (direction) done.\n")
