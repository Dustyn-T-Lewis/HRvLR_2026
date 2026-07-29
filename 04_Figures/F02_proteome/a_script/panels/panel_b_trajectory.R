# F02 Panel B: response-trajectory magnitude (RRPP permutation)
# Per-subject change-vector length ||delta|| (training T2 - T1, acute T3 - T2), HR vs LR,
# tested with RRPP (Collyer, Sekora & Adams 2015). How far the proteome moved; the
# direction geometry is in the supplement. Title drawn on the composite.

pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, RRPP)

if (!exists("meta")) source(here("04_Figures", "F02_proteome", "a_script", "setup.R"))

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
    stat = tibble(phase = phase, p = tab$`Pr(>F)`[1], label = sprintf("RRPP %s", fmt_p(tab$`Pr(>F)`[1])))
  )
}

fits <- lapply(c("Training", "Acute"), mag_fit)
relabel <- function(x) factor(recode(x, Training = phase_levels[1], Acute = phase_levels[2]), levels = phase_levels)
mag_df <- bind_rows(lapply(fits, `[[`, "mag")) |> mutate(phase = relabel(phase))
mag_stat <- bind_rows(lapply(fits, `[[`, "stat")) |> mutate(phase = relabel(phase))

pB <- ggplot(mag_df, aes(Group, magnitude, color = Group, fill = Group)) +
  geom_boxplot(alpha = 0.18, outlier.shape = NA, width = 0.55, linewidth = 0.4) +
  # Seeded: geom_jitter draws its offsets at render time, so an unseeded jitter
  # redraws the points on every ggsave and the panel never reproduces. RRPP leaves
  # the RNG in a state we do not control, so the seed must live on the layer.
  geom_point(
    position = position_jitter(width = 0.12, height = 0, seed = 42),
    size = 1.6, alpha = 0.85
  ) +
  geom_text(
    data = mag_stat, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
    vjust = 1.4, size = FIG_GEOM_TEXT, fontface = "bold", color = "grey25"
  ) +
  facet_wrap(~phase) +
  scale_color_manual(values = GROUP_COLORS, guide = "none") +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(x = NULL, y = "||delta|| (log2 units)") +
  FIG_THEME +
  theme(legend.position = "none")

save_png(pB, file.path(RPT_DIR, "panels", "panel_b_trajectory"), PB_W, PB_H)
F02_AUDIT[["panel_B_trajectory"]] <- mag_df
cat("F02 Panel B done.\n")

# --- Supplement (owned by this panel): response-trajectory DIRECTION.
# Panel B reads the magnitude of the change vectors; this reads their geometry.
pacman::p_load(here, dplyr, tidyr, tibble, ggplot2, patchwork, RRPP)


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
    stat = tibble(phase = phase, p = vec_tab$`Pr(>F)`[1], label = sprintf("vector %s\nangle = %.0f°", fmt_p(vec_tab$`Pr(>F)`[1]), vector_angle(dm, grp)))
  )
}

fits <- lapply(c("Training", "Acute"), dir_fit)
relabel <- function(x) factor(recode(x, Training = phase_levels[1], Acute = phase_levels[2]), levels = phase_levels)
dir_score <- bind_rows(lapply(fits, `[[`, "scores")) |> mutate(phase = relabel(phase))
dir_stat <- bind_rows(lapply(fits, `[[`, "stat")) |> mutate(phase = relabel(phase))

pDir <- ggplot(dir_score, aes(PC1, PC2, color = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group), level = 0.80, linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.7, alpha = 0.9) +
  geom_text(
    data = dir_stat, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
    hjust = -0.05, vjust = 1.15, size = FIG_GEOM_TEXT, fontface = "bold", color = "grey25"
  ) +
  facet_wrap(~phase, scales = "free") +
  scale_color_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    title = "Response Direction (Δ-score PCA)",
    subtitle = "geometry of the change vectors -- RRPP vector test",
    x = "PC1", y = "PC2", tag = "O"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

save_png(pDir, file.path(supp_dir, "panel_o_direction"), PD_W, PD_H)
F02_AUDIT[["panel_O_direction"]] <- dir_score
cat("F02 Supp Panel O (direction) done.\n")
