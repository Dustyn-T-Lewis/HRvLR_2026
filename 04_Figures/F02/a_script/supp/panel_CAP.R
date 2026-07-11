# F02 Supp Panel N: constrained ordination (CAP) with cross-validated allocation
# The unsupervised PCA does not separate responders. CAP actively seeks the axis that
# best discriminates HR from LR on the sample distance matrix, conditioned on timepoint,
# but stays honest via a permutation test and leave-one-out misclassification
# (Anderson & Willis 2003). A near-chance LOO rate is the legitimate negative that
# corroborates the null PCA.

pacman::p_load(here, dplyr, tibble, ggplot2, vegan, MASS)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

supp_dir <- file.path(RPT_DIR, "supp")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

PN_W <- 150
PN_H <- PANEL_H

ids <- meta$Col_ID[meta$Col_ID %in% colnames(imp_mat)]
md <- meta[match(ids, meta$Col_ID), ]
d <- dist(t(imp_mat[, ids]))

cap <- capscale(d ~ Group + Condition(Timepoint), data = md)
set.seed(42)
cap_p <- anova(cap, permutations = 999)$`Pr(>F)`[1]

n_axes <- min(10, length(ids) - 2)
pco <- cmdscale(d, k = n_axes)
loo <- lda(pco, grouping = md$Group, CV = TRUE)
loo_err <- mean(loo$class != md$Group)

sc <- as.data.frame(scores(cap, display = "sites", choices = 1:2, scaling = "sites"))
names(sc) <- c("CAP1", "MDS1")
sc$Group <- md$Group

pN <- ggplot(sc, aes(CAP1, MDS1, color = Group)) +
  stat_ellipse(aes(fill = Group), geom = "polygon", alpha = 0.10, level = 0.80, show.legend = FALSE) +
  stat_ellipse(aes(group = Group), level = 0.80, linewidth = 0.4, linetype = "dashed") +
  geom_point(size = 1.8, alpha = 0.9) +
  annotate("text",
    x = -Inf, y = Inf,
    label = sprintf("CAP %s\nLOO error = %.0f%%", fmt_p(cap_p), 100 * loo_err),
    hjust = -0.05, vjust = 1.15, size = FIG_GEOM_TEXT, color = "grey25", fontface = "bold"
  ) +
  scale_color_manual(values = GROUP_COLORS) +
  scale_fill_manual(values = GROUP_COLORS, guide = "none") +
  labs(
    title = "Supervised Ordination (CAP)",
    subtitle = "constrained on Group, conditioned on timepoint | permutation + leave-one-out",
    x = "CAP1 (responder axis)", y = "MDS1", tag = "N"
  ) +
  FIG_THEME +
  theme(legend.position = "bottom")

save_png(pN, file.path(supp_dir, "panel_n_cap"), PN_W, PN_H)
write.csv(sc, file.path(DAT_DIR, "audit_panel_N_cap.csv"), row.names = FALSE)
cat("F02 Supp Panel N (CAP) done.\n")
