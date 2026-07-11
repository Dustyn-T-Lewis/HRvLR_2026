# F02 Supp Panel M: variance explained by design, bias-corrected (omega^2) vs a null
# Raw eta^2 is upward-biased (median sits at the ~0.11 chance level for this 6-level,
# ~45-sample design), so most proteins read as "structured" when they are noise. This
# uses omega^2, which subtracts the null expectation and centres unstructured proteins
# near zero, overlaid on a label-permuted null so the real right tail is what exceeds
# chance (Lakens 2013; Olejnik & Algina 2003).

pacman::p_load(here, dplyr, tibble, tidyr, ggplot2)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

supp_dir <- file.path(RPT_DIR, "supp")
dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

PM_W <- 200
PM_H <- PANEL_H

x_mat <- imp_mat[, meta$Col_ID]
grp <- droplevels(meta$group)
k <- nlevels(grp)
n_samp <- ncol(x_mat)
grand <- rowMeans(x_mat)
ss_t <- rowSums((x_mat - grand)^2)

var_explained <- function(labels, corrected) {
  ind <- model.matrix(~ 0 + labels)
  n_g <- colSums(ind)
  means_g <- (x_mat %*% ind) / rep(n_g, each = nrow(x_mat))
  ss_b <- rowSums(sweep(means_g, 1, grand)^2 * rep(n_g, each = nrow(x_mat)))
  if (!corrected) {
    return(ss_b / ss_t)
  }
  ms_w <- (ss_t - ss_b) / (n_samp - k)
  (ss_b - (k - 1) * ms_w) / (ss_t + ms_w)
}

omega2 <- var_explained(grp, corrected = TRUE)

set.seed(42)
null_omega <- unlist(lapply(seq_len(100), function(i) {
  var_explained(factor(grp[sample.int(n_samp)]), corrected = TRUE)
}))
null_cut <- as.numeric(quantile(null_omega, 0.95, na.rm = TRUE))

eta_df <- tibble(uniprot_id = rownames(x_mat), omega2 = omega2) |>
  left_join(distinct(dep_df, uniprot_id, gene), by = "uniprot_id") |>
  mutate(label = if_else(is.na(gene) | gene == "", uniprot_id, gene))
n_hi <- sum(eta_df$omega2 > null_cut, na.rm = TRUE)
top_lab <- eta_df |> slice_max(omega2, n = 8, with_ties = FALSE)

null_df <- tibble(omega2 = null_omega)
bw <- 0.02

pM <- ggplot(eta_df, aes(omega2)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = bw, fill = "#6BAED6", color = "white", linewidth = 0.15) +
  geom_density(data = null_df, aes(omega2), color = "grey25", linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = null_cut, color = "#08519C", linewidth = 0.4) +
  annotate("text",
    x = null_cut, y = Inf, label = sprintf("%d proteins > null 95th (%.2f)", n_hi, null_cut),
    hjust = -0.05, vjust = 1.5, size = FIG_GEOM_TEXT, color = "#08519C", fontface = "bold"
  ) +
  ggrepel::geom_text_repel(
    data = top_lab, aes(x = omega2, y = 0, label = label),
    size = FIG_GEOM_TEXT, fontface = "bold", color = "#08519C", direction = "y",
    max.overlaps = 20, seed = 42, min.segment.length = 0, segment.size = 0.2
  ) +
  labs(
    title = expression(bold("Variance Explained by Group_Time (bias-corrected " * omega^2 * ")")),
    subtitle = "solid histogram = observed; dashed = label-permuted null | line = null 95th percentile",
    x = expression(omega^2), y = "density", tag = "M"
  ) +
  FIG_THEME

save_png(pM, file.path(supp_dir, "panel_m_omega2"), PM_W, PM_H)
write.csv(eta_df |> arrange(desc(omega2)),
  file.path(DAT_DIR, "audit_panel_M_omega2.csv"),
  row.names = FALSE
)
cat("F02 Supp Panel M done.\n")
