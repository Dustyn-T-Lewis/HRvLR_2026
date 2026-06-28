# F02 Panel M: proteome-wide variance structured by Group x Time (eta-squared)
# Distribution of eta2 across all proteins: how much of each protein's variance the
# responder x timepoint design explains. The right tail = design-structured proteins.

pacman::p_load(here, dplyr, tibble, ggplot2, ggrepel)

if (!exists("meta")) source(here("04_Figures", "F02", "a_script", "HRvLR_F02_setup.R"))

PM_W <- 200
PM_H <- PANEL_H

int <- readRDS(here("02_Normalization", "c_data", "00_report_intermediates.rds"))
eta_df <- tibble(uniprot_id = names(int$eta2_vals), eta2 = as.numeric(int$eta2_vals)) |>
  left_join(distinct(dep_df, uniprot_id, gene), by = "uniprot_id") |>
  mutate(label = if_else(is.na(gene) | gene == "", uniprot_id, gene))

med_eta <- median(eta_df$eta2, na.rm = TRUE)
hi_cut <- 0.5
n_hi <- sum(eta_df$eta2 > hi_cut)
top_lab <- eta_df |> slice_max(eta2, n = 8, with_ties = FALSE)

pM <- ggplot(eta_df, aes(eta2)) +
  geom_histogram(aes(y = after_stat(count)),
    bins = 60,
    fill = "#6BAED6", color = "white", linewidth = 0.15
  ) +
  geom_vline(xintercept = med_eta, linetype = "dashed", color = "grey30", linewidth = 0.4) +
  annotate("rect",
    xmin = hi_cut, xmax = Inf, ymin = 0, ymax = Inf,
    fill = "#08519C", alpha = 0.06
  ) +
  geom_rug(data = top_lab, aes(eta2), color = "#08519C", linewidth = 0.4) +
  ggrepel::geom_text_repel(
    data = top_lab, aes(x = eta2, y = 0, label = label),
    size = FIG_GEOM_TEXT, fontface = "bold", color = "#08519C", direction = "y",
    nudge_y = max(table(cut(eta_df$eta2, 60))) * 0.4, segment.size = 0.2,
    max.overlaps = 20, seed = 42
  ) +
  annotate("text",
    x = med_eta, y = Inf, label = sprintf("median = %.3f", med_eta),
    hjust = -0.1, vjust = 1.5, size = FIG_GEOM_TEXT, color = "grey30", fontface = "bold"
  ) +
  annotate("text",
    x = hi_cut, y = Inf, label = sprintf("%d proteins > %.1f", n_hi, hi_cut),
    hjust = -0.1, vjust = 3, size = FIG_GEOM_TEXT, color = "#08519C", fontface = "bold"
  ) +
  scale_x_continuous(limits = c(0, max(eta_df$eta2) * 1.05), expand = expansion(mult = c(0, 0))) +
  labs(
    title = expression(bold("Proteome-Wide Variance Structured by Group " %*% " Time (" * eta^2 * ")")),
    subtitle = sprintf(
      "%s proteins | fraction of each protein's variance explained by Group_Time",
      format(nrow(eta_df), big.mark = ",")
    ),
    x = expression(eta^2), y = "Proteins", tag = "M"
  ) +
  FIG_THEME

save_png(pM, file.path(RPT_DIR, "panels", "panel_m_eta2"), PM_W, PM_H)
write.csv(eta_df |> arrange(desc(eta2)),
  file.path(DAT_DIR, "audit_panel_M_eta2.csv"),
  row.names = FALSE
)
cat("F02 Panel M done.\n")
