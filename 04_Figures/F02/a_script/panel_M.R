# F02 Panel M: variance structured by Group x Time (eta-squared) ----
# Top proteins whose abundance variance is explained by the responder x timepoint
# design - the proteins carrying the HR-vs-LR signal. eta2 from Stage 02 intermediates.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

PM_W <- 150
PM_H <- 120
N_TOP <- 20

int <- readRDS("02_Normalization/c_data/00_report_intermediates.rds")
eta_df <- tibble(uniprot_id = names(int$eta2_vals), eta2 = as.numeric(int$eta2_vals)) |>
  left_join(distinct(dep_df, uniprot_id, gene), by = "uniprot_id") |>
  mutate(label = if_else(is.na(gene) | gene == "", uniprot_id, gene))

med_eta <- median(eta_df$eta2, na.rm = TRUE)
n_high <- sum(eta_df$eta2 > 0.5, na.rm = TRUE)

top_df <- eta_df |>
  slice_max(eta2, n = N_TOP, with_ties = FALSE) |>
  mutate(label = factor(label, levels = rev(label)))

pM <- ggplot(top_df, aes(eta2, label)) +
  geom_segment(aes(x = 0, xend = eta2, yend = label), color = "grey75", linewidth = 0.5) +
  geom_point(aes(fill = eta2), shape = 21, size = 3, color = "grey30") +
  scale_fill_gradient(low = "#C6DBEF", high = "#08519C", guide = "none") +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = expression(bold("Variance Structured by Group " %*% " Time (" * eta^2 * ")")),
    subtitle = sprintf(
      "Top %d proteins | median eta² = %.3f | %d proteins > 0.5",
      N_TOP, med_eta, n_high
    ),
    x = expression(eta^2 ~ "(Group_Time-explained variance)"), y = NULL, tag = "M"
  ) +
  FIG_THEME +
  theme(axis.text.y = element_text(size = 7, face = "bold"))

save_panel(pM, file.path(RPT_DIR, "panel_m_eta2"), PM_W, PM_H)
write.csv(eta_df |> arrange(desc(eta2)),
  file.path(DAT_DIR, "audit_panel_M_eta2.csv"),
  row.names = FALSE
)
cat("F02 Panel M done.\n")
