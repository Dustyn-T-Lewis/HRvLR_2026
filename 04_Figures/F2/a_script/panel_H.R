# panel_H.R — Training Scatter: Training_HR vs Training_LR logFC
# Protein-level concordance of training response between responder groups.

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(boot)
})

here::i_am(".here")
source(here::here("04_Figures/F2/a_script/style_f2.R"))

PW  <- 200
RPT <- here::here("04_Figures/F2/b_reports")
DAT <- here::here("04_Figures/F2/c_data")
dir.create(file.path(DAT, "panel_H"), recursive = TRUE, showWarnings = FALSE)

dep_df <- read_csv(here::here("03_DEP/c_data/03_combined_results.csv"),
                   show_col_types = FALSE)

scatter_df <- dep_df |>
  transmute(
    gene,
    logFC_HR       = logFC_Training_HR,
    logFC_LR       = logFC_Training_LR,
    pi_HR          = pi_score_Training_HR,
    pi_LR          = pi_score_Training_LR,
    pi_Interaction = pi_score_Training_Interaction
  ) |>
  filter(!is.na(logFC_HR), !is.na(logFC_LR)) |>
  mutate(
    significance = classify_proteins_scatter(pi_HR, pi_LR, pi_Interaction),
    quadrant = case_when(
      logFC_HR > 0 & logFC_LR > 0 ~ "Concordant Up",
      logFC_HR < 0 & logFC_LR < 0 ~ "Concordant Down",
      logFC_HR > 0 & logFC_LR < 0 ~ "HR\u2191 LR\u2193",
      TRUE                         ~ "HR\u2193 LR\u2191"
    ),
    point_size   = ifelse(significance == "NS", 1.8, 2.3),
    bubble_alpha = case_when(
      significance == "NS"          ~ 0.30,
      significance == "Interaction" ~ 0.55,
      significance == "Sig Both"    ~ 0.75,
      TRUE                          ~ 0.85
    )
  )

# ── Statistics ────────────────────────────────────────────────────────────────
n_obs   <- nrow(scatter_df)
cor_r   <- cor.test(scatter_df$logFC_HR, scatter_df$logFC_LR, method = "pearson")
cor_rho <- cor.test(scatter_df$logFC_HR, scatter_df$logFC_LR, method = "spearman")
rho_ci  <- tanh(atanh(cor_rho$estimate) + c(-1, 1) * qnorm(0.975) / sqrt(n_obs - 3))

conc_set <- scatter_df |> filter(abs(logFC_HR) > 0.2 | abs(logFC_LR) > 0.2)
sign_conc <- mean(sign(conc_set$logFC_HR) == sign(conc_set$logFC_LR)) * 100

set.seed(42)
boot_sc <- boot::boot(
  data = conc_set,
  statistic = function(d, i)
    mean(sign(d$logFC_HR[i]) == sign(d$logFC_LR[i])) * 100,
  R = 10000)
boot_ci <- tryCatch(
  boot::boot.ci(boot_sc, type = "bca", conf = 0.95)$bca[4:5],
  error = function(e) quantile(boot_sc$t, c(0.025, 0.975)))

stats_df <- tibble(
  metric   = c("Pearson_r", "Spearman_rho", "Sign_concordance_pct"),
  estimate = c(cor_r$estimate, cor_rho$estimate, sign_conc),
  ci_lower = c(cor_r$conf.int[1], rho_ci[1], boot_ci[1]),
  ci_upper = c(cor_r$conf.int[2], rho_ci[2], boot_ci[2]),
  p_value  = c(cor_r$p.value, cor_rho$p.value, NA_real_),
  n        = c(n_obs, n_obs, nrow(conc_set)),
  note     = c("95% CI from cor.test()",
               "95% CI via Fisher z-transformation",
               "95% BCa bootstrap CI (10000 replicates, |logFC|>0.2 filter)")
)
write_csv(stats_df, file.path(DAT, "panel_H", "training_scatter_stats.csv"))

# ── Quadrant counts ───────────────────────────────────────────────────────────
axis_max <- max(abs(c(scatter_df$logFC_HR, scatter_df$logFC_LR)), na.rm = TRUE) * 1.05
xlim_rng <- c(-axis_max, axis_max)
ylim_rng <- c(-axis_max, axis_max)

q_df     <- scatter_df |> mutate(q = case_when(
  logFC_HR > 0 & logFC_LR > 0 ~ "Q1",
  logFC_HR < 0 & logFC_LR < 0 ~ "Q3",
  logFC_HR > 0 & logFC_LR < 0 ~ "Q4",
  TRUE ~ "Q2"))
q_counts <- q_df |> count(q) |> deframe()
q_sig    <- q_df |> filter(significance != "NS") |> count(q) |> deframe()
for (qq in c("Q1","Q2","Q3","Q4")) if (is.na(q_sig[qq])) q_sig[qq] <- 0

# ── Label top proteins ────────────────────────────────────────────────────────
txt_gene <- scale_text(BASE_GENE, PW)
txt_quad <- scale_text(BASE_QUADRANT, PW)

label_df <- scatter_df |>
  filter(significance != "NS") |>
  group_by(significance) |>
  arrange(desc(abs(logFC_HR) + abs(logFC_LR))) |>
  slice_head(n = 5) |>
  ungroup() |>
  mutate(
    label_fill = SIG_LABEL_FILL_SCATTER[as.character(significance)],
    label_col  = SIG_LABEL_TEXT_SCATTER[as.character(significance)]
  )

ns_df  <- scatter_df |> filter(significance == "NS")
sig_df <- scatter_df |> filter(significance != "NS")

# ── Plot ──────────────────────────────────────────────────────────────────────
pH <- ggplot(mapping = aes(x = logFC_HR, y = logFC_LR)) +
  annotate("rect", xmin = 0,   xmax =  Inf, ymin = 0,    ymax =  Inf, fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0,   ymin = -Inf, ymax =  0,   fill = "#FFE0E0", alpha = 0.55) +
  annotate("rect", xmin = 0,   xmax =  Inf, ymin = -Inf, ymax =  0,   fill = "#DCEEFF", alpha = 0.55) +
  annotate("rect", xmin = -Inf, xmax = 0,   ymin = 0,    ymax =  Inf, fill = "#DCEEFF", alpha = 0.55) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", linewidth = 0.3) +
  geom_point(data = ns_df,  aes(x = logFC_HR, y = logFC_LR),
             color = "grey80", fill = "grey85", shape = 21,
             size = 1.0, alpha = 0.3, stroke = 0.2) +
  geom_point(data = sig_df, aes(x = logFC_HR, y = logFC_LR, fill = significance),
             shape = 21, size = sig_df$point_size,
             color = "grey75", alpha = sig_df$bubble_alpha, stroke = 0.8) +
  scale_fill_manual(values = SIG_COLORS_SCATTER, name = "Significance") +
  geom_label_repel(
    data = label_df, aes(x = logFC_HR, y = logFC_LR, label = gene),
    fill = label_df$label_fill, color = label_df$label_col,
    size = txt_gene, fontface = "italic",
    max.overlaps = 40, segment.size = 0.2, segment.color = "grey50",
    min.segment.length = 0, show.legend = FALSE,
    box.padding = 0.6, force = 3, force_pull = 0.5,
    label.padding = unit(1.5, "pt"), label.r = unit(1, "pt"),
    label.size = 0.15, seed = 42) +
  annotate("label", x = Inf,  y = Inf,  hjust = 1, vjust = 1,
           label = sprintf("Concordant Up\nn = %s / %s", q_sig["Q1"], q_counts["Q1"]),
           size = txt_quad, fontface = "bold", color = "#DC2626",
           fill = scales::alpha("white", 0.9), label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = -Inf, hjust = 0, vjust = 0,
           label = sprintf("Concordant Down\nn = %s / %s", q_sig["Q3"], q_counts["Q3"]),
           size = txt_quad, fontface = "bold", color = "#DC2626",
           fill = scales::alpha("white", 0.9), label.padding = unit(2.5, "pt")) +
  annotate("label", x = Inf,  y = -Inf, hjust = 1, vjust = 0,
           label = sprintf("HR\u2191 LR\u2193\nn = %s / %s", q_sig["Q4"], q_counts["Q4"]),
           size = txt_quad, fontface = "bold", color = "#2563EB",
           fill = scales::alpha("white", 0.9), label.padding = unit(2.5, "pt")) +
  annotate("label", x = -Inf, y = Inf,  hjust = 0, vjust = 1,
           label = sprintf("HR\u2193 LR\u2191\nn = %s / %s", q_sig["Q2"], q_counts["Q2"]),
           size = txt_quad, fontface = "bold", color = "#2563EB",
           fill = scales::alpha("white", 0.9), label.padding = unit(2.5, "pt")) +
  coord_fixed(ratio = 1, xlim = xlim_rng, ylim = ylim_rng, expand = FALSE) +
  labs(
    title    = "Training Response: HR vs LR — Protein-Level Concordance",
    subtitle = sprintf(
      "n = %s | r = %.3f [%.3f, %.3f] | \u03c1 = %.3f [%.3f, %.3f] | concordance = %.0f%% [%.0f, %.0f]",
      format(n_obs, big.mark = ","),
      cor_r$estimate, cor_r$conf.int[1], cor_r$conf.int[2],
      cor_rho$estimate, rho_ci[1], rho_ci[2],
      sign_conc, boot_ci[1], boot_ci[2]),
    x = expression(log[2]*FC ~ "(Training HR: T2 \u2212 T1)"),
    y = expression(log[2]*FC ~ "(Training LR: T2 \u2212 T1)")
  ) +
  FIG_THEME +
  theme(legend.position = "none")

save_panel(pH, file.path(RPT, "panel_H_training_scatter"), PW, PW)

scatter_df |>
  transmute(gene,
            logFC_Training_HR    = round(logFC_HR, 4),
            logFC_Training_LR    = round(logFC_LR, 4),
            pi_HR                = round(pi_HR, 6),
            pi_LR                = round(pi_LR, 6),
            pi_Interaction       = round(pi_Interaction, 6),
            significance         = as.character(significance),
            quadrant) |>
  arrange(significance, desc(abs(logFC_Training_HR) + abs(logFC_Training_LR))) |>
  write_csv(file.path(DAT, "panel_H", "training_scatter.csv"))

message("Panel H done: Training scatter")
