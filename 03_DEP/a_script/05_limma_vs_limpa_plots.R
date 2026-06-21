# HRvLR — Concordance plots: limma+missForest vs LIMPA

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

setwd(rprojroot::find_rstudio_root_file())

cfg <- list(
  detail_csv  = "03_DEP/c_data/limpa_compare/03_per_protein_detail.csv",
  summary_csv = "03_DEP/c_data/limpa_compare/01_compare_summary.csv",
  out_dir     = "03_DEP/b_reports/limpa_compare",
  fdr_cut     = 0.10
)
dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)

detail  <- read_csv(cfg$detail_csv,  show_col_types = FALSE)
summary <- read_csv(cfg$summary_csv, show_col_types = FALSE)

detail <- detail |>
  mutate(
    cat = case_when(
      FDR_lm <  cfg$fdr_cut & FDR_lp <  cfg$fdr_cut ~ "both",
      FDR_lm <  cfg$fdr_cut & FDR_lp >= cfg$fdr_cut ~ "limma only",
      FDR_lm >= cfg$fdr_cut & FDR_lp <  cfg$fdr_cut ~ "limpa only",
      TRUE ~ "neither"
    ),
    sign_disagree = sign(logFC_lm) != sign(logFC_lp)
  )

labels <- summary |>
  transmute(contrast,
            lab = sprintf("ρ=%.2f  fdr10: lm=%d / lp=%d / ∩=%d",
                          spearman_logFC, fdr10_lm, fdr10_lp, fdr10_overlap))

cat_palette <- c(both = "#2E7D32", `limma only` = "#1976D2",
                 `limpa only` = "#E65100", neither = "grey75")
cat_alpha   <- c(both = 0.95, `limma only` = 0.85,
                 `limpa only` = 0.85, neither = 0.25)

contrast_levels <- c("Training_HR", "Training_LR",
                     "Acute_HR", "Acute_LR",
                     "Baseline_HRvLR", "Trained_HRvLR", "Acute_HRvLR",
                     "Training_Interaction", "Acute_Interaction")
detail$contrast <- factor(detail$contrast, levels = contrast_levels)
labels$contrast <- factor(labels$contrast, levels = contrast_levels)

p_lfc <- detail |>
  ggplot(aes(logFC_lm, logFC_lp, colour = cat, alpha = cat)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.2) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.2) +
  geom_point(size = 0.55) +
  geom_point(data = ~ filter(.x, sign_disagree),
             colour = "red", shape = 4, size = 1, alpha = 0.6) +
  geom_text(data = labels, aes(x = -Inf, y = Inf, label = lab),
            inherit.aes = FALSE, hjust = -0.05, vjust = 1.4, size = 2.6) +
  scale_colour_manual(values = cat_palette, name = paste0("FDR<", cfg$fdr_cut)) +
  scale_alpha_manual(values = cat_alpha, guide = "none") +
  facet_wrap(~ contrast, scales = "free", ncol = 3) +
  labs(x = "logFC (limma + missForest)", y = "logFC (LIMPA)",
       title = "HRvLR — logFC concordance: limma vs LIMPA",
       subtitle = "Diagonal = perfect agreement. Red ✕ = sign disagreement.") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"))

ggsave(file.path(cfg$out_dir, "logFC_scatter.pdf"),
       p_lfc, width = 11, height = 11)

p_p <- detail |>
  mutate(mlp_lm = -log10(P_lm), mlp_lp = -log10(P_lp)) |>
  ggplot(aes(mlp_lm, mlp_lp, colour = cat, alpha = cat)) +
  geom_abline(slope = 1, intercept = 0, colour = "black", linewidth = 0.4) +
  geom_point(size = 0.55) +
  scale_colour_manual(values = cat_palette, name = paste0("FDR<", cfg$fdr_cut)) +
  scale_alpha_manual(values = cat_alpha, guide = "none") +
  facet_wrap(~ contrast, scales = "free", ncol = 3) +
  labs(x = expression(-log[10]~P~(limma)), y = expression(-log[10]~P~(LIMPA)),
       title = "HRvLR — significance concordance: limma vs LIMPA") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"))

ggsave(file.path(cfg$out_dir, "neglog10P_scatter.pdf"),
       p_p, width = 11, height = 11)

cat(sprintf("Saved -> %s\n", cfg$out_dir))
