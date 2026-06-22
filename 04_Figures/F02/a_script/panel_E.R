# F02 Panel E: Intra-Individual Proteomic Variability (Imputed) ----
# One boxplot per subject, faceted by HR/LR, ordered by median log2FC.
# Annotated with per-subject imputation fractions.
# Adapted from YvO F02 panel_E pattern.

setwd(rprojroot::find_rstudio_root_file())
if (!exists("meta")) source("04_Figures/F02/a_script/HRvLR_F02_setup.R")

library(dplyr)
library(tidyr)
library(ggplot2)

PE_W <- 160
PE_H <- 90
dir.create(RPT_DIR, recursive = TRUE, showWarnings = FALSE)

n_proteins <- nrow(imp_mat)

# Imputation mask (TRUE = was missing, i.e. imputed) derived from the
# normalized matrix: NAs there are exactly the values missForest filled in.
norm_dal <- readRDS("02_Normalization/c_data/DAList_normalized.rds")
mask_src <- is.na(as.matrix(norm_dal$data)) # rownames = uniprot_id
uid2gene <- setNames(imp_dal$annotation$gene, imp_dal$annotation$uniprot_id)
mask_df <- as_tibble(mask_src) |>
  mutate(gene = unname(uid2gene[rownames(mask_src)]), .before = 1)
has_mask <- TRUE
mask_mat <- as.matrix(mask_df[, intersect(imp_samps, names(mask_df))])
rownames(mask_mat) <- mask_df$gene

# MAR/MNAR classification no longer produced upstream; MNAR annotation off.
has_mnar <- FALSE

# Compute per-subject log2FC (T2 - T1) for training effect
subjects <- unique(meta$Subject_ID)

lfc_list <- lapply(subjects, function(s) {
  t1_id <- meta$Col_ID[meta$Subject_ID == s & meta$Timepoint == "T1"]
  t2_id <- meta$Col_ID[meta$Subject_ID == s & meta$Timepoint == "T2"]
  if (length(t1_id) != 1 || length(t2_id) != 1) {
    return(NULL)
  }

  lfc <- imp_mat[, t2_id] - imp_mat[, t1_id] # log2(T2) - log2(T1)

  # Count logFC values depending on at least one imputed value
  n_imputed_lfc <- 0
  pct_imputed <- 0
  n_mnar_lfc <- 0
  pct_mnar <- 0

  if (has_mask) {
    t1_imp <- if (t1_id %in% colnames(mask_mat)) mask_mat[, t1_id] else rep(FALSE, n_proteins)
    t2_imp <- if (t2_id %in% colnames(mask_mat)) mask_mat[, t2_id] else rep(FALSE, n_proteins)
    either_imputed <- t1_imp | t2_imp
    n_imputed_lfc <- sum(either_imputed)
    pct_imputed <- 100 * n_imputed_lfc / n_proteins

    if (has_mnar) {
      genes <- imp_df$gene
      mnar_and_imputed <- either_imputed & (genes %in% mnar_genes)
      n_mnar_lfc <- sum(mnar_and_imputed)
      pct_mnar <- 100 * n_mnar_lfc / n_proteins
    }
  }

  tibble(
    subject       = s,
    resp_group    = as.character(meta$Group[meta$Subject_ID == s][1]),
    lfc           = as.numeric(lfc),
    n_imputed_lfc = n_imputed_lfc,
    pct_imputed   = pct_imputed,
    n_mnar_lfc    = n_mnar_lfc,
    pct_mnar      = pct_mnar
  )
})

lfc_long <- bind_rows(lfc_list)
lfc_long$resp_group <- factor(lfc_long$resp_group, levels = c("HR", "LR"))

subj_summary <- lfc_long |>
  group_by(
    subject, resp_group, n_imputed_lfc, pct_imputed,
    n_mnar_lfc, pct_mnar
  ) |>
  summarise(
    median_lfc = median(lfc, na.rm = TRUE),
    mad_lfc    = mad(lfc, na.rm = TRUE),
    sd_lfc     = sd(lfc, na.rm = TRUE),
    iqr_lfc    = IQR(lfc, na.rm = TRUE),
    q25        = quantile(lfc, 0.25, na.rm = TRUE),
    q75        = quantile(lfc, 0.75, na.rm = TRUE),
    n_proteins = n(),
    .groups    = "drop"
  )

subj_summary <- subj_summary |>
  arrange(resp_group, median_lfc) |>
  mutate(subj_order = factor(subject, levels = unique(subject)))

lfc_long <- lfc_long |>
  mutate(subj_order = factor(subject, levels = levels(subj_summary$subj_order)))

group_summary <- subj_summary |>
  group_by(resp_group) |>
  summarise(
    mean_median = mean(median_lfc),
    sd_median   = sd(median_lfc),
    n           = n(),
    .groups     = "drop"
  )

wt <- wilcox.test(median_lfc ~ resp_group, data = subj_summary)
n1 <- sum(subj_summary$resp_group == "HR")
n2 <- sum(subj_summary$resp_group == "LR")
r_rb <- 1 - 2 * wt$statistic / (n1 * n2)

mean_pct_imp <- mean(subj_summary$pct_imputed)
mean_pct_mnar <- mean(subj_summary$pct_mnar)
subtitle_text <- if (has_mnar) {
  sprintf(
    "%s proteins (imputed) | %.0f%% imputed, %.0f%% MNAR\nWilcoxon %s",
    format(n_proteins, big.mark = ","), mean_pct_imp, mean_pct_mnar, fmt_p(wt$p.value)
  )
} else {
  sprintf(
    "%s proteins (imputed) | %.0f%% imputed\nWilcoxon %s",
    format(n_proteins, big.mark = ","), mean_pct_imp, fmt_p(wt$p.value)
  )
}

pE <- ggplot(lfc_long, aes(x = subj_order, y = lfc, fill = resp_group)) +
  geom_boxplot(
    width = 0.5, linewidth = 0.3, color = "black",
    outlier.shape = NA, alpha = 0.5
  ) +
  facet_grid(~resp_group, scales = "free_x", space = "free_x") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(
    x = "Subject",
    y = expression(bold(Delta ~ log[2] * "FC (T2/T1)")),
    title = "Intra-Individual Proteomic Variability",
    subtitle = subtitle_text,
    tag = "E"
  ) +
  FIG_THEME +
  theme(
    legend.position = "none",
    panel.spacing = unit(3, "mm"),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5,
      size = FIG_AXIS_TEXT - 1.5
    )
  )

write.csv(
  subj_summary |>
    select(
      subject, resp_group, median_lfc, mad_lfc, sd_lfc,
      iqr_lfc, q25, q75, n_proteins,
      n_imputed_lfc, pct_imputed, n_mnar_lfc, pct_mnar
    ),
  file.path(DAT_DIR, "audit_panel_E_imputed.csv"),
  row.names = FALSE
)

write.csv(group_summary,
  file.path(DAT_DIR, "audit_panel_E_wilcoxon.csv"),
  row.names = FALSE
)

save_panel(pE, file.path(RPT_DIR, "panel_E_imputed_SUPP"), PE_W, PE_H)
cat("F02 Panel E done.\n")
