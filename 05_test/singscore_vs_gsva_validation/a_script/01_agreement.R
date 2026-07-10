# Agreement between the single-sample singscore and the cohort-relative GSVA:
# a per-pathway Spearman across samples, the pooled Spearman over every
# pathway-sample pair, and a rank scatter. Where they disagree is where the
# cohort-relative rescaling of GSVA moves a sample against its own rank profile.
if (!exists("sing")) source(here::here("05_test", "singscore_vs_gsva_validation", "a_script", "setup.R"))
pacman::p_load(purrr)

per_pathway <- tibble(
  pathway = rownames(sing),
  rho = vapply(rownames(sing), function(p) {
    suppressWarnings(cor(sing[p, ], gsva[p, ], method = "spearman"))
  }, numeric(1))
) |>
  mutate(database = classify_database(pathway)) |>
  arrange(rho)

pooled_rho <- suppressWarnings(cor(
  as.vector(sing), as.vector(gsva),
  method = "spearman"
))

write_csv(
  per_pathway |> mutate(pooled_rho = pooled_rho),
  file.path(dat, "score_agreement.csv")
)

long <- tibble(
  pathway = rep(rownames(sing), times = ncol(sing)),
  sample = rep(colnames(sing), each = nrow(sing)),
  singscore = as.vector(sing), gsva = as.vector(gsva)
) |>
  group_by(pathway) |>
  mutate(sing_rank = rank(singscore), gsva_rank = rank(gsva)) |>
  ungroup()

med_rho <- median(per_pathway$rho)

fig_scatter <- ggplot(long, aes(sing_rank, gsva_rank)) +
  geom_point(size = 0.7, alpha = 0.25, color = unname(GROUP_COLORS["HR"])) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey55") +
  labs(
    x = "singscore rank (within pathway)", y = "GSVA rank (within pathway)",
    title = "singscore and GSVA order samples almost identically",
    subtitle = sprintf(
      "Within pathways they agree (median rho = %.2f over %d pathways); across pathways the absolute scales differ (pooled rho = %.2f)",
      med_rho, nrow(sing), pooled_rho
    )
  ) +
  FIG_THEME
save_panel(fig_scatter, file.path(rpt, "fig1_rank_scatter"), 150, 130)

fig_dist <- ggplot(per_pathway, aes(rho)) +
  geom_histogram(binwidth = 0.02, fill = "#41AB5D", color = "white") +
  geom_vline(xintercept = med_rho, linetype = "dashed", color = "grey30") +
  annotate(
    "text",
    x = med_rho, y = Inf, vjust = 1.5, hjust = 1.05,
    label = sprintf("median rho = %.3f", med_rho), size = 3
  ) +
  labs(
    x = "per-pathway Spearman rho (across samples)", y = "pathways",
    title = "Per-pathway agreement across samples",
    subtitle = sprintf(
      "%d of %d pathways at rho > 0.8; lowest = %s (rho = %.2f)",
      sum(per_pathway$rho > 0.8), nrow(per_pathway),
      clean_pathway_name(per_pathway$pathway[1], 26), per_pathway$rho[1]
    )
  ) +
  FIG_THEME
save_panel(fig_dist, file.path(rpt, "fig2_rho_distribution"), 160, 120)

message(sprintf(
  "agreement: pooled rho %.3f, median per-pathway rho %.3f, min %.3f (%s)",
  pooled_rho, med_rho, per_pathway$rho[1], per_pathway$pathway[1]
))
