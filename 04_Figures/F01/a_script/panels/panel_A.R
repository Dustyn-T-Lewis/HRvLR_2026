# F01 Panel A: Accumulated Volume Load
pacman::p_load(here, dplyr, readr, ggsignif)
source(here("04_Figures", "functions", "style.R"))

PW <- 150
PH <- 100

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE)

vl_df <- meta %>%
  filter(Timepoint == "T3", !is.na(ACCUM_VL)) %>%
  select(Subject_ID, Group, ACCUM_VL) %>%
  mutate(Group = factor(Group, levels = c("HR", "LR")))

stats_A <- t.test(ACCUM_VL ~ Group, data = vl_df)

bar_colors <- c(
  HR = unname(GROUP_FILL["HR_T2"]),
  LR = unname(GROUP_FILL["LR_T2"])
)

pA <- ggplot(vl_df, aes(x = Group, y = ACCUM_VL, fill = Group)) +
  geom_bar(
    stat = "summary", fun = mean, width = 0.45,
    color = "black", linewidth = 0.3
  ) +
  geom_errorbar(
    stat = "summary", fun.data = mean_se,
    width = 0.2, linewidth = 0.4
  ) +
  geom_point(
    position = position_jitter(width = 0.15, height = 0, seed = 42),
    size = 1.5, alpha = 0.5, shape = 16, color = "grey30"
  ) +
  geom_signif(
    comparisons = list(c("HR", "LR")),
    annotations = fmt_p(stats_A$p.value),
    textsize = KEY_TEXT, tip_length = 0.02,
    y_position = bracket_pos(vl_df$ACCUM_VL)
  ) +
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(
    HR = sprintf("HR (n = %d)", sum(vl_df$Group == "HR")),
    LR = sprintf("LR (n = %d)", sum(vl_df$Group == "LR"))
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = scales::label_comma()
  ) +
  labs(
    title = "Accumulated Volume Load", y = "Volume load (kg)", x = NULL,
    tag = "A"
  ) +
  coord_cartesian(clip = "off") +
  FIG_THEME +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

audit_A <- vl_df %>%
  group_by(Group) %>%
  summarise(
    n = n(), mean = mean(ACCUM_VL), sd = sd(ACCUM_VL),
    sem = sd(ACCUM_VL) / sqrt(n()), .groups = "drop"
  ) %>%
  mutate(t_test_p = stats_A$p.value)
write_csv(audit_A, file.path(DAT, "audit_panel_A_volume_load.csv"))

save_panel(pA, file.path(RPT, "panels", "panel_a_volume_load"), PW, PH)
