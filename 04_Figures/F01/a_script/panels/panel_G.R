# F01 Panel G: Hypertrophy Composite
pacman::p_load(here, dplyr, readr)
source(here("04_Figures", "functions", "style.R"))

PW <- 150
PH <- 100

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE) |>
  filter(Timepoint == "T2", !is.na(COMP.HYPERTROPHY)) |>
  transmute(
    Subject_ID,
    Group = factor(Group, levels = c("HR", "LR")),
    hypertrophy_pct = readr::parse_number(COMP.HYPERTROPHY)
  )

bar_colors <- c(
  HR = unname(GROUP_FILL["HR_T2"]),
  LR = unname(GROUP_FILL["LR_T2"])
)

pG <- ggplot(meta, aes(x = Group, y = hypertrophy_pct, fill = Group)) +
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
  scale_fill_manual(values = bar_colors) +
  scale_x_discrete(labels = c(
    HR = sprintf("HR (n = %d)", sum(meta$Group == "HR")),
    LR = sprintf("LR (n = %d)", sum(meta$Group == "LR"))
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15)),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    title = "Hypertrophy Composite", subtitle = "Groups defined by this composite",
    y = "Composite hypertrophy (%)", x = NULL, tag = "G"
  ) +
  coord_cartesian(clip = "off") +
  FIG_THEME +
  theme(
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

audit_G <- meta |>
  group_by(Group) |>
  summarise(
    n = n(),
    mean = mean(hypertrophy_pct),
    sd = sd(hypertrophy_pct),
    sem = sd(hypertrophy_pct) / sqrt(n()),
    .groups = "drop"
  )
write_csv(audit_G, file.path(DAT, "audit_panel_G_hypertrophy_composite.csv"))

save_panel(pG, file.path(RPT, "panels", "panel_g_hypertrophy_composite"), PW, PH)
