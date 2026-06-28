# F01 Panel I: mCSA (Pre/Post x Group + Delta)
pacman::p_load(here, dplyr, tidyr, readr, patchwork, ggsignif, rstatix)
source(here("04_Figures", "functions", "style.R"))

PW <- 150
PH <- 100

RPT <- here("04_Figures", "F01", "b_reports")
DAT <- here("04_Figures", "F01", "c_data")
dir.create(file.path(RPT, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

meta <- read_csv(here("00_input", "HRvLR_meta.csv"), show_col_types = FALSE) %>%
  filter(Timepoint %in% c("T1", "T2"), !is.na(mCSA_Pre)) %>%
  mutate(
    subject_key = sub("_T[123]$", "", Col_ID),
    Group = factor(Group, levels = c("HR", "LR")),
    Timepoint = factor(Timepoint, levels = c("T1", "T2")),
    Group_Time = factor(Group_Time, levels = c("HR_T1", "HR_T2", "LR_T1", "LR_T2"))
  ) %>%
  rename(value = mCSA_Pre)

pheno_wide <- meta %>%
  select(subject_key, Group, Timepoint, value) %>%
  pivot_wider(names_from = Timepoint, values_from = value) %>%
  mutate(delta = T2 - T1)

stats_anova <- rstatix::anova_test(
  data = meta, dv = value, wid = subject_key,
  between = Group, within = Timepoint
)

hr_wide <- pheno_wide %>% filter(Group == "HR")
lr_wide <- pheno_wide %>% filter(Group == "LR")
stats_paired_hr <- t.test(hr_wide$T2, hr_wide$T1, paired = TRUE)
stats_paired_lr <- t.test(lr_wide$T2, lr_wide$T1, paired = TRUE)
stats_delta <- t.test(delta ~ Group, data = pheno_wide)

anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf(
  "Group %s   Time %s   Interaction %s",
  fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
  fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
  fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"])
)

fill_4 <- GROUP_FILL[c("HR_T1", "HR_T2", "LR_T1", "LR_T2")]
sig_y_left <- bracket_pos(meta$value)

pI_left <- ggplot(meta, aes(x = Group_Time, y = value, fill = Group_Time)) +
  geom_bar(
    stat = "summary", fun = mean, width = 0.65,
    color = "black", linewidth = 0.3
  ) +
  geom_errorbar(
    stat = "summary", fun.data = mean_se,
    width = 0.2, linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.12, size = 1.2, alpha = 0.5,
    shape = 21, color = "black", stroke = 0.3
  ) +
  geom_signif(
    comparisons = list(c("HR_T1", "HR_T2")),
    annotations = fmt_p(stats_paired_hr$p.value),
    y_position = sig_y_left, textsize = KEY_TEXT, tip_length = 0.01
  ) +
  geom_signif(
    comparisons = list(c("LR_T1", "LR_T2")),
    annotations = fmt_p(stats_paired_lr$p.value),
    y_position = sig_y_left, textsize = KEY_TEXT, tip_length = 0.01
  ) +
  annotate("text",
    x = 1.5, y = -Inf, label = "HR",
    vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25"
  ) +
  annotate("text",
    x = 3.5, y = -Inf, label = "LR",
    vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25"
  ) +
  scale_fill_manual(values = fill_4) +
  scale_x_discrete(labels = c(
    HR_T1 = "Pre", HR_T2 = "Post",
    LR_T1 = "Pre", LR_T2 = "Post"
  )) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.22)),
    labels = scales::label_comma()
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "mCSA", subtitle = anova_sub,
    y = expression("mCSA (" * mm^2 * ")"), x = NULL, tag = "I"
  ) +
  FIG_THEME +
  theme(
    plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
    plot.margin = margin(5, 5, 20, 5), legend.position = "none"
  )

delta_colors <- c(
  HR = unname(GROUP_FILL["HR_T2"]),
  LR = unname(GROUP_FILL["LR_T2"])
)

pI_right <- ggplot(pheno_wide, aes(x = Group, y = delta, fill = Group)) +
  geom_bar(
    stat = "summary", fun = mean, width = 0.55,
    color = "black", linewidth = 0.3
  ) +
  geom_errorbar(
    stat = "summary", fun.data = mean_se,
    width = 0.15, linewidth = 0.4
  ) +
  geom_jitter(
    width = 0.12, size = 1.2, alpha = 0.5,
    shape = 21, color = "black", stroke = 0.3
  ) +
  geom_signif(
    comparisons = list(c("HR", "LR")),
    annotations = fmt_p(stats_delta$p.value),
    textsize = KEY_TEXT, tip_length = 0.02,
    y_position = bracket_pos(pheno_wide$delta)
  ) +
  scale_fill_manual(values = delta_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.25)),
    labels = scales::label_comma()
  ) +
  labs(y = expression(Delta ~ "mCSA (" * mm^2 * ")"), x = NULL) +
  FIG_THEME +
  theme(legend.position = "none")

pI <- (pI_left | pI_right) + plot_layout(widths = c(0.65, 0.35))

audit_I <- pheno_wide %>%
  group_by(Group) %>%
  summarise(
    n = sum(!is.na(delta)),
    delta_mean = mean(delta, na.rm = TRUE),
    delta_sd = sd(delta, na.rm = TRUE),
    delta_sem = delta_sd / sqrt(n),
    .groups = "drop"
  ) %>%
  mutate(delta_t_test_p = stats_delta$p.value)
write_csv(audit_I, file.path(DAT, "audit_panel_I_mcsa.csv"))

save_panel(pI, file.path(RPT, "panels", "panel_i_mcsa"), PW, PH)
