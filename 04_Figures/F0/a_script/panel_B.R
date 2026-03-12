# F0 Panel B: 1RM Leg Press (Pre/Post x Group + Delta) ----
setwd(rprojroot::find_rstudio_root_file())
source("04_Figures/shared/style.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(ggsignif)
  library(rstatix)
})

PW <- 170
PH <- 80

RPT <- "04_Figures/F0/b_reports"
DAT <- "04_Figures/F0/c_data"
dir.create(RPT, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

# Data ----
meta <- read_csv("00_input/HRvLR_meta.csv", show_col_types = FALSE) %>%
  filter(Timepoint %in% c("T1", "T2")) %>%
  mutate(
    subject_key = sub("_T[123]$", "", Col_ID),
    Group       = factor(Group, levels = c("HR", "LR")),
    Timepoint   = factor(Timepoint, levels = c("T1", "T2")),
    Group_Time  = factor(Group_Time,
                         levels = c("HR_T1", "HR_T2", "LR_T1", "LR_T2"))
  ) %>%
  rename(value = X1RM_Leg_Pre)

pheno_wide <- meta %>%
  select(subject_key, Group, Timepoint, value) %>%
  pivot_wider(names_from = Timepoint, values_from = value) %>%
  mutate(delta = T2 - T1)

stats_anova <- rstatix::anova_test(data = meta, dv = value,
                                    wid = subject_key,
                                    between = Group, within = Timepoint)

hr_wide <- pheno_wide %>% filter(Group == "HR")
lr_wide <- pheno_wide %>% filter(Group == "LR")
stats_paired_hr <- t.test(hr_wide$T2, hr_wide$T1, paired = TRUE)
stats_paired_lr <- t.test(lr_wide$T2, lr_wide$T1, paired = TRUE)

# Plot ----
anova_tbl <- as.data.frame(stats_anova)
anova_sub <- sprintf("Group %s   Time %s   Interaction %s",
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                     fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

fill_4 <- GROUP_FILL[c("HR_T1", "HR_T2", "LR_T1", "LR_T2")]
sig_y_left <- bracket_pos(meta$value)

pB_left <- ggplot(meta, aes(x = Group_Time, y = value, fill = Group_Time)) +
  geom_bar(stat = "summary", fun = mean, width = 0.65,
           color = "black", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.2, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  geom_signif(comparisons = list(c("HR_T1", "HR_T2")),
              annotations = fmt_p(stats_paired_hr$p.value),
              y_position = sig_y_left, textsize = KEY_TEXT, tip_length = 0.01) +
  geom_signif(comparisons = list(c("LR_T1", "LR_T2")),
              annotations = fmt_p(stats_paired_lr$p.value),
              y_position = sig_y_left, textsize = KEY_TEXT, tip_length = 0.01) +
  annotate("text", x = 1.5, y = -Inf, label = "HR",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  annotate("text", x = 3.5, y = -Inf, label = "LR",
           vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
  scale_fill_manual(values = fill_4) +
  scale_x_discrete(labels = c(HR_T1 = "Pre", HR_T2 = "Post",
                               LR_T1 = "Pre", LR_T2 = "Post")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
  coord_cartesian(clip = "off") +
  labs(title = "1RM Leg Press", subtitle = anova_sub,
       y = "1RM (kg)", x = NULL, tag = "B") +
  FIG_THEME +
  theme(plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
        plot.margin = margin(5, 5, 20, 5), legend.position = "none")

delta_colors <- c(HR = unname(GROUP_FILL["HR_T2"]),
                  LR = unname(GROUP_FILL["LR_T2"]))

pB_right <- ggplot(pheno_wide, aes(x = Group, y = delta, fill = Group)) +
  geom_bar(stat = "summary", fun = mean, width = 0.55,
           color = "black", linewidth = 0.3) +
  geom_errorbar(stat = "summary", fun.data = mean_se,
                width = 0.15, linewidth = 0.4) +
  geom_jitter(width = 0.12, size = 1.2, alpha = 0.5,
              shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = delta_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(y = "\u0394 1RM (kg)", x = NULL) +
  FIG_THEME + theme(legend.position = "none")

pB <- (pB_left | pB_right) + plot_layout(widths = c(0.65, 0.35))

# Save ----
save_panel(pB, file.path(RPT, "panel_B_1rm_leg"), PW, PH)
