# F0 Panel F: MyoVision Fiber Counts (Mixed + Type I) ----
# Mixed: sig interaction (p=0.029) -> show delta bracket
# Type I: no sig effects -> ANOVA subtitle only
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
PH <- 140

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
  )

fill_4 <- GROUP_FILL[c("HR_T1", "HR_T2", "LR_T1", "LR_T2")]
delta_colors <- c(HR = unname(GROUP_FILL["HR_T2"]),
                  LR = unname(GROUP_FILL["LR_T2"]))

# Helper: build Pre/Post + Delta panel pair for a fiber measure ----
build_fiber_panel <- function(meta_df, col_name, title, y_lab, tag) {
  sub <- meta_df %>% filter(!is.na(.data[[col_name]])) %>%
    rename(value = !!col_name)

  stats_anova <- rstatix::anova_test(data = sub, dv = value,
                                      wid = subject_key,
                                      between = Group, within = Timepoint)
  anova_tbl <- as.data.frame(stats_anova)
  anova_sub <- sprintf("Group %s   Time %s   Interaction %s",
                       fmt_p(anova_tbl$p[anova_tbl$Effect == "Group"]),
                       fmt_p(anova_tbl$p[anova_tbl$Effect == "Timepoint"]),
                       fmt_p(anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"]))

  wide <- sub %>% select(subject_key, Group, Timepoint, value) %>%
    pivot_wider(names_from = Timepoint, values_from = value) %>%
    mutate(delta = T2 - T1)

  hr_w <- wide %>% filter(Group == "HR")
  lr_w <- wide %>% filter(Group == "LR")
  tt_hr <- t.test(hr_w$T2, hr_w$T1, paired = TRUE)
  tt_lr <- t.test(lr_w$T2, lr_w$T1, paired = TRUE)
  tt_delta <- t.test(delta ~ Group, data = wide)

  interaction_sig <- anova_tbl$p[anova_tbl$Effect == "Group:Timepoint"] < 0.05
  time_sig        <- anova_tbl$p[anova_tbl$Effect == "Timepoint"] < 0.05

  sig_y <- bracket_pos(sub$value)
  p_left <- ggplot(sub, aes(x = Group_Time, y = value, fill = Group_Time)) +
    geom_bar(stat = "summary", fun = mean, width = 0.65,
             color = "black", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se,
                  width = 0.2, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 1.0, alpha = 0.5,
                shape = 21, color = "black", stroke = 0.3)

  if ((interaction_sig || time_sig) && tt_hr$p.value < 0.05) {
    p_left <- p_left +
      geom_signif(comparisons = list(c("HR_T1", "HR_T2")),
                  annotations = fmt_p(tt_hr$p.value),
                  y_position = sig_y, textsize = KEY_TEXT, tip_length = 0.01)
  }
  if ((interaction_sig || time_sig) && tt_lr$p.value < 0.05) {
    p_left <- p_left +
      geom_signif(comparisons = list(c("LR_T1", "LR_T2")),
                  annotations = fmt_p(tt_lr$p.value),
                  y_position = sig_y, textsize = KEY_TEXT, tip_length = 0.01)
  }

  p_left <- p_left +
    annotate("text", x = 1.5, y = -Inf, label = "HR",
             vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
    annotate("text", x = 3.5, y = -Inf, label = "LR",
             vjust = 4.2, fontface = "bold", size = 3.2, color = "grey25") +
    scale_fill_manual(values = fill_4) +
    scale_x_discrete(labels = c(HR_T1 = "Pre", HR_T2 = "Post",
                                 LR_T1 = "Pre", LR_T2 = "Post")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    coord_cartesian(clip = "off") +
    labs(title = title, subtitle = anova_sub,
         y = y_lab, x = NULL, tag = tag) +
    FIG_THEME +
    theme(plot.subtitle = element_text(size = 7, color = "grey40", face = "italic"),
          plot.margin = margin(5, 5, 20, 5), legend.position = "none")

  p_right <- ggplot(wide, aes(x = Group, y = delta, fill = Group)) +
    geom_bar(stat = "summary", fun = mean, width = 0.55,
             color = "black", linewidth = 0.3) +
    geom_errorbar(stat = "summary", fun.data = mean_se,
                  width = 0.15, linewidth = 0.4) +
    geom_jitter(width = 0.12, size = 1.0, alpha = 0.5,
                shape = 21, color = "black", stroke = 0.3)

  if (interaction_sig && tt_delta$p.value < 0.05) {
    p_right <- p_right +
      geom_signif(comparisons = list(c("HR", "LR")),
                  annotations = fmt_p(tt_delta$p.value),
                  textsize = KEY_TEXT, tip_length = 0.02,
                  y_position = bracket_pos(wide$delta))
  }

  p_right <- p_right +
    scale_fill_manual(values = delta_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
    labs(y = paste0("\u0394 ", y_lab), x = NULL) +
    FIG_THEME + theme(legend.position = "none")

  (p_left | p_right) + plot_layout(widths = c(0.65, 0.35))
}

# Build panels ----
pF_mixed <- build_fiber_panel(meta, "MyoVision_fCSA_mixed_Pre",
                               "MyoVision fCSA Mixed", "Fiber count", "F")
pF_type1 <- build_fiber_panel(meta, "MyoVision_fCSA_Type_I__Pre",
                               "MyoVision fCSA Type I", "Fiber count", NULL)

pF <- pF_mixed / pF_type1

# Audit ----
fiber_long <- meta %>%
  select(subject_key, Group, Timepoint, Group_Time,
         MyoVision_fCSA_mixed_Pre, MyoVision_fCSA_Type_I__Pre) %>%
  rename(Mixed  = MyoVision_fCSA_mixed_Pre,
         Type_I = MyoVision_fCSA_Type_I__Pre) %>%
  pivot_longer(cols = c(Mixed, Type_I),
               names_to = "fiber_type", values_to = "count") %>%
  filter(!is.na(count))

audit_F <- fiber_long %>%
  group_by(fiber_type, Group, Timepoint) %>%
  summarise(n = n(), mean = mean(count), sd = sd(count),
            sem = sd(count) / sqrt(n()), .groups = "drop")
write_csv(audit_F, file.path(DAT, "audit_panel_F_fiber_counts.csv"))

# Save ----
save_panel(pF, file.path(RPT, "panel_F_fiber_counts"), PW, PH)
