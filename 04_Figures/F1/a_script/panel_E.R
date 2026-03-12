# Figure 1 — Panel E: UpSet Plot (Dual Bar + Dot Matrix) + Legend Key ----

if (!exists("meta")) source("04_Figures/F1/a_script/HRvLR_F1_setup.R")

# Panel E depends on Panel D for sig_sets, dir_map, all_genes, pi_total, fdr_total
if (!exists("sig_sets")) source("04_Figures/F1/a_script/panel_D.R")

RPT <- RPT_DIR
DAT <- DAT_DIR

PE_W <- 200

# Analysis ----

SET_LABELS_INV <- setNames(CONTRASTS, unname(CTR_FACET[CONTRASTS]))

bin_mat <- sapply(sig_sets, function(s) as.integer(all_genes %in% s))
rownames(bin_mat) <- all_genes
colnames(bin_mat) <- CTR_FACET[colnames(bin_mat)]
label_to_contrast <- SET_LABELS_INV

cm <- make_comb_mat(bin_mat, mode = "intersect")
cs <- comb_size(cm)
keep <- cs > 0 & comb_degree(cm) > 0
cm_sub <- cm[keep]

comb_names_vec    <- comb_name(cm_sub)
n_comb            <- length(comb_names_vec)
set_names_ordered <- set_name(cm_sub)
up_counts    <- numeric(n_comb)
down_counts  <- numeric(n_comb)
mixed_counts <- numeric(n_comb)

for (i in seq_len(n_comb)) {
  members <- extract_comb(cm_sub, comb_names_vec[i])
  if (length(members) == 0) next
  bits <- as.logical(as.integer(strsplit(comb_names_vec[i], "")[[1]]))
  active_contrasts <- label_to_contrast[set_names_ordered[bits]]
  gene_dirs <- sapply(members, function(g) {
    dirs <- sapply(active_contrasts, function(ctr) {
      if (g %in% names(dir_map[[ctr]])) dir_map[[ctr]][g] else NA
    })
    dirs <- dirs[!is.na(dirs)]
    if (all(dirs == "Up")) "Up" else if (all(dirs == "Down")) "Down" else "Mixed"
  })
  up_counts[i]    <- sum(gene_dirs == "Up")
  down_counts[i]  <- sum(gene_dirs == "Down")
  mixed_counts[i] <- sum(gene_dirs == "Mixed")
}

n_mixed <- sum(mixed_counts)

display_total <- up_counts + down_counts
keep_display  <- display_total > 0
comb_ord <- which(keep_display)[order(-display_total[keep_display])]
up_ord   <- up_counts[comb_ord]
down_ord <- down_counts[comb_ord]

set_display_order <- unname(CTR_FACET[DISPLAY_ORDER])
n_int          <- length(comb_ord)
comb_names_ord <- comb_names_vec[comb_ord]
set_order_ch   <- set_name(cm_sub)
set_y_levels   <- rev(set_display_order)

comb_deg_ord <- vapply(comb_names_ord, function(cn) {
  sum(as.integer(strsplit(cn, "")[[1]]))
}, integer(1))

bar_long <- tibble(
  x         = rep(seq_len(n_int), 2),
  direction = factor(rep(c("Up", "Down"), each = n_int), levels = c("Down", "Up")),
  count     = c(up_ord, down_ord),
  is_single = rep(comb_deg_ord == 1, 2)
)

dot_df <- expand_grid(x = seq_len(n_int), set = set_display_order)
dot_df$active <- vapply(seq_len(nrow(dot_df)), function(r) {
  bits <- as.integer(strsplit(comb_names_ord[dot_df$x[r]], "")[[1]])
  idx <- match(dot_df$set[r], set_order_ch)
  if (is.na(idx)) return(FALSE)
  as.logical(bits[idx])
}, logical(1))
dot_df$set  <- factor(dot_df$set, levels = set_y_levels)
dot_df$ynum <- as.numeric(dot_df$set)

seg_list <- vector("list", n_int)
for (i in seq_len(n_int)) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  ypos <- match(set_order_ch[bits], set_y_levels)
  ypos <- ypos[!is.na(ypos)]
  if (length(ypos) > 1)
    seg_list[[i]] <- tibble(x = i, ymin = min(ypos), ymax = max(ypos))
}
seg_df <- bind_rows(seg_list)

stripe_fills <- adjustcolor(
  unname(CONTRAST_COLORS[CONTRASTS[match(set_y_levels, CTR_FACET[CONTRASTS])]]),
  alpha.f = 0.20
)

bar_bg_list <- lapply(seq_len(n_int), function(i) {
  bits <- as.logical(as.integer(strsplit(comb_names_ord[i], "")[[1]]))
  if (sum(bits) != 1) return(NULL)
  ctr_name <- CONTRASTS[match(set_order_ch[bits], CTR_FACET[CONTRASTS])]
  if (is.na(ctr_name)) return(NULL)
  tibble(xmin = i - 0.5, xmax = i + 0.5,
         fill = unname(CONTRAST_COLORS[ctr_name]))
})
bar_bg <- bind_rows(compact(bar_bg_list))

lbl_sz <- scale_text(BASE_COUNT, PE_W)

# Plots ----

pE_bars <- ggplot(bar_long, aes(x, count, fill = direction)) +
  {if (nrow(bar_bg) > 0)
    geom_rect(data = bar_bg,
              aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
              fill = bar_bg$fill, alpha = 0.20,
              color = "grey70", linewidth = 0.2, inherit.aes = FALSE)} +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           color = "black", linewidth = 0.3) +
  geom_text(data = \(d) d |> filter(count > 0, is_single),
            aes(label = count, y = count / 2),
            position = position_dodge(width = 0.7), vjust = 0.5,
            size = lbl_sz - 0.6, color = "white", fontface = "bold") +
  geom_text(data = \(d) d |> filter(count > 0, !is_single),
            aes(label = count, y = count + 2),
            position = position_dodge(width = 0.7), vjust = 0,
            size = lbl_sz - 0.6, color = "black", fontface = "bold") +
  scale_fill_manual(values = c(Up = unname(DIR_COLORS["Up"]),
                                Down = unname(DIR_COLORS["Down"]))) +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(y = "Intersection\nsize",
       title = "Contrast Overlap (UpSet)",
       subtitle = bquote(
         "DEPs by" ~ Pi < 0.05 * "; bars split by direction |" ~
         Pi * ":" ~ .(pi_total) ~ "/ FDR:" ~ .(fdr_total) ~
         "total (" * .(n_mixed) ~ "mixed-direction excluded)"
       ),
       tag = "E") +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "none",
        plot.margin = margin(2, 5.5, 0, 5.5))

n_sets <- length(set_y_levels)
pE_dots <- ggplot() +
  annotate("rect",
           xmin = 0.4, xmax = n_int + 0.6,
           ymin = seq_len(n_sets) - 0.40,
           ymax = seq_len(n_sets) + 0.40,
           fill = stripe_fills) +
  {if (nrow(seg_df) > 0)
    geom_segment(data = seg_df,
                 aes(x = x, xend = x, y = ymin, yend = ymax),
                 linewidth = 0.6, color = "grey25")} +
  geom_point(data = dot_df |> filter(!active),
             aes(x = x, y = ynum), color = "grey78", size = 1.2) +
  geom_point(data = dot_df |> filter(active),
             aes(x = x, y = ynum), color = "grey15", size = 1.2) +
  scale_x_continuous(expand = expansion(add = 0.3)) +
  scale_y_continuous(breaks = seq_len(n_sets), labels = set_y_levels,
                     expand = expansion(add = 0)) +
  coord_cartesian(ylim = c(0.55, n_sets + 0.45)) +
  labs(x = NULL, y = NULL) +
  FIG_THEME +
  theme(axis.text.x  = element_blank(), axis.ticks = element_blank(),
        panel.grid   = element_blank(),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.3),
        axis.text.y  = element_text(size = 6, face = "bold",
                                     margin = margin(r = 1)),
        plot.margin  = margin(0, 5.5, 0, 5.5))

# Legend key panel ----

ctr_leg <- tibble(
  y     = seq(from = (length(CONTRASTS) - 1) * 0.45, by = -0.45,
              length.out = length(CONTRASTS)),
  label = sprintf("%s:  %s", CTR_SHORT[CONTRASTS], CTR_FACET[CONTRASTS]),
  fill  = unname(CONTRAST_COLORS[CONTRASTS])
)

dir_leg <- tibble(
  y     = c(ctr_leg$y[2], ctr_leg$y[3]),
  label = c("Up", "Down"),
  fill  = unname(DIR_COLORS[c("Up", "Down")])
)

BOX_HALF <- 0.18
CTR_XMAX <- 5.2
DIR_XMIN <- 5.6
DIR_XMAX <- 6.5

pKeys <- ggplot() +
  # Contrast definition boxes (left side)
  geom_rect(data = ctr_leg,
            aes(xmin = 0, xmax = CTR_XMAX,
                ymin = y - BOX_HALF, ymax = y + BOX_HALF),
            fill = scales::alpha(ctr_leg$fill, 0.25),
            color = "grey70", linewidth = 0.2) +
  geom_text(data = ctr_leg,
            aes(x = 0.15, y = y, label = label),
            hjust = 0, size = KEY_TEXT, color = "grey15") +
  # Section header: Contrast Definitions
  annotate("text", x = 0.15, y = max(ctr_leg$y) + 0.42,
           label = "Contrast Definitions:", hjust = 0,
           size = KEY_TEXT + 0.3, fontface = "bold", color = "grey25") +
  # Direction boxes (right side)
  geom_rect(data = dir_leg,
            aes(xmin = DIR_XMIN, xmax = DIR_XMAX,
                ymin = y - BOX_HALF, ymax = y + BOX_HALF),
            fill = dir_leg$fill, color = "black", linewidth = 0.3) +
  geom_text(data = dir_leg,
            aes(x = DIR_XMIN + 0.15, y = y, label = label),
            hjust = 0, size = KEY_TEXT, color = "white", fontface = "bold") +
  # Section header: Direction
  annotate("text", x = DIR_XMIN + 0.15, y = max(ctr_leg$y) + 0.42,
           label = "Direction:", hjust = 0,
           size = KEY_TEXT + 0.3, fontface = "bold", color = "grey25") +
  coord_cartesian(xlim = c(-0.1, DIR_XMAX + 0.3),
                  ylim = c(min(ctr_leg$y) - 0.35,
                           max(ctr_leg$y) + 0.65)) +
  theme_void() +
  theme(plot.margin = margin(2, 5.5, 2, 5.5))

# Fisher's exact pairwise overlap ----

overlap_tests <- list()
contrast_pairs <- combn(CONTRASTS, 2, simplify = FALSE)
n_bg <- length(all_genes)
for (pair in contrast_pairs) {
  a <- sig_sets[[pair[1]]]; b <- sig_sets[[pair[2]]]
  n_both <- length(intersect(a, b))
  n_a    <- length(a); n_b <- length(b)
  mat <- matrix(c(n_both, n_a - n_both, n_b - n_both,
                  n_bg - n_a - n_b + n_both), nrow = 2)
  ft <- fisher.test(mat, alternative = "greater")
  overlap_tests[[paste(pair, collapse = " & ")]] <- data.frame(
    set_A = pair[1], set_B = pair[2],
    n_A = n_a, n_B = n_b, overlap = n_both,
    expected = round(n_a * n_b / n_bg, 1),
    odds_ratio = round(ft$estimate, 2), p_value = ft$p.value
  )
}
overlap_df <- bind_rows(overlap_tests)
overlap_df$p_bh <- p.adjust(overlap_df$p_value, method = "BH")
write_csv(overlap_df, file.path(DAT, "audit_panel_E_overlap_enrichment.csv"))

pE_combined <- (pE_bars / pE_dots / pKeys) +
  plot_layout(heights = c(0.58, 0.22, 0.20))

save_panel(pE_combined, file.path(RPT, "panel_E_upset"),
           width = PE_W, height = 160)
