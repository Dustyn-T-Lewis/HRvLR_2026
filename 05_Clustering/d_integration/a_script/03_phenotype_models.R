# Module-eigengene vs muscle-phenotype models, three per engine x module
# x trait (tracks/baseline/acute), one BH family, one family-wise
# permutation null. Hypothesis-generating; n is tiny; no group_arm.
# The phenotype_models table reports the proper lmer estimate for the
# tracks model; the permutation null uses a subject-mean lm surrogate
# for tracks (and the same lm as the table for baseline/acute) so the
# null is tractable without refitting lmer thousands of times.

pacman::p_load(here, dplyr, tidyr, lmerTest, ggplot2)

traits <- c(
  "comp_hypertrophy", "d_fcsa_I", "d_fcsa_II", "d_myovision_fcsa_I",
  "d_mcsa", "d_1rm_legpress", "d_1rm_ext"
)
lead_trait <- "comp_hypertrophy"
min_n <- 6L
b_perm <- 1000L

pheno <- read.csv(here(
  "05_Clustering", "d_integration", "c_data", "phenotype.csv"
))

eigen <- bind_rows(
  read.csv(here(
    "05_Clustering", "a_wgcna_paired", "c_data", "eigengene.csv"
  )) |> mutate(
    engine = "wgcna", module = as.character(group_id),
    group_id = as.character(group_id)
  ),
  read.csv(here(
    "05_Clustering", "b_mfuzz_gap", "c_data", "eigengene.csv"
  )) |> mutate(
    engine = "mfuzz", module = paste0("c", group_id),
    group_id = as.character(group_id)
  )
)

zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

# Per-subject predictors derived from a module's eigengene series.
subject_predictors <- function(me_long) {
  me_long |>
    group_by(subject) |>
    summarise(
      mean_me = mean(ME, na.rm = TRUE),
      t1_me = ME[timepoint == "T1"][1],
      acute_shift = ME[timepoint == "T3"][1] - ME[timepoint == "T2"][1],
      .groups = "drop"
    )
}

# Fit one cell and return a one-row tibble, or NULL if skipped (n < 6).
fit_cell <- function(engine, module, trait, model, df) {
  if (model == "tracks") {
    d <- df |>
      filter(!is.na(ME), !is.na(.data[[trait]])) |>
      mutate(me_z = zscore(ME), trait_z = zscore(.data[[trait]]))
    n_subj <- dplyr::n_distinct(d$subject)
    if (n_subj < min_n || nrow(d) < min_n) {
      return(NULL)
    }
    fit <- lmerTest::lmer(me_z ~ trait_z + timepoint + (1 | subject), d)
    co <- summary(fit)$coefficients
    fitted_fixed <- as.numeric(model.matrix(fit) %*% lme4::fixef(fit))
    r2 <- var(fitted_fixed) / var(d$me_z)
    return(tibble(
      engine = engine, module = module, trait = trait, model = model,
      n = n_subj, beta_std = co["trait_z", "Estimate"], r2 = r2,
      p = co["trait_z", "Pr(>|t|)"]
    ))
  }

  pred <- if (model == "baseline") "t1_me" else "acute_shift"
  d <- df |>
    filter(!is.na(.data[[pred]]), !is.na(.data[[trait]])) |>
    mutate(x_z = zscore(.data[[pred]]), trait_z = zscore(.data[[trait]]))
  if (nrow(d) < min_n) {
    return(NULL)
  }
  fit <- lm(trait_z ~ x_z, d)
  co <- summary(fit)$coefficients
  tibble(
    engine = engine, module = module, trait = trait, model = model,
    n = nrow(d), beta_std = co["x_z", "Estimate"],
    r2 = summary(fit)$r.squared, p = co["x_z", "Pr(>|t|)"]
  )
}

cells <- list()
skipped <- list()
for (eng in unique(eigen$engine)) {
  mods <- unique(eigen$module[eigen$engine == eng])
  for (mod in mods) {
    me_long <- eigen |> filter(engine == eng, module == mod)
    sp <- subject_predictors(me_long)
    tracks_df <- me_long |> inner_join(pheno, by = "subject")
    subj_df <- sp |> inner_join(pheno, by = "subject")
    for (tr in traits) {
      for (md in c("tracks", "baseline", "acute")) {
        df <- if (md == "tracks") tracks_df else subj_df
        out <- fit_cell(eng, mod, tr, md, df)
        if (is.null(out)) {
          skipped[[length(skipped) + 1L]] <- tibble(
            engine = eng, module = mod, trait = tr, model = md
          )
        } else {
          cells[[length(cells) + 1L]] <- out
        }
      }
    }
  }
}

models <- bind_rows(cells) |>
  mutate(p_bh = p.adjust(p, method = "BH")) |>
  arrange(engine, module, trait, model)

write.csv(
  models,
  here("05_Clustering", "d_integration", "c_data", "phenotype_models.csv"),
  row.names = FALSE
)

# Family-wise permutation null. One surrogate |t| per cell, applied
# identically to observed and permuted phenotype tables: tracks uses
# subject-mean ME, baseline uses T1 ME, acute uses the T3-T2 shift.
# Each cell is an lm(trait ~ predictor) on complete cases (n >= 6).
# We permute phenotype rows across subjects (all traits kept together)
# and take the global max |t| across the whole grid as the statistic.

surrogate_cells <- list()
for (eng in unique(eigen$engine)) {
  mods <- unique(eigen$module[eigen$engine == eng])
  for (mod in mods) {
    me_long <- eigen |> filter(engine == eng, module == mod)
    sp <- subject_predictors(me_long)
    for (tr in traits) {
      for (md in c("tracks", "baseline", "acute")) {
        pred <- switch(md,
          tracks = "mean_me",
          baseline = "t1_me",
          acute = "acute_shift"
        )
        keep <- sp |>
          transmute(subject, x = .data[[pred]]) |>
          filter(!is.na(x))
        if (nrow(keep) < min_n) next
        surrogate_cells[[length(surrogate_cells) + 1L]] <- list(
          subject = keep$subject, x = keep$x, trait = tr
        )
      }
    }
  }
}

# |t| of the slope from an lm(y ~ x) given matched complete vectors.
abs_t_slope <- function(x, y) {
  ok <- !is.na(x) & !is.na(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)
  if (n < min_n) {
    return(NA_real_)
  }
  sxx <- sum((x - mean(x))^2)
  if (sxx == 0) {
    return(NA_real_)
  }
  b <- sum((x - mean(x)) * (y - mean(y))) / sxx
  resid <- y - (mean(y) - b * mean(x) + b * x)
  se <- sqrt(sum(resid^2) / (n - 2) / sxx)
  abs(b / se)
}

pheno_by_subject <- pheno[match(
  unique(pheno$subject), pheno$subject
), , drop = FALSE]
rownames(pheno_by_subject) <- pheno_by_subject$subject

grid_max_t <- function(pheno_tbl) {
  ts <- vapply(surrogate_cells, function(cell) {
    y <- pheno_tbl[cell$subject, cell$trait]
    abs_t_slope(cell$x, y)
  }, numeric(1))
  max(ts, na.rm = TRUE)
}

observed_max_t <- grid_max_t(pheno_by_subject)

subjects <- pheno_by_subject$subject
set.seed(1)
null_max_t <- vapply(seq_len(b_perm), function(i) {
  shuffled <- pheno_by_subject[sample(nrow(pheno_by_subject)), , drop = FALSE]
  rownames(shuffled) <- subjects
  grid_max_t(shuffled)
}, numeric(1))

p_perm <- mean(null_max_t >= observed_max_t)

perm_out <- tibble(
  null_max_t = null_max_t,
  observed_max_t = observed_max_t,
  p_perm = p_perm
)
write.csv(
  perm_out,
  here("05_Clustering", "d_integration", "c_data", "perm_null.csv"),
  row.names = FALSE
)

# Plots

cap <- paste0(
  "hypothesis-generating; n~15; no group_arm; pi-gated set partly ",
  "selected on the HR-vs-LR axis (effect sizes optimistically ",
  "biased) - permutation null is the guard."
)

grid_df <- models |>
  mutate(
    sig = p_bh < 0.05,
    model = factor(model, levels = c("baseline", "tracks", "acute"))
  )

p_grid <- ggplot(grid_df, aes(trait, module, fill = beta_std)) +
  geom_tile(color = "grey90") +
  geom_tile(
    data = filter(grid_df, sig),
    color = "black", linewidth = 0.8, fill = NA
  ) +
  facet_grid(engine ~ model, scales = "free_y", space = "free_y") +
  scale_fill_gradient2(
    low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0
  ) +
  labs(
    title = "Module eigengene vs muscle phenotype (signed std beta)",
    subtitle = cap, x = NULL, y = "module", fill = "beta_std"
  ) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  here("05_Clustering", "d_integration", "b_reports", "pheno_grid.pdf"),
  p_grid,
  width = 9, height = 6
)

null_band <- quantile(null_max_t, 0.95, na.rm = TRUE)

lead_df <- models |>
  filter(engine == "wgcna", trait == lead_trait) |>
  mutate(model = factor(model, levels = c("baseline", "tracks", "acute")))

p_lead <- ggplot(lead_df, aes(beta_std, module, color = model)) +
  annotate(
    "rect",
    xmin = -null_band, xmax = null_band, ymin = -Inf, ymax = Inf,
    alpha = 0.12, fill = "grey50"
  ) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  labs(
    title = paste0("WGCNA module effects on ", lead_trait),
    subtitle = cap, x = "beta_std", y = "module",
    caption = "shaded band = 95th pct of permutation-null |t| (as beta scale)"
  ) +
  theme_minimal(base_size = 9)

ggsave(
  here("05_Clustering", "d_integration", "b_reports", "lead_forest.pdf"),
  p_lead,
  width = 7, height = 4
)

message(
  "cells: ", nrow(models), " skipped: ", length(skipped),
  " BH<.05: ", sum(models$p_bh < 0.05, na.rm = TRUE),
  " observed max|t|: ", round(observed_max_t, 3),
  " p_perm: ", p_perm
)
