# Does pathway activity track or predict the hypertrophy phenotype?
#  (a) mixed model score ~ phenotype * timepoint + (1|subject) - does the
#      trajectory track the continuous response;
#  (b) baseline (T1) score -> phenotype with leave-one-out CV, referenced to a
#      permutation null. GSVA scores are cohort-relative (Bhuva 2020) and at n=16
#      the best pathway is selection-biased (Ambroise & McLachlan 2002), so the
#      null is what makes the prediction claim honest.
if (!exists("expr")) source(here::here("05_test", "a_script", "setup.R"))
scores <- read_csv(file.path(DAT, "gsva_scores.csv"), show_col_types = FALSE) |>
  tibble::column_to_rownames("pathway") |>
  as.matrix()

TRAITS <- c("comp_hypertrophy", "d_mcsa", "d_1rm_legpress")
long <- scores |>
  as.data.frame() |>
  tibble::rownames_to_column("pathway") |>
  pivot_longer(-pathway, names_to = "sample", values_to = "score") |>
  left_join(meta, by = "sample") |>
  left_join(pheno |> select(subject, all_of(TRAITS)), by = "subject")

# (a) phenotype-trajectory association
fit_pheno <- function(pw_name, trait) {
  d <- long[long$pathway == pw_name, ]
  tryCatch(
    {
      f <- as.formula(sprintf("score ~ %s * timepoint + (1|subject)", trait))
      a <- anova(lmerTest::lmer(f, data = d))
      tibble(
        pathway = pw_name, trait = trait,
        p_pheno = a[trait, "Pr(>F)"],
        p_pheno_time = a[sprintf("%s:timepoint", trait), "Pr(>F)"]
      )
    },
    error = function(e) NULL
  )
}
pheno_lmm <- bind_rows(lapply(TRAITS, function(tr) {
  bind_rows(lapply(rownames(scores), fit_pheno, trait = tr))
}))
pheno_lmm <- pheno_lmm |>
  group_by(trait) |>
  mutate(fdr_pheno_time = p.adjust(p_pheno_time, "BH")) |>
  ungroup()
write_csv(pheno_lmm, file.path(DAT, "lmm_phenotype.csv"))

# (b) baseline-prediction with LOO-CV
loo_q2 <- function(x, y) {
  d <- data.frame(x = x, y = y)
  pred <- vapply(seq_len(nrow(d)), function(i) {
    predict(lm(y ~ x, data = d[-i, ]), newdata = d[i, ])
  }, numeric(1))
  1 - sum((d$y - pred)^2) / sum((d$y - mean(d$y))^2)
}

base <- long |>
  filter(timepoint == "T1") |>
  select(pathway, subject, score, all_of(TRAITS))
pred <- bind_rows(lapply(TRAITS, function(tr) {
  bind_rows(lapply(unique(base$pathway), function(p) {
    d <- base[base$pathway == p, ]
    d <- d[!is.na(d$score) & !is.na(d[[tr]]), ]
    if (nrow(d) < 8) {
      return(NULL)
    }
    tibble(pathway = p, trait = tr, n = nrow(d), q2 = loo_q2(d$score, d[[tr]]))
  }))
}))

# permutation null: best Q2 over all pathways under shuffled phenotype
set.seed(42)
null_q95 <- vapply(TRAITS, function(tr) {
  w <- base |>
    select(pathway, subject, score) |>
    pivot_wider(names_from = subject, values_from = score)
  subj <- intersect(pheno$subject[!is.na(pheno[[tr]])], names(w)[-1])
  y <- pheno[[tr]][match(subj, pheno$subject)]
  xmat <- as.matrix(w[, subj])
  maxq <- replicate(200, {
    yp <- sample(y)
    qs <- apply(xmat, 1, function(xr) {
      ok <- !is.na(xr)
      if (sum(ok) < 8) NA else loo_q2(xr[ok], yp[ok])
    })
    max(qs, na.rm = TRUE)
  })
  unname(quantile(maxq, 0.95))
}, numeric(1))
pred$null_q95 <- null_q95[pred$trait]
pred$beats_null <- pred$q2 > pred$null_q95
write_csv(pred, file.path(DAT, "prediction.csv"))

trait_labs <- c(
  comp_hypertrophy = "Composite hypertrophy",
  d_mcsa = "delta mCSA", d_1rm_legpress = "delta 1RM leg press"
)
top_pred <- pred |>
  group_by(trait) |>
  slice_max(q2, n = 8, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    trait = factor(trait_labs[trait], levels = trait_labs),
    term = reorder_within(clean_pathway_name(pathway, 30), q2, trait)
  )
fig_pred <- ggplot(top_pred, aes(q2, term, fill = beats_null)) +
  geom_col(width = 0.7) +
  geom_vline(aes(xintercept = null_q95), linetype = "dotted", color = "grey40") +
  geom_vline(xintercept = 0, color = "grey70", linewidth = 0.3) +
  scale_fill_manual(
    values = c(`TRUE` = "#1B7837", `FALSE` = "grey75"),
    labels = c(`TRUE` = "beats null", `FALSE` = "not above chance"), name = NULL
  ) +
  facet_wrap(~trait, scales = "free", ncol = 1) +
  labs(
    x = "leave-one-out Q-squared (baseline GSVA -> phenotype)", y = NULL,
    title = "Baseline pathway activity rarely predicts the response",
    subtitle = sprintf(
      "Dotted line = 95th percentile of the permutation null (best of %d pathways)",
      nrow(scores)
    )
  ) +
  FIG_THEME +
  scale_y_reordered() +
  theme(axis.text.y = element_text(size = 6))
save_panel(fig_pred, file.path(RPT, "fig3_prediction"), 170, 190)

cat(sprintf(
  "phenotype: %d pathway-trait baseline models; beating null: %d\n",
  nrow(pred), sum(pred$beats_null)
))
