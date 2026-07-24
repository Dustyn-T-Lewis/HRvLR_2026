# Mixed-model association fitter: builds the standard response ~ fixed +
# (1 | group) formula and returns the fitted lmerTest model, so every caller
# shares one model spec and keeps its own fixef/Satterthwaite extraction.
fit_hlm <- function(data, response, fixed = "Group * Timepoint",
                    group = "subject", reml = TRUE) {
  pacman::p_load(lmerTest)
  form <- stats::as.formula(
    sprintf("%s ~ %s + (1 | %s)", response, fixed, group)
  )
  lmerTest::lmer(form, data = data, REML = reml)
}

# Tidy sample table from a DAList: one row per sample with the design factors
# the global mixed model and its variance partition both read.
hlm_meta <- function(dal) {
  m <- dal$metadata
  data.frame(
    sample = m$Col_ID,
    subject = m$Subject_ID,
    group = factor(m$Group, levels = c("LR", "HR")),
    timepoint = factor(m$Timepoint, levels = c("T1", "T2", "T3")),
    stringsAsFactors = FALSE
  )
}

HLM_CONTRASTS <- c("group_main", "T1", "T2", "T3", "training", "acute")

# The six HR-LR contrasts read off the group x timepoint cell means: the three
# per-timepoint differences, the training and acute difference-in-differences,
# and the timepoint-averaged main effect. Weights are built from the emmGrid so
# they track its row order rather than a hard-coded cell sequence.
hlm_contrast_weights <- function(emm) {
  g <- as.data.frame(emm@grid)
  cell <- function(gr, tp) as.numeric(g$group == gr & g$timepoint == tp)
  hr_lr <- function(tp) cell("HR", tp) - cell("LR", tp)
  list(
    group_main = (hr_lr("T1") + hr_lr("T2") + hr_lr("T3")) / 3,
    T1 = hr_lr("T1"),
    T2 = hr_lr("T2"),
    T3 = hr_lr("T3"),
    training = (cell("HR", "T2") - cell("HR", "T1")) -
      (cell("LR", "T2") - cell("LR", "T1")),
    acute = (cell("HR", "T3") - cell("HR", "T2")) -
      (cell("LR", "T3") - cell("LR", "T2"))
  )
}

.global_one <- function(v, meta) {
  pacman::p_load(lmerTest, emmeans)
  d <- data.frame(
    y = v, group = meta$group, timepoint = meta$timepoint,
    subject = meta$subject
  )
  d <- d[!is.na(d$y), , drop = FALSE]
  d$timepoint <- droplevels(d$timepoint)
  d$group <- droplevels(d$group)
  ok <- nlevels(d$timepoint) == 3 && nlevels(d$group) == 2 &&
    dplyr::n_distinct(d$subject) >= 4
  fit <- if (ok) {
    tryCatch(
      lmerTest::lmer(y ~ group * timepoint + (1 | subject), data = d),
      error = function(e) NULL
    )
  } else {
    NULL
  }
  na_row <- data.frame(
    contrast = HLM_CONTRASTS, estimate = NA_real_, se = NA_real_,
    df = NA_real_, t = NA_real_, p = NA_real_
  )
  if (is.null(fit) || lme4::isSingular(fit)) {
    return(na_row)
  }
  emm <- emmeans::emmeans(fit, ~ group * timepoint, lmer.df = "satterthwaite")
  con <- as.data.frame(summary(
    emmeans::contrast(emm, method = hlm_contrast_weights(emm))
  ))
  data.frame(
    contrast = con$contrast, estimate = con$estimate, se = con$SE,
    df = con$df, t = con$t.ratio, p = con$p.value
  )
}

# One protein-wise mixed model across all samples, read out as the six HR-LR
# contrasts with Satterthwaite df, BH within contrast. Association only; F05
# adjudicates generalization.
associate_global_hlm <- function(mat, meta, cores = 1L) {
  meta <- meta[match(colnames(mat), meta$sample), , drop = FALSE]
  rows <- parallel::mclapply(seq_len(nrow(mat)), function(i) {
    cbind(feature = rownames(mat)[i], .global_one(mat[i, ], meta))
  }, mc.cores = cores)
  dplyr::bind_rows(rows) |>
    dplyr::group_by(.data$contrast) |>
    dplyr::mutate(bh = stats::p.adjust(.data$p, "BH")) |>
    dplyr::ungroup()
}

# Per-protein variance partition with every design factor as a random effect,
# the canonical variancePartition usage. Reports the fraction of each protein's
# variance from subject / group / timepoint / residual; at n=16 subject
# dominates, which is the honest headline (methods doc Part 9). Needs an
# NA-free matrix aligned to meta$sample.
varpart_global <- function(mat, meta) {
  pacman::p_load(variancePartition)
  meta <- meta[match(colnames(mat), meta$sample), , drop = FALSE]
  info <- data.frame(
    subject = meta$subject, group = meta$group, timepoint = meta$timepoint,
    row.names = meta$sample
  )
  form <- ~ (1 | subject) + (1 | group) + (1 | timepoint)
  vp <- variancePartition::fitExtractVarPartModel(mat, form, info)
  tibble::rownames_to_column(as.data.frame(vp), "feature")
}
