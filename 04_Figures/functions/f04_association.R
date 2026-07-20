# Two association engines over the same features: a per-subject-summary limma
# regression (snapshot) and a three-timepoint mixed model (trajectory). Feature
# assembly for proteins and singscore pathways lives here so the four leaves
# stay thin. Association is descriptive; there is no train/test split.
pacman::p_load(dplyr, tibble, purrr, limma, lmerTest, parallel)

phase_subject_matrix <- function(mat, meta, phase,
                                 subject_col = "Subject_ID",
                                 time_col = "Timepoint") {
  tps <- switch(phase,
    baseline = "T1",
    training = c("T1", "T2"),
    acute = c("T2", "T3"),
    stop("phase must be baseline, training, or acute")
  )
  keep <- meta[[time_col]] %in% tps
  mat <- mat[, keep, drop = FALSE]
  subj <- meta[[subject_col]][keep]
  time <- meta[[time_col]][keep]
  if (phase == "baseline") {
    colnames(mat) <- subj
    return(mat)
  }
  early <- setNames(which(time == tps[1]), subj[time == tps[1]])
  late <- setNames(which(time == tps[2]), subj[time == tps[2]])
  both <- intersect(names(early), names(late))
  delta <- mat[, late[both], drop = FALSE] - mat[, early[both], drop = FALSE]
  colnames(delta) <- both
  delta
}

associate_limma <- function(mat, meta, pheno, traits,
                            phases = c("baseline", "training", "acute")) {
  purrr::map_dfr(phases, function(ph) {
    fmat <- phase_subject_matrix(mat, meta, ph)
    phr <- pheno[match(colnames(fmat), pheno$subject), , drop = FALSE]
    purrr::map_dfr(traits, function(tr) {
      y <- phr[[tr]]
      ok <- !is.na(y)
      fit <- limma::eBayes(limma::lmFit(
        fmat[, ok, drop = FALSE],
        model.matrix(~ y[ok])
      ))
      tt <- limma::topTable(fit, coef = 2, number = Inf, sort.by = "none")
      tibble::tibble(
        feature = rownames(tt), trait = tr, phase = ph,
        beta = tt$logFC, t = tt$t, p = tt$P.Value, bh = tt$adj.P.Val
      )
    })
  })
}

# One feature's interaction F-test, or NA when the feature is too sparse to fit.
# Random intercept only: a (timepoint | subject) random slope is the ideal for a
# trajectory test but is unidentifiable at 16 subjects x 3 timepoints, so the
# interaction F runs mildly anti-conservative. The HLM narratives state this.
.hlm_one <- function(v, resp, tp, subj) {
  d <- data.frame(score = v, response = resp, timepoint = tp, subject = subj)
  d <- d[!is.na(d$score) & !is.na(d$response), , drop = FALSE]
  d$timepoint <- droplevels(factor(d$timepoint))
  if (nlevels(d$timepoint) < 2 || dplyr::n_distinct(d$subject) < 4) {
    return(c(f = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_))
  }
  fit <- tryCatch(
    lmerTest::lmer(score ~ response * timepoint + (1 | subject), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(c(f = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_))
  }
  a <- suppressWarnings(stats::anova(fit))
  row <- grep("response:timepoint", rownames(a), fixed = TRUE)
  if (!length(row)) {
    return(c(f = NA_real_, df1 = NA_real_, df2 = NA_real_, p = NA_real_))
  }
  c(
    f = a[row, "F value"], df1 = a[row, "NumDF"], df2 = a[row, "DenDF"],
    p = a[row, "Pr(>F)"]
  )
}

associate_hlm <- function(mat, meta, pheno, traits,
                          subject_col = "Subject_ID", time_col = "Timepoint",
                          cores = 1L) {
  subj <- meta[[subject_col]]
  tp <- factor(meta[[time_col]])
  purrr::map_dfr(traits, function(tr) {
    resp <- setNames(pheno[[tr]], pheno$subject)[subj]
    rows <- parallel::mclapply(seq_len(nrow(mat)), function(i) {
      .hlm_one(mat[i, ], resp, tp, subj)
    }, mc.cores = cores)
    m <- do.call(rbind, rows)
    tibble::tibble(
      feature = rownames(mat), trait = tr,
      f_interaction = m[, "f"], df1 = m[, "df1"], df2 = m[, "df2"],
      p = m[, "p"], bh = stats::p.adjust(m[, "p"], "BH")
    )
  })
}
