test_that("fit_limpa is limpa's precision-weighted fit, contrasted and robustly moderated", {
  set.seed(5)
  peptides <- limpa::simProteinDataSet(n.peptides = 200, samples.per.group = 6)
  proteins <- suppressMessages(limpa::dpcQuant(peptides, "Protein", dpc = limpa::dpc(peptides)))
  group <- factor(paste0("G", peptides$targets$Group))
  subject <- rep(paste0("P", 1:6), 2)
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contrasts <- makeContrasts(G2 - G1, levels = design)

  eb <- fit_limpa(proteins, design, contrasts, subject)
  manual <- limpa::dpcDE(proteins, design, block = subject, plot = FALSE) |>
    contrasts.fit(contrasts) |>
    eBayes(robust = TRUE)
  expect_equal(eb$coefficients, manual$coefficients)
  expect_equal(eb$p.value, manual$p.value)

  changed <- unique(peptides$genes$Protein[peptides$genes$DEStatus != "NotDE"])
  p <- eb$p.value[, 1]
  expect_lt(median(p[names(p) %in% changed]), median(p[!names(p) %in% changed]))
})

# removeBatchEffect centres the covariate, so the correction leaves each feature's mean unchanged.
test_that("blood_corrected removes each feature's blood slope fitted beside the design", {
  set.seed(6)
  blood <- rnorm(10)
  time <- factor(rep(c("T1", "T2"), 5))
  design <- model.matrix(~ 0 + time)
  values <- rbind(a = 3 * blood + as.numeric(time) + rnorm(10, sd = 0.2), b = rnorm(10))
  out <- blood_corrected(values, blood, design)
  slope <- coef(lm(values["a", ] ~ 0 + time + blood))[["blood"]]
  expect_equal(unname(out["a", ]), unname(values["a", ] - slope * (blood - mean(blood))))
  expect_equal(rowMeans(out), rowMeans(values))
})
