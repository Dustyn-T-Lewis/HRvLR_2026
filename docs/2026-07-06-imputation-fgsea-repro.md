# Imputation comparison, fresh fgsea, reproducibility (2026-07-06)

Shipped on `feature/clustering-phenotype`. This record keeps the design decisions;
the code and its history are the source of truth.

## What shipped

- **Perseus imputation arm** — `02_Normalization/imputation/a_script/d_perseus.R`,
  `proteoDA::perseus_impute(shift = 1.8, width = 0.3, robust = TRUE, seed = 42)`,
  mirroring `c_missforest.R`. Wired into `03_DEP/b_imputed` as a fourth method.
- **Fresh fgsea every run** — `HRvLR_F03_setup.R` recomputes enrichment on every
  render and overwrites `F03_source_data.xlsx`; the read-to-skip branch is gone
  and `set.seed(42)` sits directly above the fgsea loop. F04/F05 still read the
  written workbook.
- **F06 Panel A FDR marker** — asterisk on `p_bh < 0.05` cells alongside the raw-p
  outline.
- **HPA gate on accession** — `01_run_filtering.R` joins the HPA presence and
  contaminant gates on isoform-stripped, comma-split UniProt accession and logs
  the dropped identities to `hpa_absent_dropped.csv`.
- **S_imputation supplement** — `04_Figures/S_imputation/`: logFC concordance
  (Spearman vs the non-imputed arm) and π-significant counts across the four
  methods on the key contrasts.

## Design decisions

- **missForest stays canonical.** F04/F05/F06 read the missForest matrix; Perseus
  is comparison-only, the MNAR contrast to missForest's MAR assumption.
- **Imputation never feeds the reported DEP.** The primary limma fit runs on the
  non-imputed normalized matrix; the imputed arms exist to show the effects hold
  across missing-value strategy.
- **The fgsea workbook is an output, not a compute shortcut.** F03 owns the
  computation; downstream figures consume the file. Render order is F03 → F04/F05.
- **Every stochastic step is seeded at its call site** (WGCNA `randomSeed = 12345`,
  the F06 bootstrap, F02 PERMANOVA, missForest, fgsea, Perseus).

## Verification

CSV and PNG are the byte-reproducibility target; cairo PDFs and openxlsx workbooks
churn on timestamps and are compared by content. The accession join re-baselines
filtering through every downstream figure, so the pipeline re-renders once end to
end.
