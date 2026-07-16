# 04_Figures — figure stages

Each figure is a self-contained `a_script/ b_reports/ c_data/` unit with its own
`01_run_<figure>.R` and narrative `.qmd`. Shared engines live in `functions/`
(style, palette, pathway and concordance helpers) and `shared/` (`pca.R`,
`clear_dir`, references); F04 owns the WGCNA module engine that F06 reads.

| Directory | Question | Engine |
| --- | --- | --- |
| `F00_PILOT/concordance_training/`, `concordance_acute/` | Where do HR and LR adapt alike over training and the acute bout, and where do they part? | Quadrant ORA, `limma::fry`, pathway NES concordance, RRHO2, bootstrap CI. |
| `F00_PILOT/summary/` | How large is the proteome response and how concordant are the two groups, before/after training and acute? | Median/p90 \|logFC\| per group × condition and Spearman ρ with Fisher-z CI. |
| `F01_phenotype/` | The phenotype: matched training, divergent growth and strength. | Phenotype atlas + linear mixed models. |
| `F02_proteome/` | Global proteome overview and QC. | PCA, DEP counts, effect sizes, set overlaps, η². |
| `F03_volcanoes/` | Per-contrast enrichment. | enrichVolcano ring-volcanoes, fgsea, EnrichmentMap dedup. |
| `F04_modules/` | Which WGCNA modules track the phenotype? | Signed WGCNA on the missForest-imputed proteome, `limma::fry`, LOSO q². |
| `F05_association/` | Which proteins and pathways associate with the continuous training responses (ΔmCSA, strength, ΔfCSA)? | Mixed models on proteins and singscore pathway scores, `feature ~ phenotype * timepoint + (1 \| subject)`. Association only. |
| `F06_prediction/prediction_responder/`, `prediction_continuous/` | Can baseline, training-response, or acute features predict HR vs LR — or the continuous phenotype — out of sample? | Elastic net (`glmnet`) + sparse PLS-DA (`mixOmics`), nested LOSO CV against a permutation null. `prediction_shared/` holds the harness both arms import. |

## The rules that hold throughout

**Prediction is scored against a permutation null, never against zero.** n = 16 makes a
single point estimate untrustworthy, so every number carries a nested-CV interval and a
permutation p. Q² goes negative when a model predicts worse than the mean, which is the
normal outcome at this sample size.

**Composite hypertrophy stays out of any model that also carries the HR/LR term**, since
the groups were defined from it. Carrying both would be conditioning on the outcome.

**Fold-specific transforms stay train-only.** The harness (`prediction_shared/_harness.R`)
z-scores and tunes every model inside the training subjects alone. singscore pathway scores
are single-sample and rank-based, so scoring the full matrix once is leakage-free; that is
why the suite uses singscore rather than a cohort-relative method. The one exception is the
F04 module eigengenes, which are fit on all 45 samples: a leave-one-out score computed on them
is transductive and optimistic rather than strictly out-of-sample. That is the documented
limitation behind the earlier positive module result, not a claim the eigengenes are honest.

## Nothing here survives multiple testing

Both prediction arms come back null after BH correction, and the concordance figures show
weak, non-significant agreement (ρ = 0.20 training, 0.15 acute; every `fry` test null). At
this sample size that null is the finding, not a failure.
