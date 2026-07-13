# 05_test — proteome–phenotype signatures

How the muscle proteome relates to the hypertrophy phenotype, beyond the categorical
HR/LR split. Each unit is a self-contained `a_script/ b_reports/ c_data/` folder.

| Directory | Question | Engine |
| --- | --- | --- |
| `modules/` | Which co-expression modules exist, what do they couple to, and does a module still track the phenotype once the modules are refit inside the fold? | WGCNA fit, module–trait coupling, and an honest leave-one-subject-out refit. Feeds `04_Figures/F06`. |
| `phenotype_mapping/` | Which proteins and pathways associate with the continuous training responses (ΔmCSA, strength, ΔfCSA)? | Mixed models on proteins and singscore pathway scores, `feature ~ phenotype * timepoint + (1 \| subject)`. Association only; prediction is out of scope. |
| `prediction_responder_class/` | Can baseline, training-response, or acute features predict HR vs LR out of sample? | Elastic net (`glmnet`) + sparse PLS-DA (`mixOmics`), nested leave-one-subject-out CV with a permutation null. |
| `prediction_continuous_phenotype/` | Can those features predict the continuous phenotype out of sample? | Same engines, Gaussian; nested LOSO CV + permutation null. |
| `prediction_shared/` | — | The CV harness both prediction units import: feature builders, the nested-LOSO loop, the permutation null, and the panel code. |

## The rules that hold throughout

**Prediction is scored against a permutation null, never against zero.** n = 16 makes a
single point estimate untrustworthy, so every number carries a nested-CV interval and a
permutation p. Q² goes negative when a model predicts worse than the mean, which is the
normal outcome at this sample size.

**Composite hypertrophy stays out of any model that also carries the HR/LR term**, since
the groups were defined from it. Carrying both would be conditioning on the outcome.

**Features must be built inside the fold.** A module eigengene fit on all 45 samples has
already seen the held-out subject, so a leave-one-out score computed on it is transductive
and optimistic rather than out-of-sample. `modules/a_script/functions/honest_refit.R` is
the refit that answers this properly.

## Where this sits

`05_test` lives under `04_Figures/` because `F06` reads the module fit directly. Paths in
these scripts resolve through `here("04_Figures", "05_test", ...)`.
