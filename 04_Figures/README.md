# extras — supporting analyses

Units that support the main figures without being one. Each is a self-contained
`a_script/ b_reports/ c_data/` folder.

| Directory | Question | Engine |
| --- | --- | --- |
| `association/` | Which proteins and pathways associate with the continuous training responses (ΔmCSA, strength, ΔfCSA)? | Mixed models on proteins and singscore pathway scores, `feature ~ phenotype * timepoint + (1 \| subject)`. Association only; prediction is out of scope. |
| `prediction_responder_class/` | Can baseline, training-response, or acute features predict HR vs LR out of sample? | Elastic net (`glmnet`) + sparse PLS-DA (`mixOmics`), nested leave-one-subject-out CV against a permutation null. |
| `prediction_continuous_phenotype/` | Can those features predict the continuous phenotype out of sample? | Same engines, Gaussian; nested LOSO CV + permutation null. |
| `prediction_shared/` | — | The CV harness both prediction units import: feature builders, the nested-LOSO loop, the permutation null, and the panel code. |
| `concordance_training/` | Where do HR and LR adapt alike over training, and where do they part? | Quadrant ORA, `limma::fry`, pathway NES concordance, RRHO2, bootstrap CI. |
| `concordance_acute/` | The same question for the acute bout. | Same driver, different contrast pair. |
| `imputation/` | Do the four imputation arms reshape the effects? | Non-imputed reference against `imp4p`, MsCoreUtils, `missForest`, and the Perseus MNAR draw. |

The WGCNA module engine these units read lives at `04_Figures/modules/`. F04 reads it too.

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

## Nothing here survives multiple testing

Both prediction arms come back null after BH correction, and the concordance figures show
weak, non-significant agreement (ρ = 0.20 training, 0.15 acute; every `fry` test null). At
this sample size that null is the finding, not a failure.
