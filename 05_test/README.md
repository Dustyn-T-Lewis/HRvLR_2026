# 05_test — proteome–phenotype signatures

Validation and exploration of how the muscle proteome relates to the hypertrophy
phenotype, beyond the categorical HR/LR split. Each analysis is a self-contained
`a_script/ b_reports/ c_data/` unit.

| Directory | Question | Engine |
| --- | --- | --- |
| `association_dream/` | Which proteins and pathways track the continuous phenotype and differ by responder class, across baseline, training, and acute states? | `variancePartition::dream` mixed models on proteins and singscore pathway scores, `feature ~ phenotype * timepoint + moderators + (1 \| subject)`, run continuously and for HR/LR |
| `prediction_responder_class/` | Can baseline, training-response, or acute features predict HR vs LR out of sample? | elastic net (`glmnet`) + sparse PLS-DA (`mixOmics`) on singscore features, nested leave-one-subject-out CV with a permutation null |
| `prediction_continuous_phenotype/` | Can those features predict the continuous phenotype (composite, ΔfCSA, ΔmCSA, Δ1RM) out of sample? | same engines, Gaussian; nested LOSO CV + permutation null |
| `singscore_vs_gsva_validation/` | Do rank-based singscore and cohort-relative GSVA scores agree? | score-agreement analysis; the earlier GSVA exploration lives in `a_script/legacy_gsva/` as the comparator |

Two rules hold throughout. Every prediction number carries a nested-CV interval
and a permutation p, because n = 16 makes single point estimates untrustworthy.
And composite hypertrophy stays out of any model that also carries the HR/LR term,
since the groups were defined from it.
