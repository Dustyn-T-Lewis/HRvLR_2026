# Input

Study measurements and participant information. The pipeline reads these files and never changes
them.

| File | Holds | Read by |
|---|---|---|
| `report.parquet` | DIA-NN 2 precursor report: 703,314 rows, 71 columns, 57 runs from 18 people, 2,729 protein groups, precursors at run q ≤ 0.01. Run names read `JBBC_CL_S03-T1`; re-injections carry `r`, `rr`, `_r` or `_4`; one run failed (S13-T1, 255 precursors) | `01_Processing/01_Precursors` |
| `metadata.csv` | One row per analysed biopsy (48): `sample_id`, `subject`, `group`, `timepoint`, `group_time` | `01_Processing/01_Precursors` |
| `phenotype.csv` | One row per person: muscle size and strength before and after training, fibre counts, composite hypertrophy score, training volume | the comparisons and continuous analysis, once rebuilt |

## People and biopsies

Sixteen people each gave three biopsies: T1 before training, T2 after training, T3 after one more
session at the end of training. `sample_id` names a biopsy (`S03_T1`), `subject` a person
(`HR_S03`), and `group_time` the group and time point (`HR_T1`). `01_Precursors` maps each report run
to a `sample_id`; `subject` links a biopsy to `phenotype.csv`. The report also holds S07 and S13,
who are not analysed.

## HR and LR are the two ends of the composite

Twenty-five people trained. `group` marks the top eight on `comp_hypertrophy` as HR and the bottom
eight as LR; the middle nine were not profiled. The composite is (2 × %ΔmCSA + %ΔfCSA) / 3, the
manuscript's formula. Muscle and fibre growth define the groups, so a growth difference between
them is built in and is not a finding.

## phenotype.csv holds pre and post side by side

Pre is the T1 biopsy and post the T2 biopsy; T3 has no phenotype. Units sit in each column name
(`fcsa_mixed_pre_um2`, `leg_press_1rm_post_kg`). `comp_hypertrophy` is the source's composite and
`volume_load_total_kg` the total kilograms lifted over the programme. `leg_ext_1rm_post_kg` is
missing for one person. The `fibres_*` columns are MyoVision fibre counts ("Number of fCSA" in the
source workbook), not outcomes.
