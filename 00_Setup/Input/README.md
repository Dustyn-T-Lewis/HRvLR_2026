# Input

Study measurements and participant information. These files are read by the pipeline and are
not modified by it. Biological reference files live in the processing stages that use them.

## Files

| File | Contents | Used by |
|---|---|---|
| `HRvLR_raw.xlsx` | Protein quantities from a DIA-NN `pg.matrix` export: five annotation columns (`uniprot_id`, `protein`, `gene`, `description`, `n_seq`) and 48 MS-run columns | `01_Processing/01_Protein_Filtering` |
| `metadata.csv` | One row per MS run: `sample_id`, `subject`, `group`, `timepoint`, `group_time` | `01_Processing/01_Protein_Filtering` |
| `phenotype.csv` | One row per participant: pre/post muscle size and strength, fibre counts, composite hypertrophy score and training volume | Grouped and continuous analyses |

## Samples and participants

Sixteen participants each have three biopsies: T1 at baseline, T2 after training and T3 after
an acute bout at the end of training. The 48 runs represent repeated measurements of 16 people.
`sample_id` identifies a run (for example, `S03_T1`), `subject` identifies a participant
(`HR_S03`), and `group_time` combines group and time point (`HR_T1`). Sample IDs connect the workbook's
measurement columns to metadata; subject IDs connect metadata to phenotype.

All original runs stay here. Filtering decisions and exclusions belong to
[`01_Protein_Filtering`](../../01_Processing/01_Protein_Filtering/README.md).

## HR and LR are the two ends of the composite

25 people trained. `group` marks the top eight on `comp_hypertrophy` as HR and the bottom eight as
LR; the middle nine were not profiled. The composite is (2 × %ΔmCSA + %ΔfCSA) / 3, the
manuscript's formula, and `02_Groups` recomputes and checks it. Muscle and fibre growth define the
groups, so a growth difference between groups is built in and is not a finding.

## phenotype.csv holds pre and post side by side

Pre is the T1 biopsy and post the T2 biopsy; T3 has no phenotype. Units sit in each column name
(`fcsa_mixed_pre_um2`, `leg_press_1rm_post_kg`). `comp_hypertrophy` is the source's composite and
`volume_load_total_kg` the total kilograms lifted over the programme, one value per subject.
`leg_ext_1rm_post_kg` is missing for one subject.

The `fibres_*` columns are MyoVision fibre counts, labelled "Number of fCSA" in the source workbook.
`phenotype_tables()` in `00_Setup/Shared_functions/analysis_reporting_helpers.R` leaves them out of every association.

## Reference data

- [Protein filtering references](../../01_Processing/01_Protein_Filtering/references/README.md): curated blood contaminants, Human Protein Atlas and red-cell annotations.
- [Gene-set references](../../01_Processing/04_Gene_Sets/references/README.md): MSigDB and GO Slim.
- [Module references](../../01_Processing/07_Module_Characterisation/references/README.md): STRING v12.
