# Proteome

If people are grouped by their proteome alone, what else differs between the groups?

## What was clustered

The trunk's module eigengenes (`01_Trunk/06_Module_Building`), with no phenotype or group label,
five times over:

| Input | Each person's values |
|---|---|
| `T1`, `T2`, `T3` | eigengenes at that time point |
| `training` | eigengenes at T2 minus T1 |
| `acute` | eigengenes at T3 minus T2 |

The gap statistic selects K = 1 at every input. K = 2 is imposed at every input with the recorded
reason: so the exploratory arc can compare the two halves the proteome suggests; read each split
with its bootstrap Jaccard. The two clusters are `C1` and `C2`, numbered by PAM, so the names
carry no direction and do not match across inputs.

The method (PAM, gap statistic, silhouette, bootstrap Jaccard) is explained in
[`../README.md`](../README.md); the details for this arc are in
[`01_Clustering`](01_Clustering/README.md).

## Validity: circular, caution, clean

Here the groups come from the proteome, and the later stages test the proteome. Whether a
comparison means anything depends on which time points it shares with the clustering.
`validity()` in [`R/inputs.R`](R/README.md) tags every comparison against the input it was run on:

| Tag | Rule | Why |
|---|---|---|
| circular | the comparison uses exactly the clustering's time points | The clustering sorted people on those values, so the clusters differ there by construction, even in pure noise. The clustering made the difference; the result describes the split and tests nothing. |
| caution | it shares a time point with the clustering but is not the same comparison | Regression to the mean. The clustering sorted people partly on chance high or low values at the shared time point. That chance enters the comparison too, and at an unshared time point it does not repeat, so the clusters can seem to differ, or to change differently, with no real difference. |
| clean | it shares none | The noise the clustering sorted on does not enter the comparison. A difference here is a test of whether the split carries over to another time point or window. |

`uses()` reads the time points from the row's `contrast`, `task` or `window` name:

| Uses | Names |
|---|---|
| T1 | `Baseline_C1vC2` (contrast and classification task); the `baseline` and `baseline_level` windows; PERMANOVA window `T1` |
| T2 | `Trained_C1vC2` (contrast and task); `post_level`; PERMANOVA window `T2` |
| T3 | `Acute_C1vC2` (contrast and task); `post_bout_level`; PERMANOVA window `T3` |
| training (T1, T2) | `Training_C1`, `Training_C2`, `Training_Interaction`; the `Training_change_C1vC2` task; the `training` windows |
| acute (T2, T3) | `Acute_C1`, `Acute_C2`, `Acute_Interaction`; the `Acute_change_C1vC2` task; the `acute` windows |
| all (T1, T2, T3) | PERMANOVA windows `all` and `subject mean` |

Against each input that gives:

| Clustered at | T1 | T2 | T3 | training | acute | all |
|---|---|---|---|---|---|---|
| `T1` | circular | clean | clean | caution | clean | caution |
| `T2` | clean | circular | clean | caution | caution | caution |
| `T3` | clean | clean | circular | clean | caution | caution |
| `training` | caution | caution | clean | circular | caution | caution |
| `acute` | clean | caution | caution | caution | circular | caution |

Read circular rows as description, caution rows with the regression-to-the-mean risk in mind, and
clean rows as the tests. All of them inherit the imposed K: the clusters are halves of one spread,
to be read with their Jaccard, so even a clean result is exploratory.

## Per-input folders and stacked workbooks

Every stage from `02_Design_Contrasts` to `07_Trait_Associations` is a copy of the HR/LR stage
wrapped in `run_input(input)`. The script runs it once per input on that input's relabelled files
(`01_Clustering/c_data/inputs/<input>/`), and each run writes to its own folder:
`03_Protein_Contrasts/T1/c_data/protein_contrasts.xlsx`,
`03_Protein_Contrasts/training/b_reports/protein_contrasts_report.pdf`, and so on.

The script then calls `stack_inputs()` once, which joins the five per-input workbooks into
`<stage>/c_data/<content>_all_inputs.xlsx`, for example
`03_Protein_Contrasts/c_data/protein_contrasts_all_inputs.xlsx`. In it:

- every result sheet holds all five inputs, with `input` as the first column;
- `validity` follows `input` on any sheet with a `contrast`, `task` or `window` column, and is
  empty where the name maps to no time point;
- `read_me`, `columns`, `input_manifest`, `package_versions` and `plot_index` are not stacked;
  the stacked workbook writes its own, and its `input_manifest` lists the five per-input
  workbooks;
- `07_Trait_Associations` leaves `within_group` out, because five inputs pass xlsx's row limit;
  each input's own workbook keeps it.

Open the stacked workbook to compare inputs; open a per-input folder for the report PDF and any
sheet left out.

## Stages

| Stage | Asks | Copy of | Stacked workbook |
|---|---|---|---|
| [`01_Clustering`](01_Clustering/README.md) | how people split at each input | new | none: `clustering.xlsx` holds every input |
| [`02_Design_Contrasts`](02_Design_Contrasts/README.md) | which comparisons | `02_HRvLR_groups/01_Design_Contrasts` | `design_contrasts_all_inputs.xlsx` |
| [`03_Protein_Contrasts`](03_Protein_Contrasts/README.md) | which proteins differ | `02_HRvLR_groups/02_Protein_Contrasts` | `protein_contrasts_all_inputs.xlsx` |
| [`04_Set_Contrasts`](04_Set_Contrasts/README.md) | which gene sets shift | `02_HRvLR_groups/03_Set_Contrasts` | `set_contrasts_all_inputs.xlsx` |
| [`05_Module_Contrasts`](05_Module_Contrasts/README.md) | which modules shift | `02_HRvLR_groups/04_Module_Contrasts` | `module_contrasts_all_inputs.xlsx` |
| [`06_Classification/01_Proteins`](06_Classification/01_Proteins/README.md) | can one protein tell C1 from C2 | `02_HRvLR_groups/06_Classification/01_Proteins` | `protein_classification_all_inputs.xlsx` |
| [`06_Classification/02_Sets`](06_Classification/02_Sets/README.md) | can one gene-set score tell C1 from C2 | `02_HRvLR_groups/06_Classification/02_Sets` | `set_classification_all_inputs.xlsx` |
| [`06_Classification/03_Modules`](06_Classification/03_Modules/README.md) | can one module tell C1 from C2 | `02_HRvLR_groups/06_Classification/03_Modules` | `module_classification_all_inputs.xlsx` |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a feature's link to a trait differ by cluster | `02_HRvLR_groups/07_Trait_Associations` | `trait_associations_all_inputs.xlsx` |
| [`08_Cluster_Phenotype`](08_Cluster_Phenotype/README.md) | do the clusters differ in size or strength | new | none |
| [`Figures/F1_Overview`](Figures/F1_Overview/README.md) | the arc on a page, once per input | `02_HRvLR_groups/Figures/F1_Overview` | none |

`run_all.R` runs 01 to 08 in order, then the F1 script.

## Not tracked in git

`.gitignore` excludes `03_Cluster_groups/Proteome/**/b_reports/` and
`03_Cluster_groups/Proteome/**/c_data/`, with the note "five clustering inputs per stage,
regenerated by run_all.R". Only the scripts and READMEs are in git. Rebuild the outputs with
`Rscript run_all.R 03` (which also reruns the two muscle arcs) or stage by stage, then the F1
script.

## Files

```
Proteome/
├── R/inputs.R                    the five inputs, validity(), stack_inputs()
├── 01_Clustering/
│   ├── a_script/clustering.R
│   ├── b_reports/clustering_report.pdf
│   └── c_data/
│       ├── clustering.xlsx
│       └── inputs/<input>/       trunk inputs relabelled C1/C2, one folder per input
├── 02_Design_Contrasts/
│   ├── a_script/design_contrasts.R
│   ├── T1/c_data/  T2/c_data/  T3/c_data/  training/c_data/  acute/c_data/
│   └── c_data/design_contrasts_all_inputs.xlsx
├── 03_Protein_Contrasts/
│   ├── a_script/protein_contrasts.R
│   ├── <input>/b_reports/, <input>/c_data/   one pair per input
│   └── c_data/protein_contrasts_all_inputs.xlsx
├── 04_Set_Contrasts/             same layout
├── 05_Module_Contrasts/          same layout
├── 06_Classification/
│   ├── 01_Proteins/              same layout
│   ├── 02_Sets/                  same layout
│   └── 03_Modules/               same layout
├── 07_Trait_Associations/        same layout
├── 08_Cluster_Phenotype/
│   ├── a_script/cluster_phenotype.R
│   ├── b_reports/cluster_phenotype_report.pdf
│   └── c_data/cluster_phenotype.xlsx
└── Figures/F1_Overview/
    ├── a_script/f1_overview.R
    ├── b_reports/f1_overview_<input>.pdf
    └── c_data/f1_overview_<input>.xlsx
```
