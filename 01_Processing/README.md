# 01_Processing

The shared steps every arc reads: from the raw protein report to cleaned, normalised, complete
protein values, gene-set scores and co-expression modules.

## The idea

The project asks its questions three ways (with the HR/LR labels, with groups found by clustering,
and with no groups), and all three must start from the same data. 01_Processing builds that data once.
Each arc reads the 01_Processing files and never rebuilds them, so a difference between arcs comes from
the question, not from the preprocessing.

01_Processing is built without the questions in mind. No step reads a phenotype or a contrast, and the
steps that shape the values do not read the HR/LR label:

- Blood is removed against time point only (`limma::removeBatchEffect`, design `~ 0 + timepoint`).
- Modules are defined on abundance centred within each subject, and never see a label.
- Gene sets are mapped and sized before any fold change exists.

One step does read the label: the detection filter in `01_Protein_Filtering` keeps a protein seen
in at least 5 samples of at least one of the six arm-by-time cells. It counts detections, not
abundance. A few descriptive outputs (eta squared in `02_Normalisation`, the blood-by-arm check
in `01_Protein_Filtering`, plots coloured by arm) use the label for display or checking only.

## The seven stages

| Stage | Question | Writes for the arcs |
|---|---|---|
| [`01_Protein_Filtering`](01_Protein_Filtering/README.md) | Which proteins and which runs are fit to analyse? | `DAList_filtered.rds` (to `02`) |
| [`02_Normalisation`](02_Normalisation/README.md) | How do we make the samples comparable? | `DAList_normalized.rds` (to `03`, `04`) |
| [`03_Quantification`](03_Quantification/README.md) | How do we give every protein a value in every sample without pretending the missing ones were measured? | `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds` |
| [`04_Gene_Sets`](04_Gene_Sets/README.md) | Which groups of proteins count as one process, decided once? | `gene_sets.rds` |
| [`05_Set_Scores`](05_Set_Scores/README.md) | How strongly does each biopsy express each set? | `singscore.rds` |
| [`06_Module_Building`](06_Module_Building/README.md) | Which proteins move together within a person, and how does each group move? | `modules.rds` |
| [`07_Module_Characterisation`](07_Module_Characterisation/README.md) | What is each module about? | `module_characterisation.xlsx` |

`run_all.R` runs them in this order, each in its own R session. A stage reads only from `00_Input`
and from earlier stages, through files on disk, so it reruns alone once its inputs exist.

`Figures/S01_QC` is a quality-control supplement. It reads `02_HRvLR_groups/02_Protein_Contrasts`
as well as 01_Processing, so `run_all.R` runs it with the figures at the end.

## How the data flow

```
00_Input/HRvLR_raw.xlsx, metadata.csv, blood and red-cell references
  └── 01_Protein_Filtering ── DAList_filtered.rds
        └── 02_Normalisation ── DAList_normalized.rds
              ├── 03_Quantification ── proteins.rds, DAList_limpa.rds, DAList_corrected.rds
              │     ├── 05_Set_Scores ── singscore.rds          (also reads gene_sets.rds)
              │     └── 06_Module_Building ── modules.rds
              │           └── 07_Module_Characterisation ── module_characterisation.xlsx
              │                                              (also reads gene_sets.rds, STRING)
              └── 04_Gene_Sets ── gene_sets.rds
```

The three matrices from `03_Quantification` differ in one way each:

| File | Values | Used where |
|---|---|---|
| `proteins.rds` | limpa's complete values with a standard error each | limma fits that weight each value by its precision (`limpa::dpcDE`) |
| `DAList_limpa.rds` | the same complete values, no blood correction | `limma::fry` set tests, which carry blood in their design, and the row order of overview heatmaps |
| `DAList_corrected.rds` | complete values with each protein's blood slope removed | methods that cannot carry a covariate: ranks, WGCNA, correlations, PCA |

## What each arc reads

| Arc | Reads from 01_Processing |
|---|---|
| `02_HRvLR_groups` | `proteins.rds` (design, protein contrasts), `DAList_limpa.rds`, `DAList_corrected.rds`, `gene_sets.rds`, `singscore.rds`, `modules.rds`, and `module_characterisation.xlsx` (figure F2) |
| `03_Cluster_groups` | Each `01_Clustering` copies `proteins.rds`, `DAList_limpa.rds`, `DAList_corrected.rds`, `singscore.rds` and `modules.rds` into its `c_data/inputs/`, relabelled with the cluster groups; later stages read those copies and `gene_sets.rds` directly. `Proteome/01_Clustering` clusters people on the module eigengenes. |
| `04_Continuous_no_group` | `DAList_limpa.rds`, `DAList_corrected.rds`, `gene_sets.rds`, `singscore.rds`, `modules.rds` |

No arc reads `DAList_filtered.rds` or `DAList_normalized.rds`.

## Folder

```
01_Processing/
├── README.md
├── 01_Protein_Filtering/
│   ├── a_script/protein_filtering.R
│   ├── b_reports/protein_filtering_report.pdf
│   └── c_data/  protein_filtering.xlsx, DAList_filtered.rds
├── 02_Normalisation/
│   ├── a_script/normalisation.R
│   ├── b_reports/  normalisation_report.pdf, norm_comparison.pdf, qc_pre.pdf, qc_post.pdf
│   └── c_data/  normalisation.xlsx, DAList_normalized.rds
├── 03_Quantification/
│   ├── a_script/quantification.R
│   ├── b_reports/quantification_report.pdf
│   └── c_data/  quantification.xlsx, proteins.rds, DAList_limpa.rds, DAList_corrected.rds
├── 04_Gene_Sets/
│   ├── a_script/gene_sets.R
│   └── c_data/  gene_sets.xlsx, gene_sets.rds, cache/ (MSigDB and GO Slim snapshots with .md5)
├── 05_Set_Scores/
│   ├── a_script/set_scores.R
│   ├── b_reports/set_scores_report.pdf
│   └── c_data/  set_scores.xlsx, singscore.rds
├── 06_Module_Building/
│   ├── a_script/module_building.R
│   ├── b_reports/module_building_report.pdf
│   └── c_data/  module_building.xlsx, modules.rds
├── 07_Module_Characterisation/
│   ├── a_script/module_characterisation.R
│   ├── b_reports/module_characterisation_report.pdf
│   └── c_data/module_characterisation.xlsx
└── Figures/
    └── S01_QC/
        ├── a_script/s01_qc.R
        ├── b_reports/s01_qc.pdf
        └── c_data/s01_qc.xlsx
```

Every stage holds one script in `a_script/`, its PDF report in `b_reports/` (`04_Gene_Sets` has
none), and in `c_data/` one workbook plus any `.rds` a later step reads. Every workbook opens on
`read_me` and `columns` (every column defined, from `R/column_dictionary.csv`) and ends with
`input_manifest` (inputs with md5) and `package_versions`.

## Run

```sh
Rscript run_all.R 01                                             # 01_Processing only
Rscript 01_Processing/03_Quantification/a_script/quantification.R     # one stage
```

`07_Module_Characterisation` needs the STRING files; `00_Input/README.md` has the download.
