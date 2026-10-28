# 06_Classification

How well does one protein, gene-set score or module eigengene tell HR from LR? For each feature and
task, the AUC (Mann-Whitney W / (8 × 8)) with DeLong's 95% interval and the Wilcoxon p. Five tasks:
the T1, T2 and T3 levels and the training (T2 − T1) and bout (T3 − T2) changes. Nothing is fitted
and no AUC is measured on held-out people.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/03_Blood_Correction/c_data/proteins_corrected.rds` | blood-corrected protein values |
| `01_Processing/05_Set_Scores/c_data/singscore.rds` | singscore per set and biopsy |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the tested sets and their databases |
| `01_Processing/06_Module_Building/c_data/modules.rds` | eigengenes |
| `02_HRvLR/helpers.R` | `person_values()` |

## Run

```sh
quarto render 02_HRvLR/06_Classification/a_script/classification.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/classification.rds` | `auc` (every feature and task), `counts`, for the figures |
| `c_data/classification.xlsx` | the same tables |
| `b_reports/classification.html` | the rendered walkthrough |

## Read with care

- One test of 17,455 passes FDR < 0.05: the GO Slim set "transmembrane transport" in the training
  change (AUC 0.97). No protein or module does.
- The most extreme protein AUC is 0 (ACSL1 at T3: every LR person above every HR person), p =
  0.00016, the smallest p 8 against 8 can give.
- Proteins at p < 0.05 per task: 56 to 127 of 2,059, against about 103 by chance.
- FDR runs within feature type and task, and within database for sets.
