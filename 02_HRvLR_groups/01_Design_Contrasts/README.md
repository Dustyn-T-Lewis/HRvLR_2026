# 01_Design_Contrasts

Which comparisons does every later model make, and how is that model built?

## The idea

A linear model writes each protein's log2 abundance in each biopsy as a sum of terms. Here the
terms are six cell means and a blood slope. A cell mean is the model's average for one group at
one time point: HR_T1, HR_T2, HR_T3, LR_T1, LR_T2, LR_T3. The formula `~ 0 + group + blood` drops
the intercept, so each group column is a cell mean rather than a difference from a reference.

A contrast is a weighted sum of cell means that asks one question. `HR_T2 - HR_T1` asks how far
HR moved over training. The nine contrasts are written once here; every protein, set and module
fit reads them from `design.rds`, so every stage answers the same nine questions.

| Contrast | Role | Definition |
|---|---|---|
| `Training_HR` | descriptive | HR_T2 - HR_T1 |
| `Training_LR` | descriptive | LR_T2 - LR_T1 |
| `Acute_HR` | descriptive | HR_T3 - HR_T2 |
| `Acute_LR` | descriptive | LR_T3 - LR_T2 |
| `Baseline_HRvLR` | floor | HR_T1 - LR_T1 |
| `Trained_HRvLR` | descriptive | HR_T2 - LR_T2 |
| `Acute_HRvLR` | descriptive | HR_T3 - LR_T3 |
| `Training_Interaction` | primary | (HR_T2 - HR_T1) - (LR_T2 - LR_T1) |
| `Acute_Interaction` | secondary | (HR_T3 - HR_T2) - (LR_T3 - LR_T2) |

An interaction is a difference of differences: one group's change minus the other's. It is zero
when the two groups moved by the same amount, whatever that amount was, so it is the direct test
of whether the groups respond differently. A protein that changes in one group and not the other
has not, on that basis, shown a group difference; only the interaction tests that.

`Baseline_HRvLR` is the floor, not a negative control. HR and LR come from a median split of the
composite hypertrophy score, an outcome of training, so a difference at T1 shows how far the
split separates people before any training. Later group differences are read against it.

Subject is a block, not a fixed term. Each person gives up to three biopsies, and biopsies from
one person resemble each other. A fixed term per subject would absorb that, but every subject sits
wholly inside HR or LR, so it would also absorb every HR against LR comparison. Instead limpa's
`dpcDE()` estimates one within-subject correlation (how alike two biopsies from the same person
are once the model is fitted), and limma uses it to weight the biopsies. Within-person contrasts
gain from the pairing; between-group contrasts still compare people.

Blood enters as a covariate: each protein gets its own slope on the blood index (haemoglobin HBB
and HBA1, set in `01_Processing`), and every contrast is read at equal blood content.

## What the code does

| Step | Function | Detail |
|---|---|---|
| Design | `model.matrix(~ 0 + group + blood, meta)` | Six cell-mean columns and the blood index; row names are sample IDs. |
| Checks | `qr(design)$rank`, `make.names()`, `stopifnot()` | Stops if the design is rank deficient (a contrast would come back as silent NAs), if a column name is not a valid R name, if metadata and limpa columns are out of order, or if a group or subject is missing. |
| Contrasts | `limma::makeContrasts()` | The nine definitions above, named from their text. |
| Roles | `dplyr::case_when()` | primary, secondary, floor, descriptive. |
| Correlation | `limpa::dpcDE(proteins, design, block = subject)$correlation` | The same call `fit_limpa()` in `R/helpers.R` makes for the protein fit. |
| Save | `saveRDS()` | `design.rds` holds `dal` (data, annotation, metadata, and `design` with `design_matrix` and `contrast_matrix`, the field names proteoDA used), the contrast text, the roles and the floor's name. |

## How to read the outputs

- The within-subject correlation is one value for all proteins (limma's consensus). Near 0 the
  biopsies of one person behave as strangers; near 1 they are near copies. Later stages
  re-estimate it on their own matrices (`03_Set_Contrasts` and `04_Module_Contrasts` with
  `duplicateCorrelation()`), so their values can differ from this one.
- Roles say how a contrast is read, not how it is tested. All nine get the same fit and the same
  BH within contrast.

## Files

| | |
|---|---|
| Reads | `01_Processing/03_Quantification/c_data/proteins.rds` (limpa object and metadata) |
| Writes | `c_data/design.rds` (read by `02_Protein_Contrasts`, `03_Set_Contrasts`, `04_Module_Contrasts`), `c_data/design_contrasts.xlsx` (read by `Figures/F2_Programmes`) |
| Run | `Rscript 02_HRvLR_groups/01_Design_Contrasts/a_script/design_contrasts.R` |

```
01_Design_Contrasts/
├── a_script/design_contrasts.R
└── c_data/
    ├── design.rds
    └── design_contrasts.xlsx
```

No report PDF.

## Workbook: `design_contrasts.xlsx`

`read_me` lists the sheets; `columns` defines every column; `input_manifest` (files read, with
md5) and `package_versions` close the workbook.

| Sheet | One row is | Key columns |
|---|---|---|
| `contrasts` | a contrast | `definition`, `contrast`, `role` |
| `correlation_strata` | the fit | `level`, `correlation` (within subject, from limpa's fit) |
