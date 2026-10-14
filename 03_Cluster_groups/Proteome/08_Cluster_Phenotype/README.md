# 08_Cluster_Phenotype

Do people the proteome put in different clusters differ in muscle size or strength?

| | |
|---|---|
| Reads | `03_Cluster_groups/Proteome/01_Clustering/c_data/clustering.xlsx`, `00_Setup/Input/phenotype.csv` |
| Writes | `c_data/cluster_phenotype.xlsx`, `b_reports/cluster_phenotype_report.pdf` |
| Run | `Rscript 03_Cluster_groups/Proteome/08_Cluster_Phenotype/a_script/cluster_phenotype.R` |

## Key result

- 100 cluster-against-trait tests (Wilcoxon) over the five inputs: 4 at p < 0.05, 0 at FDR < 0.05.

## Method

- Per input, kind and trait: `wilcox.test(value ~ cluster)` on the training changes and the T1 and T2 levels; BH within input.
- `hr_lr`: HR and LR counts in each cluster with Fisher's exact p, a description.

## Read with care

- At most sixteen people split in two: a large p is not evidence of no difference.

## Sheets

`cluster_phenotype.xlsx` opens with `read_me` and `columns` (every column defined) and closes with `input_manifest` and `package_versions`.

| Sheet | Holds |
|---|---|
| `tests` | Per input, trait kind and trait: n and median per cluster, Wilcoxon W and p, BH within input. |
| `hr_lr` | Per input: HR and LR counts in each cluster, Fisher's exact p (description only). |
