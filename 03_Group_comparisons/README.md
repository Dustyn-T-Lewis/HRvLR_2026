# 03_Group_comparisons

How do two groups of people differ? The same analyses run for every grouping from `02_Groups`:
the study's HR and LR, and the three clusterings (muscle size, muscle strength, proteome). Each
analysis is written once in `a_script/` and rendered once per grouping; the grouping folders hold
results only.

| Walkthrough | Does | Writes (per grouping) |
|---|---|---|
| `01_design_contrasts.qmd` | cell-means model with blood and person block; nine contrasts and their roles | `design.rds`, `design_contrasts.xlsx` |
| `02_protein_contrasts.qmd` | limma on limpa values: logFC with 95% CI, BH, Storey q, Π; fry concordance; PCA and PERMANOVA | `fit.rds`, `protein_contrasts.xlsx` |
| `03_set_contrasts.qmd` | fgsea and fry on every gene set; fgsea collapse | `set_tests.rds`, `set_contrasts.xlsx` |
| `04_module_contrasts.qmd` | limma on eigengenes, fry on modules, membership against protein change | `module_contrasts.xlsx` |
| `05_classification.qmd` | single-feature AUC with DeLong and bootstrap intervals, five tasks, three levels | `05_Classification/0{1,2,3}_*/*_classification.xlsx` |
| `06_trait_associations.qmd` | Spearman within each group; limma group × trait slope difference | `trait_associations.xlsx` |
| `07_module_preservation.qmd` | HR/LR only: WGCNA preservation of each group's modules in the other | `module_preservation.xlsx` |
| `07_cluster_phenotype.qmd` | Proteome only: do the proteome clusters differ in size or strength? | `cluster_phenotype.xlsx` |

The Proteome grouping runs each walkthrough once per clustering input (T1, T2, T3, training,
acute) into `Proteome/<step>/<input>/`, then stacks the five workbooks into
`Proteome/<step>/c_data/*_all_inputs.xlsx` with a `validity` column: `circular` when a comparison
uses exactly the clustering's time points, `caution` when it shares one, `clean` otherwise.

Every clustering found K = 1 (no natural groups); K = 2 was imposed, so every clustered comparison
is exploratory. HR and LR are extreme groups (top and bottom 8 of 25), so growth differences between
them are built in.

## Inputs

`02_Groups/<grouping>/c_data/groups.rds`, and the `01_Processing` outputs: `proteins.rds`,
`DAList_limpa.rds`, `DAList_corrected.rds`, `gene_sets.rds`, `singscore.rds`, `modules.rds`; plus
`00_Setup/Input/phenotype.csv`.

## Run

```sh
Rscript run_all.R 03_Group_comparisons
```

or one walkthrough for one grouping, from the project root:

```sh
quarto render 03_Group_comparisons/a_script/02_protein_contrasts.qmd -P grouping:HRvLR \
  --output-dir ../HRvLR/02_Protein_Contrasts/b_reports
quarto render 03_Group_comparisons/a_script/02_protein_contrasts.qmd -P grouping:Proteome \
  -P input:T1 --output-dir ../Proteome/02_Protein_Contrasts/T1/b_reports
```

Steps run in order: each reads the design (`01`) and, for sets and modules, the protein fit (`02`).

## Outputs

`<grouping>/<NN_Step>/b_reports/<walkthrough>.html` is the rendered walkthrough;
`<grouping>/<NN_Step>/c_data/` holds the workbooks and `.rds` files above. Every workbook opens with
`read_me` and `columns` and closes with `input_manifest` and `package_versions`.
