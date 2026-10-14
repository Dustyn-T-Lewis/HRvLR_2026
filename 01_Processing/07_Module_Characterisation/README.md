# 07_Module_Characterisation

Says what each module is: gene sets its members over-represent (clusterProfiler ORA against the
measured genes), its ten highest-kME members, and whether its members interact more than chance in
STRING v12 (combined score at least 700).

## Inputs

| File | Holds |
|---|---|
| `01_Processing/06_Module_Building/c_data/modules.rds` | modules, membership and kME |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the tested sets and the gene universe |
| `references/STRING/` | STRING v12.0 links (score ≥ 700), aliases and protein info |

`references/README.md` gives the STRING source and version.

## Run

```sh
quarto render 01_Processing/07_Module_Characterisation/a_script/module_characterisation.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/module_characterisation.html` | the walkthrough: enrichment, top members, STRING test, networks |
| `c_data/module_characterisation.xlsx` | one label row per module, full ORA, top members, STRING test, drawn edges |
