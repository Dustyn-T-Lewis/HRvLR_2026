# 04_Gene_Sets

Decides gene-set membership once for every set-level analysis: one protein per gene, MSigDB sets
cut to the measured genes, GO Slim terms built from the GO:BP hierarchy, and a 15-to-500 size rule.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/02_Normalisation/c_data/DAList_normalized.rds` | the measured proteins and their detection counts |
| `references/msigdb_2026.1.Hs.rds` | MSigDB 2026.1.Hs: Hallmark, Reactome, KEGG Legacy, GO:BP |
| `references/goslim_generic.obo` | the GO Consortium's generic GO Slim |

`references/README.md` gives the sources and dates of both snapshots.

## Run

```sh
quarto render 01_Processing/04_Gene_Sets/a_script/gene_sets.qmd --output-dir ../b_reports
```

Run from the project root, so R starts with the project's renv library.

## Outputs

| File | Holds |
|---|---|
| `b_reports/gene_sets.html` | the walkthrough: mapping, size rule, GO Slim construction |
| `c_data/gene_sets.rds` | tested sets, full catalogue, protein-to-gene map, gene universe; read by every set-level step |
| `c_data/gene_sets.xlsx` | the same catalogue and map as tables |
