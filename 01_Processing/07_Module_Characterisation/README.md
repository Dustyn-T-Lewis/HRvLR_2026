# 07_Module_Characterisation

Says what each module is: gene sets its members over-represent (clusterProfiler ORA against the
measured genes), its ten highest-kME members, and whether its members interact more than chance in
STRING v12 (combined score at least 700).

## Inputs

| File | Holds |
|---|---|
| `01_Processing/06_Module_Building/c_data/module_building.rds` | modules, membership and kME |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the tested sets (disease sets removed) and the gene universe |
| `references/STRING/` | STRING v12.0 links (score ≥ 700), aliases and protein info |

`references/README.md` gives the STRING source and version.

## Run

```sh
quarto render 01_Processing/07_Module_Characterisation/a_script/module_characterisation.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/module_characterisation.xlsx` | one label row per module, full ORA, top members, STRING test, drawn edges |
| `b_reports/module_characterisation.html` | the rendered walkthrough |

## Read with care

- The ORA universe is the measured genes, so a set is enriched relative to muscle, not the genome.
- The top set is the smallest p in the module across all collections; overlapping sets often
  share the top ranks. In three of the 10 modules (turquoise, brown, yellow) the top set does not
  reach BH < 0.05.
- STRING edges include text-mining and co-expression evidence; a high ratio says the members are
  known to be related, not that they bind.
- The STRING enrichment test calls the STRING web API; the step needs network access.
