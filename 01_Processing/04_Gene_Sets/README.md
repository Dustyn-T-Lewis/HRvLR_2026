# 04_Gene_Sets

Decides gene-set membership once for every set-level analysis: one protein per gene, MSigDB sets
cut to the measured genes, GO Slim terms built from the GO:BP hierarchy, a 15-to-500 size rule, and
disease sets removed by name.

## Inputs

| File | Holds |
|---|---|
| `01_Processing/02_Quantification/c_data/proteins.rds` | the measured proteins with gene symbol, precursor count (`NPrec`) and observed share (`PropObs`) |
| `references/msigdb_2026.1.Hs.rds` | MSigDB 2026.1.Hs: Hallmark, Reactome, KEGG Legacy, GO:BP |
| `references/goslim_generic.obo` | the GO Consortium's generic GO Slim |

`references/README.md` gives the sources and dates of both snapshots.

## Run

```sh
quarto render 01_Processing/04_Gene_Sets/a_script/gene_sets.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/gene_sets.rds` | `sets` (tested sets), `set_catalog` (every set with `size_ok`, `disease_set`, `qualifies`), `protein_map`, `gene_universe`; read by every set-level step |
| `c_data/gene_sets.xlsx` | collection summary, the disease sets removed, the catalogue, the protein-to-gene map |
| `b_reports/gene_sets.html` | the rendered walkthrough |

## Read with care

- 1,422 sets are tested: Hallmark 35, Reactome 321, KEGG Legacy 55, GO:BP 953, GO Slim 58.
- 54 sets that pass the size rule are removed as disease sets (Reactome 32, KEGG 12, GO:BP 10). The
  rule reads the set's name only, not its genes; the list is the `disease_excluded` sheet.
- The rule removes every set whose name carries "virus" or "viral", including host-defence sets
  ("response to virus") and host-machinery sets ("viral gene expression", "viral process",
  Reactome's influenza viral RNA transcription set, "symbiont entry into host"). Add a name to
  `not_disease` in the script to bring one back.
- GO Slim sets nest over GO:BP and collect genes from every term below them, including terms whose
  own sets were removed. They are an overview family, not independent tests.
- One symbol (MOCS2) has two proteins; the one with the higher observed share represents it.
