# 02_Trait_Pathways

Runs fgsea on each window-and-trait protein ranking (limma t from the change scores), one database
at a time, and picks the sets the trait pathway figure draws: candidates at FDR < 0.05 against any
trait of a window, grouped by a Jaccard tree as in `02_HRvLR/02_Set_Contrasts`.

## Inputs

| File | Holds |
|---|---|
| `03_Continuous/01_Trait_Associations/c_data/trait_associations.rds` | protein change scores per window and trait |
| `01_Processing/04_Gene_Sets/c_data/gene_sets.rds` | the tested sets and the protein-to-gene map |

## Run

```sh
quarto render 03_Continuous/02_Trait_Pathways/a_script/trait_pathways.qmd --output-dir ../b_reports
```

## Outputs

| File | Holds |
|---|---|
| `c_data/trait_pathways.rds` | `pathways` (every set, window and trait), `trait_selection` (candidates, tree cuts, representatives) |
| `c_data/trait_pathways.xlsx` | `counts`, `trait_pathways`, `trait_selection` |
| `b_reports/trait_pathways.html` | the rendered walkthrough |

## Read with care

- 3,630 set calls at FDR < 0.05 across 42 rankings, against 3 proteins at FDR < 0.05 in the same
  rankings. fgsea treats proteins as independent, so co-moving sets look stronger than 16 people
  support.
- FDR runs within one database and one ranking. GO Slim is tested but kept out of the selection.
- fgsea could not estimate p for 12 set-ranking pairs, 10 of them in the training window (energy
  metabolism and transmembrane transport sets against the composite score, fibre area and muscle
  area). Their p, NES and FDR stay missing; the report lists them.
- 308 sets represent the six windows at cut 0.85 (32 in training to 80 in the T2 level window).
