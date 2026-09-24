# 00 · Input

Study data, and the two scripts that build derived inputs from it.

| File | One row is | Read by |
|---|---|---|
| `HRvLR_raw.xlsx` | one protein, 48 MS runs as columns | `01_Preprocess/01_Filtering` |
| `HRvLR_meta.csv` | one MS sample, 48 rows | `01_Filtering`, `a_script/01_build_phenotype.R` |
| `phenotype.csv` | one subject, 16 rows | `02_Differential_Expression/03_Phenotype`, `03_Pathway_Enrichment/05`, `04_Network/05` |
| `blood_contaminants.csv` | one curated blood protein, 95 rows | `01_Filtering` |
| `HPA_annotations_full.tsv` | one HPA gene, 20,162 rows | `01_Filtering` |
| `RBC_proteome_reference.tsv` | one red-cell protein, 7,012 rows | `01_Filtering` |
| `RBC_proteins.xlsx` | five published red-cell proteomes, one sheet each | `a_script/02_build_blood_list.R` |
| `downloads/` | STRING v12 human links at score ≥ 700, aliases and protein info | `04_Network/02` |

`PROVENANCE.md` records the red-cell sources and why the curated blood list replaced a correlation
cut.

## STRING files are downloaded, not committed

Too large for git. From the repo root:

```sh
mkdir -p 00_Input/downloads && cd 00_Input/downloads
curl -LO https://stringdb-downloads.org/download/stream/protein.links.v12.0/9606.protein.links.v12.0.min700.txt.gz
for f in aliases info; do
  curl -LO https://stringdb-downloads.org/download/protein.$f.v12.0/9606.protein.$f.v12.0.txt.gz
done
```

## Two builders derive phenotype.csv and the red-cell reference

```sh
Rscript 00_Input/a_script/01_build_phenotype.R    # HRvLR_meta.csv -> phenotype.csv
Rscript 00_Input/a_script/02_build_blood_list.R   # RBC_proteins.xlsx -> RBC_proteome_reference.tsv
```

Both outputs are committed; rerun a builder only when its source changes.

## HR/LR is the median split of COMP.HYPERTROPHY

`Group` in `HRvLR_meta.csv` is the median split of `COMP.HYPERTROPHY`: top 8 HR, bottom 8 LR. It
separates that composite's fibre-area ingredients by construction.

## Every trait but two is a T2 − T1 change

The MyoVision columns count fibres. Their names say "fCSA", but the source workbook calls them
"Number of fCSA". `phenotype.csv` names them `d_nfibre_mixed` and `d_nfibre_I`.

Phenotype exists at T1 and T2 only. `01_build_phenotype.R` computes each trait as T2 − T1, except
`comp_hypertrophy`, the source's composite read from the T2 row, and `volume_load`, total kilograms
lifted, recorded once per subject.
