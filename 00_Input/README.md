# 00 · Input

Study data, plus the two scripts that build derived inputs from it.

| File | What one row is | Read by |
|---|---|---|
| `HRvLR_raw.xlsx` | one protein, 48 MS runs as columns | `01_Preprocess/01_Filtering` |
| `HRvLR_meta.csv` | one MS sample, 48 rows | `01_Filtering`, `a_script/01_build_phenotype.R` |
| `phenotype.csv` | one subject, 16 rows | `02_Differential_Expression/03_Phenotype`, `03_Pathway_Enrichment/05`, `04_Network/03` |
| `blood_contaminants.csv` | one curated blood protein, 95 rows | `01_Filtering` |
| `HPA_annotations_full.tsv` | one HPA gene, 20,162 rows | `01_Filtering` |
| `RBC_proteome_reference.tsv` | one red-cell protein, 7,012 rows | `01_Filtering` |
| `RBC_proteins.xlsx` | five published red-cell proteomes, one sheet each | `a_script/02_build_blood_list.R` |
| `downloads/` | STRING v12 human links, aliases and protein info | `04_Network/02` |

`Phenotype.xlsx` and `HPA_skeletal_muscle_annotations.tsv` are the original workbook and an
older HPA export. No stage reads them.

## Get downloads/

Too large for git. From the repo root:

```sh
mkdir -p 00_Input/downloads && cd 00_Input/downloads
for f in links aliases info; do
  curl -LO https://stringdb-downloads.org/download/protein.$f.v12.0/9606.protein.$f.v12.0.txt.gz
done
```

## Derived inputs

```sh
Rscript 00_Input/a_script/01_build_phenotype.R    # HRvLR_meta.csv -> phenotype.csv
Rscript 00_Input/a_script/02_build_blood_list.R   # RBC_proteins.xlsx -> RBC_proteome_reference.tsv
```

Both outputs are committed; rerun only when a source changes. `PROVENANCE.md` records the
red-cell sources and why the curated blood list replaced a correlation cut.

## Notes

The HR/LR label is a median split. `Group` in `HRvLR_meta.csv` is the exact median split of
`COMP.HYPERTROPHY`: top 8 HR, bottom 8 LR. It separates that composite's fibre-area ingredients by
construction.

The MyoVision columns count fibres. Their names say "fCSA", but the source workbook calls them
"Number of fCSA". `phenotype.csv` names them `d_nfibre_mixed` and `d_nfibre_I`.

Phenotype exists at T1 and T2 only. Every trait in `phenotype.csv` is a T2 − T1 change except
`volume_load`, total kilograms lifted, recorded once per subject.
