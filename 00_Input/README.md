# 00 · Input

Study data. No code.

| File | One row is | Read by |
|---|---|---|
| `HRvLR_raw.xlsx` | one protein, 48 MS runs as columns | `01_Preprocess/01_Filtering` |
| `metadata.csv` | one MS sample, 48 rows | `01_Preprocess/01_Filtering` |
| `phenotype.csv` | one subject, 16 rows | every `04_Associate` step |
| `blood_contaminants.csv` | one curated blood protein, 95 rows | `01_Preprocess/01_Filtering` |
| `HPA_annotations_full.tsv` | one Human Protein Atlas gene, 20,162 rows | `01_Preprocess/01_Filtering` |
| `RBC_proteome_reference.tsv` | one red-cell protein, 7,012 rows | `01_Preprocess/01_Filtering` |
| `downloads/` | STRING v12 human links at score ≥ 700, aliases and protein info | `04_Network/05_Characterise` |

## STRING files are downloaded, not committed

They are too large for git. From the repo root:

```sh
mkdir -p 00_Input/downloads && cd 00_Input/downloads
curl -LO https://stringdb-downloads.org/download/stream/protein.links.v12.0/9606.protein.links.v12.0.min700.txt.gz
for f in aliases info; do
  curl -LO https://stringdb-downloads.org/download/protein.$f.v12.0/9606.protein.$f.v12.0.txt.gz
done
```

## HR and LR are the median split of the composite

`arm` is the median split of `comp_hypertrophy`: top 8 HR, bottom 8 LR. The composite is built from
fibre-area change, so a fibre-area difference between arms is built in and is not a finding.

## phenotype.csv holds pre and post side by side

Pre is the T1 biopsy and post the T2 biopsy; T3 has no phenotype. Units sit in each column name.
`comp_hypertrophy` is the source's composite and `volume_load_total_kg` is total kilograms lifted
over the programme, one value per subject. The extension 1RM is missing for one subject at one
timepoint.

The `fibres_*` columns are MyoVision fibre counts. The source workbook labels them "Number of
fCSA", and they fall as fibre area rises.

## Blood proteins are removed by identity

`blood_contaminants.csv` lists 95 proteins, each with its `class` and `reason`: 41 plasma, 20
immunoglobulin, 19 erythrocyte, 10 complement, 5 leukocyte. Every entry was checked against its
UniProt name. The `blood_cor`, `ery` and `myo` columns record the evidence at curation and decided
nothing.

The list replaced a cut at `blood_cor` 0.45. A permutation null for that correlation depends on how
many samples saw the protein: its 99.9th percentile is 0.71 for proteins seen in 15 to 25 samples
and 0.48 for 41 to 48. One cut cannot serve both. Sixteen proteins the old cut removed, all
expressed in myonuclei, stay in the matrix: RBMX, SYNE2, HNRNPA3, DDI2, GDI2, STK38, CPNE3, CCS,
FBXO7, LXN, ARF1, ANP32E, SEPTIN6, RAN, DNAJB4 and FLOT1.

## The red-cell reference flags and never removes

`RBC_proteome_reference.tsv` is the union of five red-cell proteomes. `sources` names the ones that
list each protein and `n_sources` counts them. A red-cell proteome shares about two thirds of a
muscle proteome, so `01_Filtering` reports membership as `in_rbc` and removes nothing on it.

| Source | Contributes |
|---|---|
| Téletchéa et al. 2019, *PLOS ONE* 14:e0211043 (RESPIRE) | 736 curated red-cell proteins |
| UniProtKB "Erythrocyte" annotation, retrieved 2026-07-23 | 597 entries |
| Bryk & Wiśniewski 2017, *J Proteome Res* 16:2752 | 2,577 quantified proteins |
| Ravenhill et al. 2019, *Commun Biol* 2:350 | 1,559 surface proteins |
| Bai et al. 2026, *Sci Data* 13, PXD067677 | 5,264 proteins |
