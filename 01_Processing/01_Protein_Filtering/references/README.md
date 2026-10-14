# Protein filtering references

These reference files are read by `../a_script/protein_filtering.R`. The curated blood list
is part of the filtering specification; preserve it with the analysis code. The HPA and RBC
tables are retained reference snapshots.

## Blood proteins are removed by identity

`blood_contaminants.csv` lists 95 proteins, each with its `class` and `reason`: 41 plasma, 20
immunoglobulin, 19 erythrocyte, 10 complement, 5 leukocyte. Every entry was checked against its
UniProt name. `01_Protein_Filtering` reads only `uniprot_id`, `gene`, `class` and `reason`; the
`blood_cor`, `ery` and `myo` columns record the evidence at curation and decide nothing.

No single cut on `blood_cor` works: its permutation null depends on how many samples saw the
protein (sheet `blood_cor_null` in `01_Protein_Filtering`). Sixteen proteins that a cut at 0.45
would remove, all expressed in myonuclei, are not on the list: RBMX, SYNE2, HNRNPA3, DDI2, GDI2, STK38, CPNE3, CCS, FBXO7, LXN, ARF1, ANP32E,
SEPTIN6, RAN, DNAJB4 and FLOT1.

## The Human Protein Atlas table tags blood proteins the list misses

`01_Protein_Filtering` reads five HPA columns: `Protein class` (immunoglobulin genes),
`Secretome location` ("Secreted to blood"), `Blood concentration - Conc. blood MS [pg/L]`, and the
single-cell RNA levels in erythrocytes and in myonuclei (nCPM). The first three tag a protein as
blood; the myonuclei level and the blood concentration can rescue it as muscle.

## The red-cell reference flags and never removes

`RBC_proteome_reference.tsv` is the union of five red-cell proteomes. `sources` names the ones that
list each protein and `n_sources` counts them. Red cells share much of their proteome with muscle,
so `01_Protein_Filtering` reports membership as `in_rbc` and removes nothing on it.

| Tag in `sources` | Source | Proteins |
|---|---|---|
| `RESPIRE` | Téletchéa et al. 2019, *PLOS ONE* 14:e0211043 | 736 curated red-cell proteins |
| `Uniprot` | UniProtKB "Erythrocyte" annotation, retrieved 2026-07-23 | 597 entries |
| `JPR2017` | Bryk & Wiśniewski 2017, *J Proteome Res* 16:2752 | 2,577 quantified proteins |
| `CB2019` | Ravenhill et al. 2019, *Commun Biol* 2:350 | 1,559 surface proteins |
| `This study` | Bai et al. 2026, *Sci Data* 13, PXD067677 | 5,264 proteins |
