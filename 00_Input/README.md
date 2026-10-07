# 00_Input

The study data and the reference files the pipeline reads. No code lives here, no script writes
here, and no script downloads anything.

## The files

| File | One row is | Read by |
|---|---|---|
| `HRvLR_raw.xlsx` | a protein: `uniprot_id`, `protein`, `gene`, `description`, `n_seq`, then one intensity column per MS run (48 runs) | `01_Processing/01_Protein_Filtering` |
| `metadata.csv` | an MS run (48 rows): `sample_id`, `subject`, `arm`, `timepoint`, `group` | `01_Processing/01_Protein_Filtering` |
| `phenotype.csv` | a subject (16 rows): fibre area, fibre counts, muscle area and strength before and after training, the composite and the training volume | the analysis folders (contrasts, trait associations, clustering, figures); never 01_Processing |
| `blood_contaminants.csv` | a curated blood protein (95 rows) with its `class` and `reason` | `01_Processing/01_Protein_Filtering` |
| `HPA_annotations_full.tsv` | a Human Protein Atlas gene (20,162 rows) | `01_Processing/01_Protein_Filtering` |
| `RBC_proteome_reference.tsv` | a red-cell protein (7,012 rows): `acc`, `gene`, `sources`, `n_sources` | `01_Processing/01_Protein_Filtering` |
| `msigdb_2026.1.Hs.rds` | a gene in a gene set: MSigDB 2026.1.Hs Hallmark, Reactome, KEGG Legacy and GO:BP | `01_Processing/04_Gene_Sets` |
| `goslim_generic.obo` | a GO Slim term (GO Consortium generic slim, GO release 2026-07-26) | `01_Processing/04_Gene_Sets` |
| `STRING/` | STRING v12 human: links at combined score 700 or more, aliases, protein info | `01_Processing/07_Module_Characterisation` |

Sixteen people gave a biopsy at T1 (baseline), T2 (after the training programme) and T3 (after an
acute bout at the end of it): 48 MS runs, one per biopsy. `sample_id` names the run (`S03_T1`),
`subject` the person with their arm (`HR_S03`), and `group` the arm-by-time cell (`HR_T1`), six
cells of eight runs. Three runs fail on depth in `01_Protein_Filtering`, so 45 biopsies go
forward.

## Where the reference files came from

| File | Source |
|---|---|
| `msigdb_2026.1.Hs.rds` | `msigdbr::msigdbr(db_species = "HS", collection, subcollection)` with msigdbr 26.1.0, saved 2026-09-23. Membership changes between releases, so the release is fixed here rather than read from the installed package. |
| `goslim_generic.obo` | <https://current.geneontology.org/ontology/subsets/goslim_generic.obo> |
| `STRING/` | <https://string-db.org/cgi/download?species_text=Homo+sapiens>, version 12.0: `protein.links` (the `min700` file), `protein.aliases`, `protein.info`. At score 700 STRINGdb reads only the `min700` links file. |

## HR and LR are the median split of the composite

`arm` splits the 16 subjects on `comp_hypertrophy`: the top eight are HR, the bottom eight LR. The
composite tracks fibre-area change far more than strength gain (its formula is not in the
repository), so a fibre-area difference between arms is built in and is not a finding.

## phenotype.csv holds pre and post side by side

Pre is the T1 biopsy and post the T2 biopsy; T3 has no phenotype. Units sit in each column name
(`fcsa_mixed_pre_um2`, `leg_press_1rm_post_kg`). `comp_hypertrophy` is the source's composite and
`volume_load_total_kg` the total kilograms lifted over the programme, one value per subject.
`leg_ext_1rm_post_kg` is missing for one subject.

The `fibres_*` columns are MyoVision fibre counts, labelled "Number of fCSA" in the source workbook.
`phenotype_tables()` in `R/helpers.R` leaves them out of every association.

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
