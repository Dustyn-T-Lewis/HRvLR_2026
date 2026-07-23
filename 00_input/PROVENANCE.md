# Red-cell proteome reference — provenance

`RBC_proteome_reference.tsv` is the mature-erythrocyte protein list that Stage 01 uses to confirm
blood-contamination calls. Erythrocytes are enucleate, so no transcriptomic atlas (including HPA's
single-cell erythrocyte RNA) represents their proteome; the membrane skeleton — band 3 (`SLC4A1`),
spectrins, protein 4.2 — is invisible to RNA and only recognisable from red-cell mass spectrometry.
The reference collapses five published sources into one lookup so a red-cell removal rests on
independent evidence rather than a single dataset.

## Sources

`RBC_proteins.xlsx` (in this directory) holds one sheet per source. `a_script/build_blood_list.R`
parses them and writes `RBC_proteome_reference.tsv`; `sources` records which sources list each
protein and `n_sources` how many agree.

| Sheet | Source | Contributes |
|---|---|---|
| RESPIRE | Téletchéa et al. 2019, *PLOS ONE* 14(2):e0211043 (PMID 30794542). Repository of Enriched Structures of Proteins Involved in the Red Blood Cell Environment. | 736 curated red-cell proteins |
| Uniprot | UniProt "Erythrocyte" annotation query (UniProtKB, retrieved 2026-07-23). | 597 erythrocyte-annotated entries |
| JPR2017 | Bryk & Wiśniewski 2017, *J Proteome Res* 16(8):2752–2761 (doi:10.1021/acs.jproteome.7b00025). Quantitative RBC proteome with copy numbers. | 2,577 quantified proteins |
| CB2019 | Ravenhill et al. 2019, *Commun Biol* 2:350 (doi:10.1038/s42003-019-0596-y). Erythrocyte surface proteome by plasma-membrane profiling. | 1,534 surface/membrane proteins |
| This study | Bai et al. 2026, *Sci Data* 13 (doi:10.1038/s41597-026-06792-5); ProteomeXchange PXD067677. Mature-RBC membrane + cytoplasm proteome. | 5,264 proteins, largest to date |

Union: 8,011 entries (6,452 with accession, 7,937 with gene); 2,089 backed by two or more sources.

## Why membership never removes on its own

A deep red-cell proteome shares roughly 70% of any skeletal-muscle proteome — glycolysis, tubulins,
ferritins, chaperones — because both cell types carry the same cytosolic housekeeping program.
Membership therefore cannot flag contamination by itself; it only corroborates the in-sample
haemoglobin-covariation call (`blood_cor`, computed in Stage 01). A protein is removed as red-cell
only when it both tracks the haemoglobin index in these samples and appears in this reference.

## Regenerating

```r
Rscript 00_input/a_script/02_build_blood_list.R
```

Reads `RBC_proteins.xlsx`, rewrites `RBC_proteome_reference.tsv`. Rerun only when a source updates.
