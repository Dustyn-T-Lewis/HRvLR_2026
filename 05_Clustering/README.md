# 05_Clustering

Group the pi-gated proteome with three engines, summarise each group by an eigengene,
and link eigengenes to phenotype with a mixed model. Pure R (no Python). Only WGCNA
carries the inferential claim; the rest corroborate or are quarantined.

Design and rationale: `docs/design/specs/2026-06-23-hrvlr-clustering-phenotype-design.md`.

| Dir | Engine | Role |
|---|---|---|
| `a_wgcna_paired/` | WGCNA, paired design (within-subject centered) | primary — inferential |
| `b_mfuzz_gap/` | Mfuzz on HR−LR gap trajectory | companion — temporal shape |
| `c_supervised/` | sPLS / DIABLO (mixOmics) | quarantine — exploratory only |
| `d_integration/` | shared downstream | LMM, permutation null, BH, concordance, ORA |

Each engine writes the same two files to its `c_data/`, the only interface
`d_integration/` depends on:

- `membership.csv` — `protein_id, group_id, membership_weight`
- `eigengene.csv` — `sample_id, subject, group_arm, timepoint, group_id, ME`

Within each engine dir: `a_script/` code, `b_reports/` figures/tables, `c_data/`
the two artifacts above.
