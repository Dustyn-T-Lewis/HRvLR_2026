# Gene-set references

Fixed reference snapshots read by `../a_script/gene_sets.qmd`. MSigDB supplies Hallmark,
Reactome, KEGG Legacy and GO:BP membership. GO Slim supplies broad ontology terms; the
script combines these with `GO.db` and MSigDB to construct additional gene sets.

| File | Source |
|---|---|
| `msigdb_2026.1.Hs.rds` | `msigdbr::msigdbr(db_species = "HS", collection, subcollection)` with msigdbr 26.1.0, saved 2026-09-23. Membership changes between releases, so the release is fixed here rather than read from the installed package. |
| `goslim_generic.obo` | <https://current.geneontology.org/ontology/subsets/goslim_generic.obo> |

The GO Slim file records GO release 2026-07-26. Preserve these snapshots when reproducing
the analysis; the current GO download URL can change over time.
