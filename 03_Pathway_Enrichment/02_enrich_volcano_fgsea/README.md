# 02_enrich_volcano_fgsea

Protein volcanoes with collapse-surviving pathways ringed. Computes nothing.

| | |
|---|---|
| **Reads** | `set_tests.rds` |
| **Writes** | 8 volcanoes |

Point colour reads protein BH FDR; the ring reads set fgsea FDR. Red is up, blue down. Six FDR
panels (both interactions and the four within-arm contrasts) and two pi-ranked repeats for the
training contrasts. The floor is not drawn. No protein reaches BH < 0.05, so the FDR panels carry
no protein labels and every point is grey.

Leading-edge genes are translated to point labels before drawing. When two ringed sets reduce to
one display name, the collection is appended.
