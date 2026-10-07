# Muscle_strength

People clustered on their strength response (% change in leg-press and leg-extension 1RM), then the HR/LR template run on the clusters.

The gap statistic selects K = 1; K = 2 is imposed, with the reason in `01_Clustering/c_data/01_cluster.xlsx` (`choice`). K = 2 is stable (bootstrap Jaccard 0.84 and 0.88) and cuts across HR/LR. HR_S22 has no leg-extension 1RM and is unassigned, so 42 runs remain.

| Substage | Copy of |
|---|---|
| [`01_Clustering`](01_Clustering/README.md) | new: K diagnostics, assignments, relabelled inputs |
| [`02_Design_Contrasts`](02_Design_Contrasts/README.md) | `02_HRvLR_groups/01_Design_Contrasts` |
| [`03_Protein_Contrasts`](03_Protein_Contrasts/README.md) | `02_HRvLR_groups/02_Protein_Contrasts` |
| [`04_Set_Contrasts`](04_Set_Contrasts/README.md) | `02_HRvLR_groups/03_Set_Contrasts` |
| [`05_Module_Contrasts`](05_Module_Contrasts/README.md) | `02_HRvLR_groups/04_Module_Contrasts` |
| [`06_Classification/01_Proteins`](06_Classification/01_Proteins/README.md) | `02_HRvLR_groups/06_Classification/01_Proteins` |
| [`06_Classification/02_Sets`](06_Classification/02_Sets/README.md) | `02_HRvLR_groups/06_Classification/02_Sets` |
| [`06_Classification/03_Modules`](06_Classification/03_Modules/README.md) | `02_HRvLR_groups/06_Classification/03_Modules` |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | `02_HRvLR_groups/07_Trait_Associations` |
