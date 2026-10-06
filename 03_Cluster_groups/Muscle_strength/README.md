# Muscle_strength

People clustered on their strength response (% change in leg-press and leg-extension 1RM), then the HR/LR template run on the clusters.

The gap statistic selects K = 1; K = 2 is imposed, with the reason in `01_Cluster/c_data/01_cluster.xlsx` (`choice`). K = 2 is stable (bootstrap Jaccard 0.84 and 0.88) and cuts across HR/LR. HR_S22 has no leg-extension 1RM and is unassigned, so 42 runs remain.

| Substage | Copy of |
|---|---|
| [`01_Cluster`](01_Cluster/README.md) | new: K diagnostics, assignments, relabelled inputs |
| [`02_Design`](02_Design/README.md) | `02_HRvLR_groups/01_Design` |
| [`03_DEP`](03_DEP/README.md) | `02_HRvLR_groups/02_DEP` |
| [`04_Pathways`](04_Pathways/README.md) | `02_HRvLR_groups/03_Pathways` |
| [`05_Modules`](05_Modules/README.md) | `02_HRvLR_groups/04_Modules/01_Contrasts` |
| [`06_Classify/01_Proteins`](06_Classify/01_Proteins/README.md) | `02_HRvLR_groups/05_Classify/01_Proteins` |
| [`06_Classify/02_Sets`](06_Classify/02_Sets/README.md) | `02_HRvLR_groups/05_Classify/02_Sets` |
| [`06_Classify/03_Modules`](06_Classify/03_Modules/README.md) | `02_HRvLR_groups/05_Classify/03_Modules` |
| [`07_Associate`](07_Associate/README.md) | `02_HRvLR_groups/06_Associate` |
