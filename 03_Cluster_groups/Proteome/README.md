# Proteome

People clustered on their proteome, with no phenotype or group label, separately at T1, T2, T3,
training change and acute change. Features: the trunk's module eigengenes, z-scored per input.

The gap statistic selects K = 1 at every input: the proteome does not split these people into
groups. K = 2 is imposed at every input for the exploratory comparisons, and most of those splits
are unstable (bootstrap Jaccard 0.56 to 0.86). None tracks HR/LR.

| Substage | Asks | Valid as |
|---|---|---|
| [`01_Cluster`](01_Cluster/README.md) | how people split at each input | K diagnostics |
| [`02_Describe`](02_Describe/README.md) | what separates the two clusters | description: the clusters were built from these data |
| [`03_Phenotype`](03_Phenotype/README.md) | do the clusters differ in size or strength | test: phenotype played no part |
| [`04_Interactions`](04_Interactions/README.md) | do clusters respond differently later, or link differently to phenotype | test |
