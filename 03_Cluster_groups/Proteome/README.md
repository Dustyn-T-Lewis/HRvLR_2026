# Proteome

People clustered on their proteome, with no phenotype or group label, separately at T1, T2, T3,
training change and acute change. Features: the trunk's module eigengenes, z-scored per input.

The gap statistic selects K = 1 at every input: the proteome does not split these people into
groups. K = 2 is imposed at every input for the exploratory comparisons, and most of those splits
are unstable (bootstrap Jaccard 0.56 to 0.86). None tracks HR/LR.

Every stage from `02_Design_Contrasts` to `07_Trait_Associations` is the HR/LR template run once per input on that
input's clusters, C1 and C2. Each input writes to its own folder (`03_Protein_Contrasts/T1/`, `03_Protein_Contrasts/training/`,
and so on), and `stack_inputs()` (`R/inputs.R`) joins the five workbooks of a stage into
`c_data/all_inputs.xlsx`, with an `input` column on every row and a `validity` column wherever a row
names its comparison.

Validity says how far a comparison can be trusted when the groups were built from the proteome:

| Tag | Rule | Read as |
|---|---|---|
| circular | the comparison uses exactly the clustering's time points | description: the clustering made the difference |
| caution | it shares a time point with the clustering | regression to the mean can make a difference |
| clean | it shares none | a test |

| Substage | Asks | Copy of |
|---|---|---|
| [`01_Clustering`](01_Clustering/README.md) | how people split at each input | new: K diagnostics, assignments, relabelled inputs per input |
| [`02_Design_Contrasts`](02_Design_Contrasts/README.md) | which comparisons | `02_HRvLR_groups/01_Design_Contrasts` |
| [`03_Protein_Contrasts`](03_Protein_Contrasts/README.md) | which proteins differ | `02_HRvLR_groups/02_Protein_Contrasts` |
| [`04_Set_Contrasts`](04_Set_Contrasts/README.md) | which gene sets differ | `02_HRvLR_groups/03_Set_Contrasts` |
| [`05_Module_Contrasts`](05_Module_Contrasts/README.md) | which modules differ | `02_HRvLR_groups/04_Module_Contrasts` |
| [`06_Classification/01_Proteins`](06_Classification/01_Proteins/README.md), [`02_Sets`](06_Classification/02_Sets/README.md), [`03_Modules`](06_Classification/03_Modules/README.md) | does one feature separate the clusters | `02_HRvLR_groups/06_Classification` |
| [`07_Trait_Associations`](07_Trait_Associations/README.md) | does a phenotype link differ by cluster | `02_HRvLR_groups/07_Trait_Associations` |
| [`08_Cluster_Phenotype`](08_Cluster_Phenotype/README.md) | do the clusters differ in size or strength | new: phenotype played no part in the clustering, so these are tests |

The methods, figures and sheets of each copy are described in the template's README.
