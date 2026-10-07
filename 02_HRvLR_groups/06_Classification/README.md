# 06_Classification

Can a single protein, gene-set score or module eigengene tell one group from the other?

The contrast steps ask whether a group's average moved. These steps ask about people: knowing one
feature's value for a person, how well can you say which group they are in? All three levels call
`classify_step()` in `R/classify.R`.

| Step | Feature | Value per biopsy |
|---|---|---|
| [`01_Proteins`](01_Proteins/README.md) | a protein | log2 abundance, blood-corrected complete matrix |
| [`02_Sets`](02_Sets/README.md) | a gene set | singscore |
| [`03_Modules`](03_Modules/README.md) | a module | eigengene |

Five tasks compare the groups: at T1, T2 and T3, and on each person's training (T2 - T1) and acute
(T3 - T2) change. The AUC is the Mann-Whitney W / (n1 n2): the chance a random person from the
first group has the higher value, 0.5 for no separation. Every AUC has DeLong's 95% interval;
the features drawn in the report also have a stratified bootstrap interval (2,000 resamples). A
perfectly separating feature gets [1, 1] from the bootstrap and no DeLong interval.
