# 01_Associate / 02_Sets

Ties each gene set's singscore to phenotype, every subject pooled, no group label. The code is the shared
`associate_step()` in `R/associate.R`, which the protein, set and module levels all call.

| | |
|---|---|
| Reads | `01_Trunk/05_Scores/c_data/singscore.rds`, `01_Trunk/04_Gene_Sets/c_data/gene_sets.rds`, `01_Trunk/03_Quantify/c_data/DAList_corrected.rds` (sample sheet), `00_Input/phenotype.csv` |
| Writes | `c_data/02_sets.xlsx`, `b_reports/02_sets_figures.pdf` |
| Run | `Rscript 04_Continuous_no_group/01_Associate/02_Sets/a_script/associate_sets.R` |

## Method

Traits: fibre area (mixed, type I, type II), mCSA and the leg-press and leg-extension 1RMs, as
levels at each biopsy and as % change; plus `comp_hypertrophy` and `volume_load_total_kg`.

Sample model. Per trait, limma fits `value ~ timepoint + between + within` on the T1 and T2 biopsies
with subject blocked: `between` is the subject's mean trait, `within` the biopsy's deviation from
it. BH within trait and term.

Change score. One value per subject in six pairings: training change, T1 level and acute change
against the % changes, composite and volume (`baseline` is the T1 level against the change to
come); T1 and T2 levels against the trait levels at the same biopsy (`baseline_level`,
`post_level`); and the T3 level against the T2 trait levels (`post_bout_level`), since T3 follows
the same training and also carries the bout's acute effect. Per set: limma on one value per
subject (`p`), Spearman (`rho`, `p_rank`, from `cor.test`) and biweight midcorrelation (`bicor`,
WGCNA). BH and Storey q within pairing and outcome.

Two predictor: `value ~ comp_hypertrophy + strength gain`, once per 1RM. Interaction:
`value ~ growth * strength gain`: does the link to strength gain depend on how much the person
grew? A negative term means the strength slope falls as growth rises.

Each collection is its own screen: BH runs within collection.

## Figures

Overview pages: within- and between-person t for every trait beside status maps, then pooled rho
and limma status for every outcome, one page per pairing. Then every model hit, change-score pair
and two-predictor hit at p < 0.05, twelve panels to a page, ordered by p.

## Workbook

| Sheet | Holds |
|---|---|
| `counts` | Tested, p < 0.05 and BH < 0.05 per analysis and comparison. |
| `model` | Sample model: effect per trait unit, t, p, BH, Storey q per set, trait and term. |
| `change_score` | limma effect, t, p, BH, Storey q; Spearman rho and p_rank; bicor and its p. |
| `two_predictor` | Growth and strength terms, t, p, BH, Storey q. |
| `interaction` | Growth x strength term, t, p, BH, Storey q. |
| `plot_index` | PDF page of each overview, and page and panel of every p < 0.05 row. |
