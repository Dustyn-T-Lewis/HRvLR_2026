# 03_Classify

Scores how well each module eigengene separates the groups on the protein level's eight tasks.

| | |
|---|---|
| Reads | `01_Build_Modules/c_data/modules.rds` |
| Writes | `c_data/03_classify.xlsx`, `b_reports/03_classify_figures.pdf` |
| Run | `Rscript 04_Network/03_Classify/a_script/03_classify.R` |
| Cost | about 5 s |

The tasks, AUC and tests are those of `02_Differential_Expression/03_Classify`, on twelve
eigengenes. BH within task. `chance_expectation` gives each task's group sizes, smallest attainable
p (HR training 0.031, HR acute 0.016, LR 0.0078), tested, nominal and BH counts; chance is 0.6
modules per task. The figure PDF opens on the chance ratio per task, then draws an ROC curve for
every module at p < 0.05 on each task. The floor task gets no curves.
