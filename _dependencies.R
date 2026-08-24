# Every package this project uses, named so renv can see them.
#
# This file is never sourced. It exists because renv's dependency scan reads
# library() and pkg:: calls and cannot see inside pacman::p_load(), so a
# snapshot taken without it silently drops most of the stack (rstudio/renv#143).
#
# ragg earns its place despite appearing in no script. ggplot2::ggsave() picks
# ragg::agg_png() when ragg is installed and grDevices::png() when it is not,
# and the two shape text differently, so a library without it re-renders every
# committed figure to different bytes without erroring.

library(AnnotationDbi)
library(GO.db)
library(MASS)
library(MsCoreUtils)
library(WGCNA)
library(broom)
library(digest)
library(dplyr)
library(fgsea)
library(forcats)
library(ggplot2)
library(ggrepel)
library(here)
library(imp4p)
library(imputeLCMD)
library(limma)
library(lme4)
library(lintr)
library(matrixStats)
library(mclust)
library(missForest)
library(msigdbr)
library(openxlsx)
library(org.Hs.eg.db)
library(pacman)
library(patchwork)
library(proteoDA)
library(purrr)
library(ragg)
library(readr)
library(readxl)
library(scales)
library(singscore)
library(stringr)
library(styler)
library(testthat)
library(withr)
library(tibble)
library(tidyr)
