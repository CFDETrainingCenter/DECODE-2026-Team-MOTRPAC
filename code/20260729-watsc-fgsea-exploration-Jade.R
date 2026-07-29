# Goal ----
# Examine the code used for the sexual-dimorphism analysis in the WAT-SC
# Article: https://www.nature.com/articles/s42255-023-00959-9
# GitHub: https://github.com/MoTrPAC/MotrpacRatTraining6moWATData/
# Try to replicate FGSEA --- output: transcriptomic part of Fig 2B, 2C


# From examining the GitHub repo, it looks like FGSEA involves the following files
# (1) `data-raw/msigdb_pathways.R` which generates `data/MSIGDB_PAYWAYS.rda`
#     it also uses the function `msigdbr2` from MotrpacRatTraining6moWAT
# (2) `vignettes/MAIN_FGSEA.Rmd` which uses `fgsea2` from MotrpacRatTraining6moWAT
# (3) `vignettes/artciles/plot_FGSEA_male_vs_female.Rmd` which creates Fig 2B and 2C
#     along with proteomics data



# Load libraries -----
library(tidyverse)
library(here)


# Import data ----
## DA results
read_csv(here("data_processed", "watsc_transcriptomics_da.csv"))


