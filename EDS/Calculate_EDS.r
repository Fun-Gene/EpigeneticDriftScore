# Calculate EDS

source("FUN_EDS_POS.r")
source("FUN_EDS_NEG.r")

input <- readRDS("input.rds")

eds_pos <- EDS_POS(input, ref_stat, ref_drift, ref_coef)
eds_neg <- EDS_NEG(input, ref_stat, ref_drift, ref_coef)

eds <- list(eds_pos, eds_neg)
