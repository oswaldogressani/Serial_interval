# 2019-nCov dataset of Li et al. (2019).
# Estimation of serial interval with EpiLPS.
# Oswaldo Gressani. All rights reserved.

rm(list = ls())
library("EpiLPS")

# Load data
data <- read.table("SIdatnCoV2019.txt")

# Nonparametric SI estimation
set.seed(2019)
SIfit <- estimSI(data)
round(SIfit$estim,1)

