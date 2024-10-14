#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 12: n=100; Gaussian target; no censoring                      #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S12 <- sisimul(n = 100, dist = "gaussian", censor = "none", seedval = 1911)

# Extract coverage performance
S12coverage_n100 <- data.frame(CP90 = S12$Simperf[,4], CP95 = S12$Simperf[,5])
save(S12coverage_n100, file = "S12coverage_n100.RData")

# Print table for Latex input
print(xtable(S12$Simperf, digits = 3, label = "Scenario 12"), file="S12-table.txt")

# Print results in console
round(S12$Simperf,3)


