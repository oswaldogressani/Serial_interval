#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 10: n=20; Gaussian target; no censoring                       #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S10 <- sisimul(n = 20, dist = "gaussian", censor = "none", seedval = 1909)

# Extract coverage performance
S10coverage_n20 <- data.frame(CP90 = S10$Simperf[,4], CP95 = S10$Simperf[,5])
save(S10coverage_n20, file = "S10coverage_n20.RData")

# Print table for Latex input
print(xtable(S10$Simperf, digits = 3, label = "Scenario 10"), file="S10-table.txt")

# Print results in console
round(S10$Simperf,3)


