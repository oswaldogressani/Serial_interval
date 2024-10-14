#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 11: n=50; Gaussian target; no censoring                       #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S11 <- sisimul(n = 50, dist = "gaussian", censor = "none", seedval = 1910)

# Extract coverage performance
S11coverage_n50 <- data.frame(CP90 = S11$Simperf[,4], CP95 = S11$Simperf[,5])
save(S11coverage_n50, file = "S11coverage_n50.RData")

# Print table for Latex input
print(xtable(S11$Simperf, digits = 3, label = "Scenario 11"), file="S11-table.txt")

# Print results in console
round(S11$Simperf,3)


