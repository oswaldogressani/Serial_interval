#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 9: n=10; Gaussian target; no censoring                        #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S9 <- sisimul(n = 10, dist = "gaussian", censor = "none", seedval = 1908)

# Extract coverage performance
S9coverage_n10 <- data.frame(CP90 = S9$Simperf[,4], CP95 = S9$Simperf[,5])
save(S9coverage_n10, file = "S9coverage_n10.RData")

# Print table for Latex input
print(xtable(S9$Simperf, digits = 3, label = "Scenario 9"), file="S9-table.txt")

# Print results in console
round(S9$Simperf,3)


