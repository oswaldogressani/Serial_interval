#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 7: n=50; Gaussian target; single interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S7 <- sisimul(n = 50, dist = "gaussian", censor = "single", seedval = 1906)

# Extract coverage performance
S7coverage_n50 <- data.frame(CP90 = S7$Simperf[,4], CP95 = S7$Simperf[,5])
save(S7coverage_n50, file = "S7coverage_n50.RData")

# Print table for Latex input
print(xtable(S7$Simperf, digits = 3, label = "Scenario 7"), file="S7-table.txt")

# Print results in console
round(S7$Simperf,3)


