#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 8: n=100; Gaussian target; single interval-censored           #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S8 <- sisimul(n = 100, dist = "gaussian", censor = "single", seedval = 1907)

# Extract coverage performance
S8coverage_n100 <- data.frame(CP90 = S8$Simperf[,4], CP95 = S8$Simperf[,5])
save(S8coverage_n100, file = "S8coverage_n100.RData")

# Print table for Latex input
print(xtable(S8$Simperf, digits = 3, label = "Scenario 8"), file="S8-table.txt")

# Print results in console
round(S8$Simperf,3)


