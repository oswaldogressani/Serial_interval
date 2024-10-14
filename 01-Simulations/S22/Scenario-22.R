#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 22: n=20; Weibull target; no censoring                        #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S22 <- sisimul(n = 20, dist = "weibull", censor = "none", seedval = 1921)

# Extract coverage performance
S22coverage_n20 <- data.frame(CP90 = S22$Simperf[,4], CP95 = S22$Simperf[,5])
save(S22coverage_n20, file = "S22coverage_n20.RData")

# Print table for Latex input
print(xtable(S22$Simperf, digits = 3, label = "Scenario 22"), file="S22-table.txt")

# Print results in console
round(S22$Simperf,3)


