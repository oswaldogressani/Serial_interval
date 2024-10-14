#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 6: n=20; Gaussian target; single interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S6 <- sisimul(n = 20, dist = "gaussian", censor = "single", seedval = 1905)

# Extract coverage performance
S6coverage_n20 <- data.frame(CP90 = S6$Simperf[,4], CP95 = S6$Simperf[,5])
save(S6coverage_n20, file = "S6coverage_n20.RData")

# Print table for Latex input
print(xtable(S6$Simperf, digits = 3, label = "Scenario 6"), file="S6-table.txt")

# Print results in console
round(S6$Simperf,3)


