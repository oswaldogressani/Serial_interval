#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 5: n=10; Gaussian target; single interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S5 <- sisimul(n = 10, dist = "gaussian", censor = "single", seedval = 1904)

# Extract coverage performance
S5coverage_n10 <- data.frame(CP90 = S5$Simperf[,4], CP95 = S5$Simperf[,5])
save(S5coverage_n10, file = "S5coverage_n10.RData")

# Print table for Latex input
print(xtable(S5$Simperf, digits = 3, label = "Scenario 5"), file="S5-table.txt")

# Print results in console
round(S5$Simperf,3)


