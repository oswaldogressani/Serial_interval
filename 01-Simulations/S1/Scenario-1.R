#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 1: n=10; Gaussian target; doubly interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")
library("EpiLPS")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S1 <- sisimul(n = 10, dist = "gaussian", censor = "double", seedval = 1900)

# Extract coverage performance
S1coverage_n10 <- data.frame(CP90 = S1$Simperf[,4], CP95 = S1$Simperf[,5])
save(S1coverage_n10, file = "S1coverage_n10.RData")

# Print table for Latex input
print(xtable(S1$Simperf, digits = 3, label = "Scenario 1"), file="S1-table.txt")

# Print results in console
round(S1$Simperf,3)


