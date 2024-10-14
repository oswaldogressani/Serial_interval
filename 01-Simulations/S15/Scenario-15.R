#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 15: n=50; Weibull target; doubly interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S15 <- sisimul(n = 50, dist = "weibull", censor = "double", seedval = 1914)

# Extract coverage performance
S15coverage_n50 <- data.frame(CP90 = S15$Simperf[,4], CP95 = S15$Simperf[,5])
save(S15coverage_n50, file = "S15coverage_n50.RData")

# Print table for Latex input
print(xtable(S15$Simperf, digits = 3, label = "Scenario 15"), file="S15-table.txt")

# Print results in console
round(S15$Simperf,3)


