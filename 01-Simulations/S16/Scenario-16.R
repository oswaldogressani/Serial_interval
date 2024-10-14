#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 16: n=100; Weibull target; doubly interval-censored           #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S16 <- sisimul(n = 100, dist = "weibull", censor = "double", seedval = 1915)

# Extract coverage performance
S16coverage_n100 <- data.frame(CP90 = S16$Simperf[,4], CP95 = S16$Simperf[,5])
save(S16coverage_n100, file = "S16coverage_n100.RData")

# Print table for Latex input
print(xtable(S16$Simperf, digits = 3, label = "Scenario 16"), file="S16-table.txt")

# Print results in console
round(S16$Simperf,3)


