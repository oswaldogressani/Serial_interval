#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 4: n=100; Gaussian target; doubly interval-censored           #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S4 <- sisimul(n = 100, dist = "gaussian", censor = "double", seedval = 1903)

# Extract coverage performance
S4coverage_n100 <- data.frame(CP90 = S4$Simperf[,4], CP95 = S4$Simperf[,5])
save(S4coverage_n100, file = "S4coverage_n100.RData")

# Print table for Latex input
print(xtable(S4$Simperf, digits = 3, label = "Scenario 4"), file="S4-table.txt")

# Print results in console
round(S4$Simperf,3)


