#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 2: n=20; Gaussian target; doubly interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S2 <- sisimul(n = 20, dist = "gaussian", censor = "double", seedval = 1989)

# Extract coverage performance
S2coverage_n20 <- data.frame(CP90 = S2$Simperf[,4], CP95 = S2$Simperf[,5])
save(S2coverage_n20, file = "S2coverage_n20.RData")

# Print table for Latex input
print(xtable(S2$Simperf, digits = 3, label = "Scenario 2"), file="S2-table.txt")

# Print results in console
round(S2$Simperf,3)


