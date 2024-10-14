#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 14: n=20; Weibull target; doubly interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S14 <- sisimul(n = 20, dist = "weibull", censor = "double", seedval = 1913)

# Extract coverage performance
S14coverage_n20 <- data.frame(CP90 = S14$Simperf[,4], CP95 = S14$Simperf[,5])
save(S14coverage_n20, file = "S14coverage_n20.RData")

# Print table for Latex input
print(xtable(S14$Simperf, digits = 3, label = "Scenario 14"), file="S14-table.txt")

# Print results in console
round(S14$Simperf,3)


