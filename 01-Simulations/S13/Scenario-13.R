#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 13: n=10; Weibull target; doubly interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S13 <- sisimul(n = 10, dist = "weibull", censor = "double", seedval = 1912)

# Extract coverage performance
S13coverage_n10 <- data.frame(CP90 = S13$Simperf[,4], CP95 = S13$Simperf[,5])
save(S13coverage_n10, file = "S13coverage_n10.RData")

# Print table for Latex input
print(xtable(S13$Simperf, digits = 3, label = "Scenario 13"), file="S13-table.txt")

# Print results in console
round(S13$Simperf,3)


