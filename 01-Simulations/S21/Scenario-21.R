#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 21: n=10; Weibull target; no censoring                        #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S21 <- sisimul(n = 10, dist = "weibull", censor = "none", seedval = 1920)

# Extract coverage performance
S21coverage_n10 <- data.frame(CP90 = S21$Simperf[,4], CP95 = S21$Simperf[,5])
save(S21coverage_n10, file = "S21coverage_n10.RData")

# Print table for Latex input
print(xtable(S21$Simperf, digits = 3, label = "Scenario 21"), file="S21-table.txt")

# Print results in console
round(S21$Simperf,3)


