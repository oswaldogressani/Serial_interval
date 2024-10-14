#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 17: n=10; Weibull target; single interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S17 <- sisimul(n = 10, dist = "weibull", censor = "single", seedval = 1916)

# Extract coverage performance
S17coverage_n10 <- data.frame(CP90 = S17$Simperf[,4], CP95 = S17$Simperf[,5])
save(S17coverage_n10, file = "S17coverage_n10.RData")

# Print table for Latex input
print(xtable(S17$Simperf, digits = 3, label = "Scenario 17"), file="S17-table.txt")

# Print results in console
round(S17$Simperf,3)


