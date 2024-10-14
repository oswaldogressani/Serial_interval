#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 23: n=50; Weibull target; no censoring                        #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S23 <- sisimul(n = 50, dist = "weibull", censor = "none", seedval = 1922)

# Extract coverage performance
S23coverage_n50 <- data.frame(CP90 = S23$Simperf[,4], CP95 = S23$Simperf[,5])
save(S23coverage_n50, file = "S23coverage_n50.RData")

# Print table for Latex input
print(xtable(S23$Simperf, digits = 3, label = "Scenario 23"), file="S23-table.txt")

# Print results in console
round(S23$Simperf,3)


