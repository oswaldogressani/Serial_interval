#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 24: n=100; Weibull target; no censoring                       #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S24 <- sisimul(n = 100, dist = "weibull", censor = "none", seedval = 1923)

# Extract coverage performance
S24coverage_n100 <- data.frame(CP90 = S24$Simperf[,4], CP95 = S24$Simperf[,5])
save(S24coverage_n100, file = "S24coverage_n100.RData")

# Print table for Latex input
print(xtable(S24$Simperf, digits = 3, label = "Scenario 24"), file="S24-table.txt")

# Print results in console
round(S24$Simperf,3)


