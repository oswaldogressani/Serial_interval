#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 19: n=50; Weibull target; single interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S19 <- sisimul(n = 50, dist = "weibull", censor = "single", seedval = 1918)

# Extract coverage performance
S19coverage_n50 <- data.frame(CP90 = S19$Simperf[,4], CP95 = S19$Simperf[,5])
save(S19coverage_n50, file = "S19coverage_n50.RData")

# Print table for Latex input
print(xtable(S19$Simperf, digits = 3, label = "Scenario 19"), file="S19-table.txt")

# Print results in console
round(S19$Simperf,3)


