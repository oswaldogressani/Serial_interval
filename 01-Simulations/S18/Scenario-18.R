#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 18: n=20; Weibull target; single interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S18 <- sisimul(n = 20, dist = "weibull", censor = "single", seedval = 1917)

# Extract coverage performance
S18coverage_n20 <- data.frame(CP90 = S18$Simperf[,4], CP95 = S18$Simperf[,5])
save(S18coverage_n20, file = "S18coverage_n20.RData")

# Print table for Latex input
print(xtable(S18$Simperf, digits = 3, label = "Scenario 18"), file="S18-table.txt")

# Print results in console
round(S18$Simperf,3)


