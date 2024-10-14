#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 20: n=100; Weibull target; single interval-censored           #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S20 <- sisimul(n = 100, dist = "weibull", censor = "single", seedval = 1919)

# Extract coverage performance
S20coverage_n100 <- data.frame(CP90 = S20$Simperf[,4], CP95 = S20$Simperf[,5])
save(S20coverage_n100, file = "S20coverage_n100.RData")

# Print table for Latex input
print(xtable(S20$Simperf, digits = 3, label = "Scenario 20"), file="S20-table.txt")

# Print results in console
round(S20$Simperf,3)


