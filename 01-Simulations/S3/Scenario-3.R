#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#       Scenario 3: n=50; Gaussian target; doubly interval-censored            #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

library("xtable")

rm(list = ls())
source("sisim.R")
source("sisimul.R")

# Run scenario
S3 <- sisimul(n = 50, dist = "gaussian", censor = "double", seedval = 1902)

# Extract coverage performance
S3coverage_n50 <- data.frame(CP90 = S3$Simperf[,4], CP95 = S3$Simperf[,5])
save(S3coverage_n50, file = "S3coverage_n50.RData")

# Print table for Latex input
print(xtable(S3$Simperf, digits = 3, label = "Scenario 3"), file="S3-table.txt")

# Print results in console
round(S3$Simperf,3)


