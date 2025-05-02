#------------------------------------------------------------------------------#
#                                                                              #
# Nonparametric analysis of illness onset data for 2019-nCoV in Wuhan, China   #
# from Li et al. (2020) (n=6 transmission pairs)                              #
# Oswaldo Gressani. All rights reserved.                                       #
#                                                                              #
#------------------------------------------------------------------------------# 


rm(list = ls())

library("EpiDelays")

# Load data
data <- read.table("RawData/Cov19-Li-2020.txt")
x <- data[,c(2,3)]

# Nonparametric SI estimation
set.seed(123)
SIfit <- estimSI(x, nboot = 2000)
round(SIfit$npestim, 1)









