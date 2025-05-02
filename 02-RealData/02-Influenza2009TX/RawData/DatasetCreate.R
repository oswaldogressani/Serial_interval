#------------------------------------------------------------------------------#
#                                                                              #
# Nonparametric analysis of the 2009 H1N1 influenza pandemic in                #
# San Antonio, Texas, USA                                                      #
# Oswaldo Gressani. All rights reserved.                                       #
#                                                                              #
#------------------------------------------------------------------------------# 


rm(list = ls())
library("EpiEstim")

data(Flu2009)
datain <- data.frame(Caseid = seq_len(16), Flu2009$si_data[,-5])
n <- nrow(datain) # Sample size

## EL/ER show lower/upper bound of the symptoms onset time of infector
## SL/SR show lower/upper bound of the symptoms onset time of infectee

# Infector
tS1L <- datain$EL
tS1R <- datain$ER
all(tS1L < tS1R)

# Infectee
tS2L <- datain$SL
tS2R <- datain$SR
all(tS2L < tS2R)

sl <- tS2L - tS1R
sr <- tS2R - tS1L

# Data out
dataH1N1 <- data.frame(Pair_index = seq_len(n), sl = sl, sr = sr)
write.table(dataH1N1, file = "H1N1-EpiEstim-2009.txt")














