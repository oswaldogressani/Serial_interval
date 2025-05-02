#------------------------------------------------------------------------------#
#                                                                              #
# Nonparametric analysis of illness onset data for COVID-19                    #
# from Nishiura et al. (2020) (n=28 transmission pairs)                        #
# Oswaldo Gressani. All rights reserved.                                       #
#                                                                              #
#------------------------------------------------------------------------------# 

rm(list = ls())

# Source: Nishiura, H., Linton, N. M., & Akhmetzhanov, A. R. (2020). Serial 
# interval of novel coronavirus (COVID-19) infections. International journal 
# of infectious diseases, 93, 284-286.

# Source of .csv file: 
# https://github.com/aakhmetz/COVID19SerialInterval/tree/master/data
# Accessed April 16, 2025 at 14h25 CET.

datain <- read.csv(file = "Nishiura-COVID19.csv")

n <- nrow(datain) # Total number of transmission pairs

# We use index 1 to refere to infector and index 2 for infectee

# EL and ER are left and right bound of symptom onset times for infector
tS1L <- datain$EL
tS1R <- datain$ER
all(tS1L < tS1R)

#SL and SR are left and right bound of symptom onset times for infectee
tS2L <- datain$SL
tS2R <- datain$SR
all(tS2L < tS2R)

# Compute serial interval bounds
sl <- tS2L - tS1R
sr <- tS2R - tS1L

# Data out
data2019nCoV <- data.frame(Pair_index = seq_len(n), sl = sl, sr = sr)

# Extract dataset in txt format
write.table(data2019nCoV, file = "Nishiura-COVID19.txt")