#------------------------------------------------------------------------------#
#                                                                              #
# Nonparametric analysis of Omicron and Delta variant of SARS-CoV-2            #
# Data is from Kremer et al. 2022.                                             #
# Oswaldo Gressani. All rights reserved.                                       #
#                                                                              #
#------------------------------------------------------------------------------# 

# Source: Kremer, C., Braeye, T., et al. (2022). Serial intervals for SARS-CoV-2
# omicron and delta variants, Belgium, November 19–December 31, 2021. 
# Emerging infectious diseases, 28(8), 1699

# Source of dataset:
# https://github.com/cecilekremer/serial_interval/blob/main/data_final.csv
# Last accessed: April 16, 2025 at 15h05 CET.

rm(list = ls())

library("lubridate") # To use dmy() function
library("EpiDelays")

# Read data
data <- read.csv("RawData/data_final.csv", header = T, sep = ";")

n_tot <- nrow(data) # Total sample size

# Calendar dates for symptom onset times related to Omicron variant
date_infector_Omicron <- dmy(data$date.onsetCC[data$Variant=="Omicron"])
date_infectee_Omicron <- dmy(data$date.onset.as.index[data$Variant=="Omicron"])
prop_Omicron <- round(length(date_infector_Omicron)/n_tot * 100,2)
n_Omicron <- length(date_infector_Omicron)

# Calendar dates for symptom onset times related to Delta variant
date_infector_Delta <- dmy(data$date.onsetCC[data$Variant=="Delta"])
date_infectee_Delta <- dmy(data$date.onset.as.index[data$Variant=="Delta"])
prop_Delta <- round(length(date_infector_Delta)/n_tot * 100,2)
n_Delta <- length(date_infector_Delta)

# Serial interval for Omicron
SI_Omicron <- as.numeric(date_infectee_Omicron-date_infector_Omicron)
SI_Delta <- as.numeric(date_infectee_Delta-date_infector_Delta)
# all(SI_Omicron == data$serialInterval[data$Variant=="Omicron"])
# all(SI_Delta == data$serialInterval[data$Variant=="Delta"])


par(mfrow = c(2,2), mai = c(1,1,1,1))

ftsize <- 1.6
ftsize2 <- 1.3

#----- Plot A and B (histograms of Omicron and Delta)
hist(SI_Omicron, xlim = c(-6,14), 
     xlab = "Serial interval (days)", 
     ylab = "No.cases", main = "",
     xaxt = "n", cex.lab = ftsize, cex.axis = ftsize)
axis(1, at = seq(-6,14,by = 2), cex.axis = ftsize)
title(main = "A (Omicron)",adj = 0, cex.main = 2)

hist(SI_Delta, xlim = c(-6,12), 
     xlab = "Serial interval (days)", 
     ylab = "No.cases", main = "",
     breaks = 20, xaxt = "n", cex.lab = ftsize, cex.axis = ftsize)
axis(1, at = seq(-6,12, by = 2), cex.axis = ftsize)
title(main = "B (Delta)",adj = 0, cex.main = 2)

#------------------------------------------------------------------------------#
#       Estimation with nonparametric approach for the Omicron variant         #
#------------------------------------------------------------------------------#

set.seed(2022)

# Within a day coarsening for infector and infectee
t_Omicron_infector_left <- as.numeric(date_infector_Omicron)
t_Omicron_infector_right <- t_Omicron_infector_left + 1
t_Omicron_infectee_left <- as.numeric(date_infectee_Omicron)
t_Omicron_infectee_right<- t_Omicron_infectee_left + 1
sl <- t_Omicron_infectee_left - t_Omicron_infector_right
sr <- t_Omicron_infectee_right - t_Omicron_infector_left
x <- data.frame(sl = sl, sr = sr)
SIfit <- estimSI(x = x, nboot = 2000)
# round(SIfit$npestim, 2)
n <- n_Omicron
SIfit_Omicron <- SIfit

#----- Plot C (Visualize estimates for Omicron)
sl <- x$sl
sr <- x$sr
sw <- sr - sl
Fhat <- function(s) (1 / n) * sum((s - sl) / sw * (s >= sl & s <= sr) + (s > sr))
sf <- seq(min(sl) - 1, max(sr) + 1, length = 1000)
txtoff <- 1.1
B <- SIfit$nboot

plot(sf, sapply(sf, Fhat), type = "l", lwd = 2, 
     xlab = "Serial interval (days)", 
     ylab = "Estimated cdf", xaxt = "n",
     cex.axis = ftsize, cex.lab = ftsize)
axis(1, at = seq((min(x$sl) - 1),((max(x$sr) + 1)), by = 2), cex.axis = ftsize)

# Add CIs for 5th, 25th, 50th and 95th quantiles

#5th quantile
lines(c(SIfit$npestim$q0.05[5], SIfit$npestim$q0.05[6]), c(0.05,0.05),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.05[1], 0.05, type = "p", pch = 16)
text(x = SIfit$npestim$q0.05[6] + txtoff, y = 0.05, label = "5th quantile",
     cex = ftsize2)
title(main = "C", adj = 0, cex.main = 3)

#25th quantile
lines(c(SIfit$npestim$q0.25[5], SIfit$npestim$q0.25[6]), c(0.25,0.25),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.25[1], 0.25, type = "p", pch = 16)
text(x = SIfit$npestim$q0.25[6] + txtoff, y = 0.25, label = "25th quantile",
     cex = ftsize2)

#50th quantile (median)
lines(c(SIfit$npestim$q0.5[5], SIfit$npestim$q0.5[6]), c(0.50, 0.50),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.5[1], 0.5, type = "p", pch = 16)
text(x = SIfit$npestim$q0.5[6] + txtoff, y = 0.50, label = "50th quantile",
     cex = ftsize2)

#75th quantile
lines(c(SIfit$npestim$q0.75[5], SIfit$npestim$q0.75[6]), c(0.75,0.75),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.75[1], 0.75, type = "p", pch = 16)
text(x = SIfit$npestim$q0.75[5] - txtoff, y = 0.75, label = "75th quantile",
     cex = ftsize2)

#95th quantile
lines(c(SIfit$npestim$q0.95[5], SIfit$npestim$q0.95[6]), c(0.95,0.95),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.95[1], 0.95, type = "p", pch = 16)
text(x = SIfit$npestim$q0.95[5] - txtoff, y = 0.95, label = "95th quantile",
     cex = ftsize2)

#------------------------------------------------------------------------------#
#       Estimation with nonparametric approach for the Delta variant           #
#------------------------------------------------------------------------------#

# Within a day coarsening for infector and infectee
t_Delta_infector_left  <- as.numeric(date_infector_Delta)
t_Delta_infector_right <- t_Delta_infector_left + 1
t_Delta_infectee_left  <- as.numeric(date_infectee_Delta)
t_Delta_infectee_right <- t_Delta_infectee_left + 1
sl <- t_Delta_infectee_left - t_Delta_infector_right
sr <- t_Delta_infectee_right - t_Delta_infector_left
x <- data.frame(sl = sl, sr = sr)
SIfit <- estimSI(x = x, nboot = 2000)
# round(SIfit$npestim,2)
n <- n_Delta
SIfit_Delta <- SIfit

#----- Plot D (Visualize estimates for Delta)
sl <- x$sl
sr <- x$sr
sw <- sr - sl
Fhat <- function(s) (1 / n) * sum((s - sl) / sw * (s >= sl & s <= sr) + (s > sr))
sf <- seq(min(sl) - 1, max(sr) + 1, length = 1000)
txtoff <- 1.1
B <- SIfit$nboot

plot(sf, sapply(sf, Fhat), type = "l", lwd = 2, 
     xlab = "Serial interval (days)", 
     ylab = "Estimated cdf", xaxt = "n",
     cex.axis = ftsize, cex.lab = ftsize)
axis(1, at = seq((min(x$sl) - 1),((max(x$sr) + 1)), by = 2), cex.axis = ftsize)

# Add CIs for 5th, 25th, 50th and 95th quantiles

#5th quantile
lines(c(SIfit$npestim$q0.05[5], SIfit$npestim$q0.05[6]), c(0.05,0.05),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.05[1], 0.05, type = "p", pch = 16)
text(x = SIfit$npestim$q0.05[6] + txtoff, y = 0.05, label = "5th quantile",
     cex = ftsize2)
title(main = "D", adj = 0, cex.main = 3)

#25th quantile
lines(c(SIfit$npestim$q0.25[5], SIfit$npestim$q0.25[6]), c(0.25,0.25),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.25[1], 0.25, type = "p", pch = 16)
text(x = SIfit$npestim$q0.25[6] + txtoff, y = 0.25, label = "25th quantile",
     cex = ftsize2)

#50th quantile (median)
lines(c(SIfit$npestim$q0.5[5], SIfit$npestim$q0.5[6]), c(0.50, 0.50),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.5[1], 0.5, type = "p", pch = 16)
text(x = SIfit$npestim$q0.5[6] + txtoff, y = 0.50, label = "50th quantile",
     cex = ftsize2)

#75th quantile
lines(c(SIfit$npestim$q0.75[5], SIfit$npestim$q0.75[6]), c(0.75,0.75),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.75[1], 0.75, type = "p", pch = 16)
text(x = SIfit$npestim$q0.75[5] - txtoff, y = 0.75, label = "75th quantile",
     cex = ftsize2)

#95th quantile
lines(c(SIfit$npestim$q0.95[5], SIfit$npestim$q0.95[6]), c(0.95,0.95),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.95[1], 0.95, type = "p", pch = 16)
text(x = SIfit$npestim$q0.95[5] - txtoff, y = 0.95, label = "95th quantile",
     cex = ftsize2)

# Print final results

round(SIfit_Omicron$npestim, 2)
round(SIfit_Delta$npestim, 2)












