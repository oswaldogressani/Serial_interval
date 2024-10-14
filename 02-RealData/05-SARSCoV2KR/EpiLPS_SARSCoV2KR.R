# Estimation of SI of Omicron and Delta variant of SARS-CoV-2 with EpiLPS
# Oswaldo Gressani. All rights reserved.

library("lubridate") # To use dmy() function
library("EpiLPS")

rm(list = ls())
source("estimSI_boot.R")

# Read data
data <- read.csv("data_final.csv", header = T, sep = ";")

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

pdf(file = "Fig4.pdf", width = 24.5, height = 14.5)

par(mfrow = c(2,2), mai = c(1,1,1,1))

#----- Plot A and B (histograms of Omicron and Delta)
hist(SI_Omicron, xlim = c(-6,14), col = "#54C8F0", 
     xlab = "Serial interval (days)", 
     ylab = "No.cases", main = "",
     xaxt = "n", cex.lab = 2, cex.axis = 2)
axis(1, at = seq(-6,14,by = 2),cex.axis = 2)
title(main = "A (Omicron)",adj = 0, cex.main = 2.5)

hist(SI_Delta, xlim = c(-6,12), col = "#00F0E1", 
     xlab = "Serial interval (days)", 
     ylab = "No.cases", main = "",
     breaks = 20, xaxt = "n", cex.lab = 2, cex.axis = 2)
axis(1, at = seq(-6,12, by = 2),cex.axis = 2)
title(main = "B (Delta)",adj = 0, cex.main = 2.5)

# Estimation with EpiLPS for data on the Omicron variant
xOmicron <- data.frame(sL = SI_Omicron - 0.5, sR = SI_Omicron + 0.5)

set.seed(2022)

fitOmicron <- estimSI_boot(x = xOmicron)

round(fitOmicron$estim, 2)

xO <- xOmicron
SIfitO <- fitOmicron

#----- Plot C (Visualize smoothed bootstraped cdfs and estimates for Omicron)
bsLO <- min(xO$sL) - 0.5
bsRO <- max(xO$sR) + 0.5
sfineO <- seq(bsLO, bsRO, length = 100)
dsfineO <- sfineO[2] - sfineO[1]
BO <- nrow(SIfitO$bootsamples)
sboot_cdfO <- matrix(0, nrow = BO, ncol = 100)

for(b in 1:BO){
  fsbootO <- histosmooth(SIfitO$bootsamples[b,], xl = bsLO, xr = bsRO, K = 12)
  sboot_densO <- sapply(sfineO, fsbootO$fdens)
  sboot_densO <- sboot_densO/sum(sboot_densO * dsfineO)
  sboot_cdfO[b,] <- cumsum(sboot_densO * dsfineO)
  print(b)
}

plot(sfineO, sapply(sfineO, SIfitO$Fhat), type = "l", col = "darkblue",
     xlab = "Serial interval (days)",
     ylab = "Cumulative distribution function",
     xaxt = 'n', xlim = c(min(xO$sL)-1, max(xO$sR)+1),
     cex.lab = 2, cex.axis = 2)
title(main = "C (Omicron)", adj = 0, cex.main = 2.5)
axis(1, at=seq(-6, 14, by = 2), cex.axis = 2)
for(b in 1:BO){
  lines(sfineO, sboot_cdfO[b,], type = "l", col = "#54C8F0")
}
lines(sfineO, sapply(sfineO, SIfitO$Fhat), type = "l", col = "darkblue",
      lwd = 3)

low_95 <- SIfitO$estim$CI95p_l[3:7] 
up_95 <- SIfitO$estim$CI95p_r[3:7]
perc <- c(0.05,0.25,0.50,0.75,0.95)
perctxt <- paste0(perc * 100, "th percentile")
for(j in 1:length(up_95)){
  lines(x = c(low_95[j],up_95[j]), y = c(perc[j],perc[j]), type = "l",
        col = "darkblue", lwd = 3, lty = 1)
  lines(x = SIfitO$estim$Estim[3:7], perc, type = "p", pch = 16,
        col = "darkblue")
  text(x = up_95[j]+1.4, y = perc[j]-0.05, perctxt[j], offset = 1, cex = 1.5)
}
text(x = -1.8, y = 0.9, paste0("Mean SI (95% CI): ",
                               round(SIfitO$estim$Estim[1],2)," (",
                               round(SIfitO$estim$CI95p_l[1],2),
                               "-",
                               round(SIfitO$estim$CI95p_r[1],2),")"), cex = 1.5)
text(x = -2.1, y = 0.8, paste0("SD SI (95% CI): ",
                               round(SIfitO$estim$Estim[2],2)," (",
                               round(SIfitO$estim$CI95p_l[2],2),
                               "-",
                               round(SIfitO$estim$CI95p_r[2],2),")"), cex = 1.5)


# Estimation with EpiLPS for data on the Delta variant
xDelta <- data.frame(sL = SI_Delta - 0.5, sR = SI_Delta + 0.5)

fitDelta <- estimSI_boot(x = xDelta)

round(fitDelta$estim,2)

x <- xDelta
SIfit <- fitDelta

#----- Plot C (Visualize smoothed bootstraped cdfs and estimates for Omicron)
bsL <- min(x$sL) - 1.5
bsR <- max(x$sR) + 1.5
sfine <- seq(bsL, bsR, length = 100)
dsfine <- sfine[2] - sfine[1]
B <- nrow(SIfit$bootsamples)
sboot_cdf <- matrix(0, nrow = B, ncol = 100)

for(b in 1:B){
  fsboot <- histosmooth(SIfit$bootsamples[b,], xl = bsL, xr = bsR, K = 5)
  sboot_dens <- sapply(sfine, fsboot$fdens)
  sboot_dens <- sboot_dens/sum(sboot_dens * dsfine)
  sboot_cdf[b,] <- cumsum(sboot_dens * dsfine)
  print(b)
}

plot(sfine, sapply(sfine, SIfit$Fhat), type = "l", col = "darkgreen",
     xlab = "Serial interval (days)",
     ylab = "Cumulative distribution function",
     xaxt = 'n', xlim = c(min(x$sL)-1, max(x$sR)+1),
     cex.lab = 2, cex.axis = 2)
title(main = "D (Delta)", adj = 0, cex.main = 2.5)
axis(1, at=seq(-6, 14, by = 2), cex.axis = 2)
for(b in 1:B){
  lines(sfine, sboot_cdf[b,], type = "l", col = "#00F0E1")
}
lines(sfine, sapply(sfine, SIfit$Fhat), type = "l", col = "darkgreen",
      lwd = 3)

low_95 <- SIfit$estim$CI95p_l[3:7] 
up_95 <- SIfit$estim$CI95p_r[3:7]
perc <- c(0.05,0.25,0.50,0.75,0.95)
perctxt <- paste0(perc * 100, "th percentile")
for(j in 1:length(up_95)){
  lines(x = c(low_95[j],up_95[j]), y = c(perc[j],perc[j]), type = "l",
        col = "darkgreen", lwd = 3, lty = 1)
  lines(x = SIfit$estim$Estim[3:7], perc, type = "p", pch = 16,
        col = "darkgreen")
  text(x = up_95[j]+1, y = perc[j]-0.05, perctxt[j], offset = 1, cex = 1.5)
}
text(x = -2.5, y = 0.9, paste0("Mean SI (95% CI): ",
                            format(round(SIfit$estim$Estim[1],2), nsmall = 2),
                            " (",
                            round(SIfit$estim$CI95p_l[1],2),
                            "-",
                            round(SIfit$estim$CI95p_r[1],2),")"), cex = 1.5)
text(x = -2.5, y = 0.8, paste0("SD SI (95% CI): ",
                            round(SIfit$estim$Estim[2],2)," (",
                            format(round(SIfit$estim$CI95p_l[2],2),nsmall = 2),
                            "-",
                            round(SIfit$estim$CI95p_r[2],2),")"), cex = 1.5)

dev.off()












