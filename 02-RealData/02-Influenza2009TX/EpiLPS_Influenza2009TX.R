# H1N1 symptom onset data of EpiEstim package of Cori et al. (2013)
# Estimation of serial interval with EpiLPS.
# Oswaldo Gressani. All rights reserved.


rm(list = ls())
source("estimSI_boot.R")
library("EpiLPS")

# Load data
data <- read.table("dataH1N1Texas.txt")
x <- data[, c(2,3)]

# Nonparametric SI estimation
set.seed(123)
SIfit <- estimSI_boot(x)
round(SIfit$estim, 1)

#----------------- Fit visualization

pdf(file = "Fig2.pdf", width = 24, height = 8.5)

par(mfrow = c(1,2), mai = c(1,1,1,1))

#----- Plot A (Visualize serial interval window)
n <- nrow(x)
pairidx <- rep(seq_len(n), each = 2)
siwindowdat <- cbind(matrix(pairidx, ncol = 2, byrow = T), as.matrix(x))

plot(siwindowdat[1,1:2], siwindowdat[1,3:4], type = "l", 
     ylim = c(min(x$SL) - 1, max(x$SR) + 1),
     lwd = 2, col = "black", 
     xlab = "Pair index number",
     ylab = "Serial interval windows",
     xlim = c(1,n),
     xaxt = 'n', cex.axis = 2, cex.lab = 2)
for(i in 2:n){
  lines(siwindowdat[i,1:2],siwindowdat[i,3:4], type = "l",
        lwd = 3, col = "black")
}
axis(1, at = seq(2,16,by=2), cex.axis = 2)
title(main = "A",adj = 0, cex.main = 3)


#----- Plot B (Visualize smoothed bootstraped cdfs and estimates)
bsL <- min(x$SL) - 1.5
bsR <- max(x$SR) + 1.5
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
     xaxt = 'n', xlim = c(min(x$SL)-1, max(x$SR)+1),
     cex.lab = 2, cex.axis = 2)
title(main = "B", adj = 0, cex.main = 3)
axis(1, at=seq(min(x$SL)-1, max(x$SR)+1, by = 1), cex.axis = 2)
for(b in 1:B){
  lines(sfine, sboot_cdf[b,], type = "l", col = "#1FF0C4")
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
  text(x = up_95[j], y = perc[j]-0.05, perctxt[j], offset = 1, cex = 1.5)
}
text(x = 0.75, y = 0.9, paste0("Mean SI (95% CI): ",
                               format(round(SIfit$estim$Estim[1],1),nsmall=1),
                               " (",
                               round(SIfit$estim$CI95p_l[1],1),
                               "-",
                               format(round(SIfit$estim$CI95p_r[1],1), nsmall = 1),")")
                               , cex = 1.5)
text(x = 0.90, y = 0.8, paste0("SD SI (95% CI): ",
                               round(SIfit$estim$Estim[2],1)," (",
                               round(SIfit$estim$CI95p_l[2],1),
                               "-",
                               round(SIfit$estim$CI95p_r[2],1),")"), cex = 1.5)

dev.off()