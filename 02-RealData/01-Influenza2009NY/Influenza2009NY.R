#------------------------------------------------------------------------------#
#                                                                              #
# Nonparametric analysis of the 2009 H1N1 influenza outbreak in New York City  #
# Oswaldo Gressani. All rights reserved.                                       #
#                                                                              #
#------------------------------------------------------------------------------#  


rm(list = ls())

library("EpiDelays")

# Load data
data <- read.table("RawData/H1N1-Lessler-2009.txt")
x <- data[,c(2,3)]


# Nonparametric SI estimation
set.seed(123)
SIfit <- estimSI(x, nboot = 2000)
round(SIfit$npestim, 1)

#----------------- Plot results

pdf(file = "Fig4.pdf", width = 24, height = 8.5)

par(mfrow = c(1,2), mai = c(1,1,1,1))

#----- Plot A (Observed serial interval windows)
ftsize <- 1.6
ftsize2 <- 1.3
n <- nrow(x)
pairidx <- rep(seq_len(n), each = 2)
siwindowdat <- cbind(matrix(pairidx, ncol = 2, byrow = T), cbind(x$sl, x$sr))

plot(siwindowdat[1,3:4], siwindowdat[1,1:2], type = "l", 
     xlim = c(min(x$sl) - 1, max(x$sr) + 1), ylim = c(1, n),
     lwd = 2, col = "black", 
     xlab = "Serial interval window", ylab = "Pair index number",
     xaxt = 'n', yaxt = "n", cex.axis = ftsize, cex.lab = ftsize)
abline(v = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1),
       col = "gray85", lty = 2)
for(i in 1:n){
  lines(siwindowdat[i,3:4], siwindowdat[i,1:2], type = "l",
        lwd = 2, col = "black")
}
axis(1, at = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1), cex.axis = ftsize)
axis(2, at = seq(2, n, by = 2), cex.axis = ftsize)
title(main = "A", adj = 0, cex.main = 3)


#----- Plot B (CDF for point estimates and bootstrapped CIs)
sl <- x$sl
sr <- x$sr
sw <- sr - sl
Fhat <- function(s) (1 / n) * sum((s - sl) / sw * (s >= sl & s <= sr) + (s > sr))
sf <- seq(min(sl) - 1, max(sr) + 1, length = 1000)
txtoff <- 0.65
B <- SIfit$nboot

plot(sf, sapply(sf, Fhat), type = "l", lwd = 2, 
     xlab = "Serial interval (days)", 
     ylab = "Estimated cdf", xaxt = "n",
     cex.axis = ftsize, cex.lab = ftsize)
axis(1, at = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1), cex.axis = ftsize)

# Add CIs for 5th, 25th, 50th and 95th quantiles

#5th quantile
lines(c(SIfit$npestim$q0.05[5], SIfit$npestim$q0.05[6]), c(0.05,0.05),
      type = "l", col = "black", lwd = 2)
lines(SIfit$npestim$q0.05[1], 0.05, type = "p", pch = 16)
text(x = SIfit$npestim$q0.05[6] + txtoff, y = 0.05, label = "5th quantile",
     cex = ftsize2)
title(main = "B", adj = 0, cex.main = 3)

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


dev.off()

# Results by Lessler et al. (2009)

wshape <- 2.36
wscale <- 3.18

mean_weibull <- wscale * gamma(1 + 1 / wshape)
sd_weibull <- sqrt(wscale^2 * (gamma(1+2/wshape) - (gamma(1+1/wshape))^2))

round(mean_weibull, 1)
round(sd_weibull, 1)

round(qweibull(p = c(0.50, 0.95), shape = wshape, scale = wscale), 1)


















