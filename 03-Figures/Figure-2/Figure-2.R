#------------------------------------------------------------------------------#
#                               Figure 2                                       #
#------------------------------------------------------------------------------#

rm(list = ls())
library("EpiDelays")

ftsize <- 1.8
n <- 15

pdf(file = "Figure 2.pdf", width = 17, height = 8)

par(mfrow = c(1,2))

# A ----------------------------------------------------
set.seed(1879)
x <- simSI(muS = 2.5, sdS = 3, maxcoarse = 3, 
           probcoarse = c(0.8,0.15,0.05), n = n)

x$sw
summary(x$sw)

# Plot generated serial intervals
n <- nrow(x)
pairidx <- rep(seq_len(n), each = 2)
siwindowdat <- cbind(matrix(pairidx, ncol = 2, byrow = T), cbind(x$sl, x$sr))

plot(siwindowdat[1,3:4], siwindowdat[1,1:2], type = "l", 
     xlim = c(min(x$sl) - 1, max(x$sr) + 1), ylim = c(1, (n + 2)),
     lwd = 2, col = "black", 
     xlab = "Serial interval window", ylab = "",
     xaxt = 'n', yaxt = "n", cex.axis = ftsize, cex.lab = ftsize)
abline(v = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1),
       col = "gray85", lty = 2)
for(i in 1:n){
  lines(siwindowdat[i,3:4], siwindowdat[i,1:2], type = "l",
        lwd = 2, col = "black")
  lines(x$s[i], i, type = "p", pch = 15, col = "black")
  # lines((0.5 * (siwindowdat[i,3] + siwindowdat[i,4])),
  #       i, type = "p", pch = 16, col = "black", lwd = 1)
}
axis(1, at = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1), cex.axis = ftsize)
legend("topleft", c("True serial interval"),
       pch = c(15), col = "black", bty = "n", cex = ftsize)

# B ----------------------------------------------------
set.seed(1919)
x <- simSI(muS = 2.5, sdS = 3, maxcoarse = 3, 
           probcoarse = c(0.8,0.15,0.05), n = n)

x$sw
summary(x$sw)

# Plot generated serial intervals
n <- nrow(x)
pairidx <- rep(seq_len(n), each = 2)
siwindowdat <- cbind(matrix(pairidx, ncol = 2, byrow = T), cbind(x$sl, x$sr))

plot(siwindowdat[1,3:4], siwindowdat[1,1:2], type = "l", 
     xlim = c(min(x$sl) - 1, max(x$sr) + 1), ylim = c(1, (n + 2)),
     lwd = 2, col = "black", 
     xlab = "Serial interval window", ylab = "",
     xaxt = 'n', yaxt = "n", cex.axis = ftsize, cex.lab = ftsize)
abline(v = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1),
       col = "gray85", lty = 2)
for(i in 1:n){
  lines(siwindowdat[i,3:4], siwindowdat[i,1:2], type = "l",
        lwd = 2, col = "black")
  lines(x$s[i], i, type = "p", pch = 15, col = "black")
  # lines((0.5 * (siwindowdat[i,3] + siwindowdat[i,4])),
  #       i, type = "p", pch = 16, col = "black", lwd = 2)
}
axis(1, at = seq((min(x$sl) - 1),((max(x$sr) + 1)),by=1), cex.axis = ftsize)
legend("topleft", c("True serial interval"),
       pch = c(15), col = "black", bty = "n", cex = ftsize)

dev.off()








