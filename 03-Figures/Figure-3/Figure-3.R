#------------------------------------------------------------------------------#
#                   Asymptotic bias of nonparametric estimates                 #
#------------------------------------------------------------------------------#
#   -->  Gaussian SI with mean=2.8 and sd=2.5 from Kremer et al. (2022)        #
#   -->  Average coarseness of 2 days and epsilon-coarseness                   #
#------------------------------------------------------------------------------#

rm(list = ls())

library("EpiDelays")

#---- Tune scenario setting here
nvec <- seq(6, 500, by = 2)
nvec_len <- length(nvec)
S <- 50
meanSI  <-  2.8   # Kremer et al. (2022)
sdSI    <-  2.5   # Kremer et al. (2022)
qtrue   <- qnorm(p = c(0.05,0.25,0.50,0.75,0.95), mean = meanSI, sd = sdSI)
nboot <- 2
eps <- 0.01
set.seed(2025)
#-------------------------------

#--- Hosting average estimates for S = 50 replications with coarseness >= 2 days
mu_n <- c()
sd_n <- c()
qp05_n <- c()
qp25_n <- c()
qp50_n <- c()
qp75_n <- c()
qp95_n <- c()

#--- Hosting average estimates for S = 50 replications with epsilon-coarseness 
mu_n_ec <- c()
sd_n_ec <- c()
qp05_n_ec <- c()
qp25_n_ec <- c()
qp50_n_ec <- c()
qp75_n_ec <- c()
qp95_n_ec <- c()

#--- Run simulations for coarseness >= 2 days
for(i in 1:nvec_len){
  
  n <- nvec[i]
  print(n)
  
  xsim <- list()
  
  for(s in 1:S){
    simdata <- simSI(muS = meanSI, sdS = sdSI, n = n)
    xsim[[s]] <- data.frame(sl = simdata$sl, sr = simdata$sr, sw = simdata$sw)
  }
  
  # Hosting nonparametric estimates
  mu_S   <- data.frame(matrix(0, nrow = S, ncol = 5))
  sd_S   <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp05_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp25_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp50_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp75_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp95_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  colnames(mu_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(sd_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp05_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp25_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp50_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp75_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp95_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  
  for(s in 1:S){
    x <- xsim[[s]][,c(1,2)]
    fit <- estimSI(x = x, nboot = nboot)
    
    mu_S[s, ]   <- fit$npestim$mean[-2]
    sd_S[s, ]   <- fit$npestim$sd[-2]
    qp05_S[s, ] <- fit$npestim$q0.05[-2]
    qp25_S[s, ] <- fit$npestim$q0.25[-2]
    qp50_S[s, ] <- fit$npestim$q0.5[-2]
    qp75_S[s, ] <- fit$npestim$q0.75[-2]
    qp95_S[s, ] <- fit$npestim$q0.95[-2]
    
  }
  
  mu_n[i]   <- mean(mu_S$point)
  sd_n[i]   <- mean(sd_S$point)
  qp05_n[i] <- mean(qp05_S$point)
  qp25_n[i] <- mean(qp25_S$point)
  qp50_n[i] <- mean(qp50_S$point)
  qp75_n[i] <- mean(qp75_S$point)
  qp95_n[i] <- mean(qp95_S$point)

}

#--- Run simulations for epsilon-coarseness 
for(i in 1:nvec_len){
  
  n <- nvec[i]
  print(n)
  
  xsim <- list()
  
  for(s in 1:S){
    si_true <- rnorm(n = n, mean = meanSI, sd = sdSI)
    si_left  <- si_true - (eps / 2)
    si_right <- si_true + (eps / 2)
    sw <- si_right - si_left
    xsim[[s]] <- data.frame(sl = si_left, sr = si_right, sw = sw)
  }
  
  # Hosting nonparametric estimates
  mu_S   <- data.frame(matrix(0, nrow = S, ncol = 5))
  sd_S   <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp05_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp25_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp50_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp75_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  qp95_S <- data.frame(matrix(0, nrow = S, ncol = 5))
  colnames(mu_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(sd_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp05_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp25_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp50_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp75_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  colnames(qp95_S) <- c("point", "ci90l","ci90r", "ci95l", "ci95r")
  
  for(s in 1:S){
    x <- xsim[[s]][,c(1,2)]
    fit <- estimSI(x = x, nboot = nboot)
    
    mu_S[s, ]   <- fit$npestim$mean[-2]
    sd_S[s, ]   <- fit$npestim$sd[-2]
    qp05_S[s, ] <- fit$npestim$q0.05[-2]
    qp25_S[s, ] <- fit$npestim$q0.25[-2]
    qp50_S[s, ] <- fit$npestim$q0.5[-2]
    qp75_S[s, ] <- fit$npestim$q0.75[-2]
    qp95_S[s, ] <- fit$npestim$q0.95[-2]

  }
  
  mu_n_ec[i]   <- mean(mu_S$point)
  sd_n_ec[i]   <- mean(sd_S$point)
  qp05_n_ec[i] <- mean(qp05_S$point)
  qp25_n_ec[i] <- mean(qp25_S$point)
  qp50_n_ec[i] <- mean(qp50_S$point)
  qp75_n_ec[i] <- mean(qp75_S$point)
  qp95_n_ec[i] <- mean(qp95_S$point)
  
}


pdf(file = "Figure-3.pdf", width = 15, height = 9)

# Plot asymptotic behavior of mean estimates under coarseness >=2 days

par(mfrow = c(2,3))
ftsize <- 1.8

# Estimated mean SI
plot(nvec, mu_n, type = "l", col = "black", xlab = "Sample size",
     ylab = "", main = "Mean of SI", 
     cex.axis = ftsize, cex.lab = ftsize, cex.main = ftsize, yaxt = "n")
axis(2, at = meanSI, label = expression(mu[S]), cex.axis = ftsize)
abline(h = meanSI, lty = 3)
lines(nvec, mu_n, type = "l", lwd = 2)
legend("topright", lty = 1, lwd = 2, c("Mean estimate"), cex = ftsize, 
       bty = "n")

# Estimated SI standard deviation
plot(nvec, sd_n, type = "l", col = "black", xlab = "Sample size",
     ylab = "", main = "Standard deviation of SI",
     cex.axis = ftsize, cex.lab = ftsize, cex.main = ftsize, yaxt = "n")
axis(2, at = sdSI, label = expression(sigma[S]), cex.axis = ftsize)
abline(h = sdSI, lty = 3)
lines(nvec, sd_n, type = "l", lwd = 2)
legend("bottomright", lty = 1, lwd = 2, c("Mean estimate"), cex = ftsize, 
       bty = "n")

# Estimated quantiles
plot(nvec, qp05_n, type = "l", col = "black", xlab = "Sample size",
     ylab = "", main = "Selected SI quantiles",
     ylim = c(-3,10), cex.axis = ftsize, cex.lab = ftsize, 
     yaxt = "n", cex.main = ftsize)
abline(h = qtrue, lty = 3)
axis(2, at = qtrue, c(expression(q[0.05]), expression(q[0.25]),
                      expression(q[0.50]), expression(q[0.75]),
                      expression(q[0.95])), cex.axis = ftsize)
lines(nvec, qp05_n, type = "l", lwd = 2)
lines(nvec, qp25_n, type = "l", lwd = 2)
lines(nvec, qp50_n, type = "l", lwd = 2)
lines(nvec, qp75_n, type = "l", lwd = 2)
lines(nvec, qp95_n, type = "l", lwd = 2)
legend("topright", lty = 1, lwd = 2, c("Mean estimate"), cex = ftsize, 
       bty = "n")

# Plot asymptotic behavior of mean estimates under epsilon-coarseness

# Estimated mean SI
plot(nvec, mu_n_ec, type = "l", col = "black", xlab = "Sample size",
     ylab = "", main = "Mean of SI", 
     cex.axis = ftsize, cex.lab = ftsize, cex.main = ftsize, yaxt = "n")
axis(2, at = meanSI, label = expression(mu[S]), cex.axis = ftsize)
abline(h = meanSI, lty = 3)
lines(nvec, mu_n_ec, type = "l", lwd = 2)
legend("topright", lty = 1, lwd = 2, c("Mean estimate"), cex = ftsize, 
       bty = "n")

# Estimated SI standard deviation
plot(nvec, sd_n_ec, type = "l", col = "black", xlab = "Sample size",
     ylab = "", main = "Standard deviation of SI",
     cex.axis = ftsize, cex.lab = ftsize, cex.main = ftsize, yaxt = "n")
axis(2, at = sdSI, label = expression(sigma[S]), cex.axis = ftsize)
abline(h = sdSI, lty = 3)
lines(nvec, sd_n_ec, type = "l", lwd = 2)
legend("bottomright", lty = 1, lwd = 2, c("Mean estimate"), cex = ftsize, 
       bty = "n")

# Estimated quantiles
plot(nvec, qp05_n_ec, type = "l", col = "black", xlab = "Sample size",
     ylab = "", main = "Selected SI quantiles",
     ylim = c(-3,10), cex.axis = ftsize, cex.lab = ftsize, 
     yaxt = "n", cex.main = ftsize)
abline(h = qtrue, lty = 3)
axis(2, at = qtrue, c(expression(q[0.05]), expression(q[0.25]),
                      expression(q[0.50]), expression(q[0.75]),
                      expression(q[0.95])), cex.axis = ftsize)
lines(nvec, qp05_n_ec, type = "l", lwd = 2)
lines(nvec, qp25_n_ec, type = "l", lwd = 2)
lines(nvec, qp50_n_ec, type = "l", lwd = 2)
lines(nvec, qp75_n_ec, type = "l", lwd = 2)
lines(nvec, qp95_n_ec, type = "l", lwd = 2)
legend("topright", lty = 1, lwd = 2, c("Mean estimate"), cex = ftsize, 
       bty = "n")

dev.off()















