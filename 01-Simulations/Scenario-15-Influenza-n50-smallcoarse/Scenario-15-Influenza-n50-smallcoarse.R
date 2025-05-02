#------------------------------------------------------------------------------#
#                         Scenario 15 (Influenza A)                            #
#------------------------------------------------------------------------------#
#   -->  Gaussian SI with mean=2.1 and sd=1.2 from Vink et al. (2014)          #
#   -->  Sample size n = 50                                                    #
#   -->  Virtually no coarseness (epsilon-coarseness)                          #
#------------------------------------------------------------------------------#

rm(list = ls())

library("EpiDelays")
library("xtable")

#---- Tune scenario setting here
S <- 1000
set.seed(2017)
n <- 50
meanSI  <- 2.1  # Vink et al. (2014)
sdSI    <- 1.2  # Vink et al. (2014)
qtrue   <- qnorm(p = c(0.05,0.25,0.50,0.75,0.95), mean = meanSI, sd = sdSI)
nboot <- 2000
#-------------------------------

xsim <- list()

for(s in 1:S){
  si_true <- rnorm(n = n, mean = meanSI, sd = sdSI)
  si_left  <- si_true - 0.005
  si_right <- si_true + 0.005
  sw <- si_right - si_left
  xsim[[s]] <- data.frame(sl = si_left, sr = si_right, sw = sw)
}

# Summary statistics of observed serial interval windows
swobs <- data.frame(sw = as.numeric(unlist(lapply(xsim, "[", 3))))
freqswobs <- table(swobs)
meancoarseness <- round(mean(swobs$sw))
barplot(freqswobs, main = paste0("Mean coarseness: ", meancoarseness, " days"),
        ylim = c(0, max(freqswobs)),
        xlab = "SI window width", ylab = "Frequency")


# Estimation

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


progbar <- utils::txtProgressBar(min = 1, max = S, initial = 1,
                                 style = 3, char =">")

for(s in 1:S){
  x <- xsim[[s]][,c(1,2)]
  fit <- estimSI(x = x, nboot = nboot)
  
  # Record nonparametric estimates
  mu_S[s, ]   <- fit$npestim$mean[-2]
  sd_S[s, ]   <- fit$npestim$sd[-2]
  qp05_S[s, ] <- fit$npestim$q0.05[-2]
  qp25_S[s, ] <- fit$npestim$q0.25[-2]
  qp50_S[s, ] <- fit$npestim$q0.5[-2]
  qp75_S[s, ] <- fit$npestim$q0.75[-2]
  qp95_S[s, ] <- fit$npestim$q0.95[-2]
  
  utils::setTxtProgressBar(progbar, s)
}

close(progbar)


# Performance
perf <- function(xS, true){
  point <- xS$point
  ci90l <- xS$ci90l
  ci90r <- xS$ci90r
  ci95l <- xS$ci95l
  ci95r <- xS$ci95r
  meanpoint <- mean(point)
  bias  <- mean(point - true)
  ese   <- sqrt((1/(S-1)) * sum((point - meanpoint)^2))
  rmse  <- sqrt(mean((point - true)^2))
  cp90  <- mean((ci90l <= true) & (true <= ci90r)) * 100
  cp95  <- mean((ci95l <= true) & (true <= ci95r)) * 100
  ci90w <- median(ci90r - ci90l)
  ci95w <- median(ci95r - ci95l)
  perfout <- data.frame(bias = bias, ese = ese, rmse = rmse,
                        cp90 = cp90, cp95 = cp95, ci90w = ci90w,
                        ci95w = ci95w)
  return(perfout)
}


# Performance of nonparametric method
perfnonparam <- matrix(0, nrow = 7, ncol = 7)
rownames(perfnonparam) <- c("Mean", "SD","q05","q25","q50","q75","q95")
colnames(perfnonparam) <- c("Bias", "ESE","RMSE","90%CP","95%CP",
                            "90%CIw","95%CIw")
perfnonparam[1, ] <- as.numeric(perf(mu_S, true = meanSI))
perfnonparam[2, ] <- as.numeric(perf(sd_S, true = sdSI))
perfnonparam[3, ] <- as.numeric(perf(qp05_S, true = qtrue[1]))
perfnonparam[4, ] <- as.numeric(perf(qp25_S, true = qtrue[2]))
perfnonparam[5, ] <- as.numeric(perf(qp50_S, true = qtrue[3]))
perfnonparam[6, ] <- as.numeric(perf(qp75_S, true = qtrue[4]))
perfnonparam[7, ] <- as.numeric(perf(qp95_S, true = qtrue[5]))

# Output results in txt format
round(perfnonparam, 3)

print(xtable(perfnonparam, digits = 3,
             label = "Scenario-15-Influenza-n50-smallcoarse"), 
      file="Scenario-15-Influenza-n50-smallcoarse.txt")


















