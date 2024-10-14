#------------------------------------------------------------------------------#
#       Assessing performance of nonparametric SI estimation with EpiLPS       #
#                  Function for simulation of various scenarios.               #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#


sisimul <- function(n = 10, dist = c("gaussian","weibull"),
                    censor = c("double","single","none"), 
                    seedval = 123){
  
  set.seed(seedval)
  
  #----------------------- Simulations
  
  J <- 500   # Number of simulations
  
  # Record mean
  si_mean_vec <- c()
  si_mean_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_mean_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_mean_width <- c()
  CI95_mean_width <- c()
  
  # Record sd
  si_sd_vec <- c()
  si_sd_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_sd_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_sd_width <- c()
  CI95_sd_width <- c()
  
  # Record q05
  si_q05_vec <- c()
  si_q05_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_q05_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_q05_width <- c()
  CI95_q05_width <- c()
  
  # Record q25
  si_q25_vec <- c()
  si_q25_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_q25_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_q25_width <- c()
  CI95_q25_width <- c()
  
  # Record q50
  si_q50_vec <- c()
  si_q50_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_q50_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_q50_width <- c()
  CI95_q50_width <- c()
  
  # Record q75
  si_q75_vec <- c()
  si_q75_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_q75_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_q75_width <- c()
  CI95_q75_width <- c()
  
  # Record q95
  si_q95_vec <- c()
  si_q95_CI90p_mat <- matrix(0, nrow = J, ncol = 2)
  si_q95_CI95p_mat <- matrix(0, nrow = J, ncol = 2)
  CI90_q95_width <- c()
  CI95_q95_width <- c()
  
  # Check constraints of generated data
  ck <- c()
  
  progbar <- utils::txtProgressBar(min = 1, max = J, initial = 1,
                                   style = 3, char =">")
  
  for(j in 1:J){
    
    SIdist <- match.arg(dist)
    
    # Generate data
    sidata <- sisim(n = n, dist = SIdist)
    
    censoringtype <- match.arg(censor)
    
    if(censoringtype == "double"){
      x <- sidata$double
    } else if (censoringtype == "single"){
      x <- sidata$single
    } else if (censoringtype == "none"){
      xtemp <- sidata$exact
      x <- data.frame(sL = xtemp - 0.5, sR = xtemp + 0.5)
    }
    
    ck[j] <- sidata$checkval
    
    fitSI <- EpiLPS::estimSI(x = x, verbose = FALSE)
    fitSI <- fitSI$estim
    
    si_mean_vec[j] <- fitSI$Estim[1]
    si_mean_CI90p_mat[j,] <- c(fitSI$CI90p_l[1], fitSI$CI90p_r[1])
    si_mean_CI95p_mat[j,] <- c(fitSI$CI95p_l[1], fitSI$CI95p_r[1])
    CI90_mean_width[j] <- diff(si_mean_CI90p_mat[j,])
    CI95_mean_width[j] <- diff(si_mean_CI95p_mat[j,])
    
    si_sd_vec[j] <- fitSI$Estim[2]
    si_sd_CI90p_mat[j,] <- c(fitSI$CI90p_l[2], fitSI$CI90p_r[2])
    si_sd_CI95p_mat[j,] <- c(fitSI$CI95p_l[2], fitSI$CI95p_r[2])
    CI90_sd_width[j] <- diff(si_sd_CI90p_mat[j,])
    CI95_sd_width[j] <- diff(si_sd_CI95p_mat[j,])
    
    si_q05_vec[j] <- fitSI$Estim[3]
    si_q05_CI90p_mat[j,] <- c(fitSI$CI90p_l[3], fitSI$CI90p_r[3])
    si_q05_CI95p_mat[j,] <- c(fitSI$CI95p_l[3], fitSI$CI95p_r[3])
    CI90_q05_width[j] <- diff(si_q05_CI90p_mat[j,])
    CI95_q05_width[j] <- diff(si_q05_CI95p_mat[j,])
    
    si_q25_vec[j] <- fitSI$Estim[4]
    si_q25_CI90p_mat[j,] <- c(fitSI$CI90p_l[4], fitSI$CI90p_r[4])
    si_q25_CI95p_mat[j,] <- c(fitSI$CI95p_l[4], fitSI$CI95p_r[4])
    CI90_q25_width[j] <- diff(si_q25_CI90p_mat[j,])
    CI95_q25_width[j] <- diff(si_q25_CI95p_mat[j,])
    
    si_q50_vec[j] <- fitSI$Estim[5]
    si_q50_CI90p_mat[j,] <- c(fitSI$CI90p_l[5], fitSI$CI90p_r[5])
    si_q50_CI95p_mat[j,] <- c(fitSI$CI95p_l[5], fitSI$CI95p_r[5])
    CI90_q50_width[j] <- diff(si_q50_CI90p_mat[j,])
    CI95_q50_width[j] <- diff(si_q50_CI95p_mat[j,])
    
    si_q75_vec[j] <- fitSI$Estim[6]
    si_q75_CI90p_mat[j,] <- c(fitSI$CI90p_l[6], fitSI$CI90p_r[6])
    si_q75_CI95p_mat[j,] <- c(fitSI$CI95p_l[6], fitSI$CI95p_r[6])
    CI90_q75_width[j] <- diff(si_q75_CI90p_mat[j,])
    CI95_q75_width[j] <- diff(si_q75_CI95p_mat[j,])
    
    si_q95_vec[j] <- fitSI$Estim[7]
    si_q95_CI90p_mat[j,] <- c(fitSI$CI90p_l[7], fitSI$CI90p_r[7])
    si_q95_CI95p_mat[j,] <- c(fitSI$CI95p_l[7], fitSI$CI95p_r[7])
    CI90_q95_width[j] <- diff(si_q95_CI90p_mat[j,])
    CI95_q95_width[j] <- diff(si_q95_CI95p_mat[j,])
    
    utils::setTxtProgressBar(progbar, j)
  }
  close(progbar)
  
  # Compute performance metrics
  
  #--- Performance mean
  Bias_mean_si <- mean(si_mean_vec - sidata$Sfeat$mean)
  ESE_mean_si <- sqrt((1/(J-1)) * sum((si_mean_vec-mean(si_mean_vec))^2))
  RMSE_mean_si <- sqrt(mean(((si_mean_vec-sidata$Sfeat$mean)^2)))
  CP90p_mean_si <- mean((sidata$Sfeat$mean >= si_mean_CI90p_mat[,1]) & 
                          (sidata$Sfeat$mean <= si_mean_CI90p_mat[,2])) * 100
  CP95p_mean_si <- mean((sidata$Sfeat$mean >= si_mean_CI95p_mat[,1]) & 
                          (sidata$Sfeat$mean <= si_mean_CI95p_mat[,2])) * 100
  Median_mean_CI90 <- median(CI90_mean_width)
  Median_mean_CI95 <- median(CI95_mean_width)
  
  #--- Performance sd
  Bias_sd_si <- mean(si_sd_vec - sidata$Sfeat$sd)
  ESE_sd_si <- sqrt((1/(J-1)) * sum((si_sd_vec-mean(si_sd_vec))^2))
  RMSE_sd_si <- sqrt(mean(((si_sd_vec-sidata$Sfeat$sd)^2)))
  CP90p_sd_si <- mean((sidata$Sfeat$sd >= si_sd_CI90p_mat[,1]) & 
                        (sidata$Sfeat$sd <= si_sd_CI90p_mat[,2])) * 100
  CP95p_sd_si <- mean((sidata$Sfeat$sd >= si_sd_CI95p_mat[,1]) & 
                        (sidata$Sfeat$sd <= si_sd_CI95p_mat[,2])) * 100
  Median_sd_CI90 <- median(CI90_sd_width)
  Median_sd_CI95 <- median(CI95_sd_width)
  
  #--- Performance q05
  tarq05 <- sidata$Sfeat$q05
  Bias_q05_si <- mean(si_q05_vec - tarq05)
  ESE_q05_si <- sqrt((1/(J-1)) * sum((si_q05_vec-mean(si_q05_vec))^2))
  RMSE_q05_si <- sqrt(mean(((si_q05_vec-tarq05)^2)))
  CP90p_q05_si <- mean((tarq05 >= si_q05_CI90p_mat[,1]) & 
                         (tarq05 <= si_q05_CI90p_mat[,2])) * 100
  CP95p_q05_si <- mean((tarq05 >= si_q05_CI95p_mat[,1]) & 
                         (tarq05 <= si_q05_CI95p_mat[,2])) * 100
  Median_q05_CI90 <- median(CI90_q05_width)
  Median_q05_CI95 <- median(CI95_q05_width)
  
  #--- Performance q25
  tarq25 <- sidata$Sfeat$q025
  Bias_q25_si <- mean(si_q25_vec - tarq25)
  ESE_q25_si <- sqrt((1/(J-1)) * sum((si_q25_vec-mean(si_q25_vec))^2))
  RMSE_q25_si <- sqrt(mean(((si_q25_vec-tarq25)^2)))
  CP90p_q25_si <- mean((tarq25 >= si_q25_CI90p_mat[,1]) & 
                         (tarq25 <= si_q25_CI90p_mat[,2])) * 100
  CP95p_q25_si <- mean((tarq25 >= si_q25_CI95p_mat[,1]) & 
                         (tarq25 <= si_q25_CI95p_mat[,2])) * 100
  Median_q25_CI90 <- median(CI90_q25_width)
  Median_q25_CI95 <- median(CI95_q25_width)
  
  #--- Performance q50
  tarq50 <- sidata$Sfeat$q50
  Bias_q50_si <- mean(si_q50_vec - tarq50)
  ESE_q50_si <- sqrt((1/(J-1)) * sum((si_q50_vec-mean(si_q50_vec))^2))
  RMSE_q50_si <- sqrt(mean(((si_q50_vec-tarq50)^2)))
  CP90p_q50_si <- mean((tarq50 >= si_q50_CI90p_mat[,1]) & 
                         (tarq50 <= si_q50_CI90p_mat[,2])) * 100
  CP95p_q50_si <- mean((tarq50 >= si_q50_CI95p_mat[,1]) & 
                         (tarq50 <= si_q50_CI95p_mat[,2])) * 100
  Median_q50_CI90 <- median(CI90_q50_width)
  Median_q50_CI95 <- median(CI95_q50_width)
  
  #--- Performance q75
  tarq75 <- sidata$Sfeat$q075
  Bias_q75_si <- mean(si_q75_vec - tarq75)
  ESE_q75_si <- sqrt((1/(J-1)) * sum((si_q75_vec-mean(si_q75_vec))^2))
  RMSE_q75_si <- sqrt(mean(((si_q75_vec-tarq75)^2)))
  CP90p_q75_si <- mean((tarq75 >= si_q75_CI90p_mat[,1]) & 
                         (tarq75 <= si_q75_CI90p_mat[,2])) * 100
  CP95p_q75_si <- mean((tarq75 >= si_q75_CI95p_mat[,1]) & 
                         (tarq75 <= si_q75_CI95p_mat[,2])) * 100
  Median_q75_CI90 <- median(CI90_q75_width)
  Median_q75_CI95 <- median(CI95_q75_width)
  
  #--- Performance q95
  tarq95 <- sidata$Sfeat$q095
  Bias_q95_si <- mean(si_q95_vec - tarq95)
  ESE_q95_si <- sqrt((1/(J-1)) * sum((si_q95_vec-mean(si_q95_vec))^2))
  RMSE_q95_si <- sqrt(mean(((si_q95_vec-tarq95)^2)))
  CP90p_q95_si <- mean((tarq95 >= si_q95_CI90p_mat[,1]) & 
                         (tarq95 <= si_q95_CI90p_mat[,2])) * 100
  CP95p_q95_si <- mean((tarq95 >= si_q95_CI95p_mat[,1]) & 
                         (tarq95 <= si_q95_CI95p_mat[,2])) * 100
  Median_q95_CI90 <- median(CI90_q95_width)
  Median_q95_CI95 <- median(CI95_q95_width)
  
  
  Simperf <- matrix(0, nrow = 7, ncol = 7)
  rownames(Simperf) <- c("Mean", "SD","q05","q25","q50","q75","q95")
  colnames(Simperf) <- c("Bias", "ESE","RMSE","90%CP_p","95%CP_p",
                         "90%CI width","95%CI width")
  Simperf[,1] <- c(Bias_mean_si, Bias_sd_si, Bias_q05_si,Bias_q25_si, 
                   Bias_q50_si, Bias_q75_si, Bias_q95_si)
  Simperf[,2] <- c(ESE_mean_si, ESE_sd_si, ESE_q05_si, ESE_q25_si, 
                   ESE_q50_si, ESE_q75_si, ESE_q95_si)
  Simperf[,3] <- c(RMSE_mean_si, RMSE_sd_si, RMSE_q05_si, RMSE_q25_si, 
                   RMSE_q50_si, RMSE_q75_si, RMSE_q95_si)
  Simperf[,4] <- c(CP90p_mean_si, CP90p_sd_si, CP90p_q05_si,
                   CP90p_q25_si, CP90p_q50_si, CP90p_q75_si, CP90p_q95_si)
  Simperf[,5] <- c(CP95p_mean_si, CP95p_sd_si, CP95p_q05_si,
                   CP95p_q25_si, CP95p_q50_si, CP95p_q75_si, CP95p_q95_si)
  Simperf[,6] <- c(Median_mean_CI90, Median_sd_CI90, Median_q05_CI90,
                   Median_q25_CI90, Median_q50_CI90, Median_q75_CI90,
                   Median_q95_CI90)
  Simperf[,7] <- c(Median_mean_CI95, Median_sd_CI95, Median_q05_CI95,
                   Median_q25_CI95, Median_q50_CI95, Median_q75_CI95,
                   Median_q95_CI95)
  
  outlist <- list(Simperf = Simperf, ck = ck, n=n, 
                  dist = SIdist, censor = censoringtype, seedval = seedval)
  
}