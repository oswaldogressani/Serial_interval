#------------------------------------------------------------------------------#
#                  Estimation of serial intervals with EpiLPS                  #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#

estimSI_boot <- function(x, B = 5000){

sL <- x[, 1]
sR <- x[, 2]
ds <- sR - sL
n <- nrow(x)
ninv <- 1 / n

Fhat <- function(s) ninv *  sum((s - sL) / ds * (s >= sL & s < sR) +  (s >= sR))

# The bootstrap
sorder <- sort(c(sL,sR))
Fsorder <- sapply(sorder, Fhat)
diffF <- diff(Fsorder)

Sfeat <- data.frame(meanboot = double(), sdboot = double(), q005boot = double(),
                    q025boot = double(), q050boot = double(), 
                    q075boot = double(), q095boot = double())

# Record bootstrap samples
bootmat <- matrix(0, nrow = B, ncol = n)

# Write loop below in C++
for (b in 1:B) {
  bootsample <- c()
  for(i in 1:n){
    u <- runif(1)
    idxr <-  min(which(u<=Fsorder))
    idxl <- idxr-1
    deltaF <- diffF[idxl]
    if(deltaF>0){
      bootsample[i] <- ((Fsorder[idxr]-u)*sorder[idxl]+(u-Fsorder[idxl])*sorder[idxr])/deltaF
    }else{
      bootsample[i] <- sorder[idxl]
    }
  }
  
  bootmat[b,] <- bootsample
  Sfeat[b,] <- c(mean(bootsample), sd(bootsample),
                 quantile(bootsample, probs = c(0.05,0.25,0.5,0.75,0.95)))
}

#----- Point and interval estimate of mean of SI
si_mean <- mean(Sfeat$meanboot)
si_mean_CI90_p <- quantile(Sfeat$meanboot, probs = c(0.05,0.95))
si_mean_CI95_p <- quantile(Sfeat$meanboot, probs = c(0.025,0.975))

#----- Point and interval estimate of sd of SI
si_sd <- mean(Sfeat$sdboot)
si_sd_CI90_p <- quantile(Sfeat$sdboot, probs = c(0.05,0.95))
si_sd_CI95_p <- quantile(Sfeat$sdboot, probs = c(0.025,0.975))

#----- Point and interval estimate of 5th percentile of SI 
si_q005 <- mean(Sfeat$q005boot)
si_q005_CI90_p <- quantile(Sfeat$q005boot, probs = c(0.05, 0.95))
si_q005_CI95_p <- quantile(Sfeat$q005boot, probs = c(0.025, 0.975))

#----- Point and interval estimate of 25th percentile of SI 
si_q025 <- mean(Sfeat$q025boot)
si_q025_CI90_p <- quantile(Sfeat$q025boot, probs = c(0.05, 0.95))
si_q025_CI95_p <- quantile(Sfeat$q025boot, probs = c(0.025, 0.975))

#----- Point and interval estimate of 50th percentile of SI 
si_q050 <- mean(Sfeat$q050boot)
si_q050_CI90_p <- quantile(Sfeat$q050boot, probs = c(0.05, 0.95))
si_q050_CI95_p <- quantile(Sfeat$q050boot, probs = c(0.025, 0.975))

#----- Point and interval estimate of 75th percentile of SI 
si_q075 <- mean(Sfeat$q075boot)
si_q075_CI90_p <- quantile(Sfeat$q075boot, probs = c(0.05, 0.95))
si_q075_CI95_p <- quantile(Sfeat$q075boot, probs = c(0.025, 0.975))

#----- Point and interval estimate of 95th percentile of SI 
si_q095 <- mean(Sfeat$q095boot)
si_q095_CI90_p <- quantile(Sfeat$q095boot, probs = c(0.05, 0.95))
si_q095_CI95_p <- quantile(Sfeat$q095boot, probs = c(0.025, 0.975))


outlist <- list(estim = data.frame(Estim = double(), 
                      CI90p_l = double(), CI90p_r = double(),
                      CI95p_l = double(), CI95p_r = double()),
                bootsamples = bootmat,
                Fhat = Fhat)
outlist$estim[1,] <- c(si_mean, si_mean_CI90_p, si_mean_CI95_p)
outlist$estim[2,] <- c(si_sd, si_sd_CI90_p, si_sd_CI95_p)
outlist$estim[3,] <- c(si_q005, si_q005_CI90_p, si_q005_CI95_p)
outlist$estim[4,] <- c(si_q025, si_q025_CI90_p, si_q025_CI95_p)
outlist$estim[5,] <- c(si_q050, si_q050_CI90_p, si_q050_CI95_p)
outlist$estim[6,] <- c(si_q075, si_q075_CI90_p, si_q075_CI95_p)
outlist$estim[7,] <- c(si_q095, si_q095_CI90_p, si_q095_CI95_p)
rownames(outlist$estim) <- c("mean","sd","q05", "q25","q50","q75", "q95")

return(outlist)
}
