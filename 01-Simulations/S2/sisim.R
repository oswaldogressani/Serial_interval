#------------------------------------------------------------------------------#
#                      Simulation of illness onset data                        #
#                    Oswaldo Gressani. All rights reserved.                    #
#------------------------------------------------------------------------------#


sisim <- function(n = 10, dist = c("gaussian","weibull"),
                  elast = c(a = 4, b = 4), c = 10, excen = 5){
  
  # Draw a sample of size n from the population SI distribution
  if(match.arg(dist) == "gaussian") {
    mean <- 2.8
    sd <- 2.5
    s <- rnorm(n = n, mean = mean, sd = sd)
    perc <- qnorm(p=c(0.05,0.25,0.5,0.75,0.95), mean = mean, sd = sd)
    Sfeat <- list(mean = mean, sd = sd, q05 = perc[1], q025 = perc[2],
                  q50 = perc[3], q075 = perc[4], q095 = perc[5])
  } else if (match.arg(dist) == "weibull") {
    shape <- 2.36
    scale <- 3.18
    s <- rweibull(n = n, shape = shape, scale = scale)
    perc <- qweibull(p=c(0.05,0.25,0.5,0.75,0.95), shape = shape, scale = scale)
    weibmean <- scale * gamma(1+1/shape)
    weibsd <- sqrt((scale^2)* (gamma(1+2/shape)-(gamma(1+1/shape))^2))
    Sfeat <- list(mean = weibmean, sd = weibsd, q05 = perc[1], 
                  q025 = perc[2], q50 = perc[3], q075 = perc[4], q095 = perc[5])
  }
  
  # Compute illness onset time for infector (1) and infectee (2)
  t1 <- runif(n = n, min = abs(s) + c, abs(s) + c + 1) # Infector
  t2 <- t1 + s                                         # Infectee
  
  # Generation of interval censored data
  
  # Left and right boundaries of illness onset data for infector
  Delta1 <- rgamma(n = n, shape = elast[1], rate = elast[2])
  rho1 <- rbeta(n = n, shape1 = excen, shape2 = excen)
  t1L <- t1 - rho1 * Delta1
  t1R <- t1 + (1 - rho1) * Delta1
  
  # Left and right boundaries of illness onset data for infectee
  Delta2 <- rgamma(n = n, shape = elast[1], rate = elast[2])
  rho2 <- rbeta(n = n, shape1 = excen, shape2 = excen)
  t2L <- t2 - rho2 * Delta2
  t2R <- t2 + (1 - rho2) * Delta2
  
  # Single interval censoring 
  sLs <- t2 - t1R
  sRs <- t2 - t1L
  
  # Double interval censoring
  sLd <- t2L-t1R
  sRd <- t2R-t1L
  
  # Check constraints on illness onset time
  ck <- all(t1 >= 0 & t2 >= 0) & all(t1L >= 0 & t2L >= 0) & 
    all(t1R > t1L & t2R > t2L)
  if(isTRUE(ck)){
    constr <- "All constraints satisfied"
    constrval <- 1
  } else{
    constr <- "Constraints not satisfied"
    constrval <- 0
  }
  
  data <- list(exact = s, single = data.frame(sL = sLs, sR = sRs),
       double = data.frame(sL = sLd, sR = sRd),
       Sfeat = Sfeat,
       check = constr, checkval = constrval)
  
  return(data)
}
  
  
  








