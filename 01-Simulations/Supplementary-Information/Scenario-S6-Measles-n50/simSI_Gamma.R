#------------------------------------------------------------------------------#
#               Simulation of artificial serial interval data                  #
#                   Oswaldo Gressani. All rights reserved.                     #
#------------------------------------------------------------------------------#


simSI_Gamma <- function(muS = 9.9, sdS = 2.4, maxcoarse = 3, probcoarse = c(0.8,0.15,0.05), n = 10){
 as <- (muS / sdS)^2
 bs <- muS / (sdS^2)
 s <- rgamma(n = n, shape = as, rate = bs)
 coarseness <- seq_len(maxcoarse)
 c <- sample(coarseness, size = n, replace = TRUE, prob = probcoarse)
 u <- runif(n, min = 0, max = 1)
 sl <- floor((s - u * c))
 sr <- ceiling((s + (1 - u) * c))
 sw <- sr - sl
 x <- data.frame(s = s, sl = sl, sr = sr, sw = sw)
 return(x)
}