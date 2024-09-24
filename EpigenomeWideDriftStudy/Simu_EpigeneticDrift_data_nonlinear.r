# simu
library("lmtest")
library("sandwich")
library(dplyr)
library(tidyr)
# 1. Age
Age <- rnorm(10000, 50, 10)

# 2. CpG
SimCpGdata <- function(type, nSamples, Age, refcpg, b0, b1, b2, minX, maxX, meanResid, sdResid, meanError, sdError, nIterations){
                    data <- matrix(NA, nSamples, nIterations+1)
                    set.seed(111)
                    # sample Age x
                    x <- sample(Age[which(Age> minX & Age < maxX)], nSamples)
                    data[,1] <- x
                    if (type == "drift") {  #drif y
                                 for(i in (2:(nIterations/2+1))){
                                    resid <- rnorm(nSamples,meanResid, sd=sdResid+b2*x)
                                    y <- b0 + resid + rnorm(nSamples, meanError, sdError)
                                    beta_y <- (y-min(y))/(max(y)-min(y))
  							        data[,i] <- beta_y
                                    }}
                                 for(i in (nIterations/2+2):(nIterations+1)){
                                    y <- rep(NA, nSamples)
									# negative
									index1 <- which(x<=45)
									x1 <- x[index1]
									resid <- rnorm(length(x1),meanResid, sd=100-b2*10*x1)
                                    y1 <- b0 + resid + rnorm(length(x1), meanError, sdError)
									beta_y1 <- (y1-min(y1))/(max(y1)-min(y1))
									y[index1] <- beta_y1
									
									# null
									index2 <- which(x>45 & x<55)
									x2 <- x[index2]
									resid <- rnorm(length(x2),meanResid, sd=sdResid)
                                    y2 <- b0 + resid + rnorm(length(x2), meanError, sdError)
									beta_y2 <- (y2-min(y2))/(max(y2)-min(y2))
									y[index2] <- beta_y2				
									
									# positive
									index3 <- which(x>=55) 
									x3 <- x[index3]
									resid <- rnorm(length(x3),meanResid, sd=sdResid+b2*10*x3)
                                    y3 <- b0 + resid + rnorm(length(x3), meanError, sdError)
									beta_y3 <- (y3-min(y3))/(max(y3)-min(y3))
									y[index3] <- beta_y3
									
  							        data[,i] <- y
                                    }
 return(data)
}						
