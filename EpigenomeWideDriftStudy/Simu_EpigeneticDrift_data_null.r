# Simulate
library("lmtest")
library("sandwich")
library(dplyr)
library(tidyr)

SimCpGdata <- function(type, nSamples, Age, refcpg, b0, b1, b2, minX, maxX, meanResid, sdResid, nIterations){
                    data <- matrix(NA, nSamples, nIterations+1)
                    set.seed(111)
                    # sample age x
                    x <- sample(Age[which(Age> minX & Age < maxX)], nSamples)
                    data[,1] <- x                      
                               if(type == "null"){# null y
                                  for(i in 2:(nIterations+1)){
                                    y <- sample(refcpg[,sample(1:ncol(refcpg),1)], nSamples, replace=T)
                                    data[,i] <- y
                                    }}

                               if (type == "clock") { # clock y
                                 for(i in 2:(nIterations+1)){
                                    y <- b0 + b1*x + rnorm(nSamples,meanResid ,sdResid)
                                    beta_y <- (y-min(y))/(max(y)-min(y))
                                    data[,i] <- beta_y
                                    }}

							   if (type == "drift") {  # drift y
                                 for(i in 2:(nIterations+1)){
                                    resid <- rnorm(nSamples,meanResid, sd=sdResid+b2*x)
                                    y <- b0 + resid
                                    beta_y <- (y-min(y))/(max(y)-min(y))
  							        data[,i] <- beta_y
                                    }}
								if (type == "drift2age") {  # drift2 y
                                 for(i in 2:(nIterations+1)){
                                    resid <- rnorm(nSamples,meanResid, sd=sdResid+b2*x*x)
                                    y <- b0 + resid
                                    beta_y <- (y-min(y))/(max(y)-min(y))
  							        data[,i] <- beta_y
                                    }}

                               if (type == "clokdrift") { # clock and drift y
                                 for(i in 2:(nIterations+1)){
                                    resid <- rnorm(nSamples,meanResid, sd=sdResid+b2*x)
                                    y <- b0 + b1*x + resid
                                    beta_y <- (y-min(y))/(max(y)-min(y))
                                    data[,i] <- beta_y
                                    }}
 return(data)
}
