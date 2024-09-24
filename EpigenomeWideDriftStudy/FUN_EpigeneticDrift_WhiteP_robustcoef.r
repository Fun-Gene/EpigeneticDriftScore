library("lmtest")
library("sandwich")
library(dplyr)
library(tidyr)

#
getWhitePdriftCpG <- function(method, alpha, nSamples, nIterations, age, cpgdata, covdata){
                    whiteP <-  betaX <- seX <- tX <- pX <- betaX2 <- seX2 <- tX2 <- pX2 <- betaInter <- seInter <- tInter <- pInter <- rep(NA,nIterations)
                    set.seed(111)
                    for(i in 1:nIterations){
			　　　　x <- age
                               y <- cpgdata[,i]
                               fo <- as.formula(paste0("y~x+", paste0(colnames(covdata), collapse="+")))
			　　　　model1 <- lm(fo, data=covdata)
			　　　if (method == "white") {model2 <- lm((model1$residuals)^2 ~ x + I(x^2))
                                 # Calculate White test heteroskedasticity
                                 r2  <- summary(model2)$r.squared
                                 # White test statistic
                                 stat <- nSamples*r2
                                 # Calculate p value
                                 whiteP[i] <- pchisq(q=stat, df=2, lower.tail=F)
                                 names(whiteP)[i] <- names(cpgdata)[i]
                                # capture output from model summary
                                 robustcoef <- coeftest(model2, vcov = vcovHC(model2, type = "HC0"))
                                 betaX[i] <- robustcoef['x', 'Estimate']
                                 seX[i] <- robustcoef['x', 'Std. Error']
                                 tX[i] <- robustcoef['x', 't value']
                                 pX[i] <- robustcoef['x', 'Pr(>|t|)']
                                 betaX2[i] <- robustcoef['I(x^2)', 'Estimate']
                                 seX2[i] <- robustcoef['I(x^2)', 'Std. Error']
                                 tX2[i] <- robustcoef['I(x^2)', 't value']
                                 pX2[i] <- robustcoef['I(x^2)', 'Pr(>|t|)']
                                 betaInter[i] <- robustcoef['(Intercept)', 'Estimate']
                                 seInter[i] <- robustcoef['(Intercept)', 'Std. Error']
                                 tInter[i] <- robustcoef['(Intercept)', 't value']
                                 pInter[i] <- robustcoef['(Intercept)', 'Pr(>|t|)']



                               }
                               }
 power <- length(whiteP[whiteP<alpha])/length(whiteP)
 return(list(power=power,whiteP=whiteP,betaX=betaX,seX=seX,tX=tX,pX=pX,betaX2=betaX2,seX2=seX2,tX2=tX2,pX2=pX2,betaInter=betaInter,seInter=seInter,tInter=tInter,pInter=pInter))
}
