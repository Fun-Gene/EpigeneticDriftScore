getWhitePdriftCpG_rawcoef <- function(method, alpha, nSamples, nIterations, age, cpgdata, covdata){
                    whiteP <-  estVec_x <- seVec_x <- pVec_x <- estVec_x2 <- seVec_x2 <- pVec_x2 <- rep(NA,nIterations)

                    set.seed(111)
                    for(i in 1:nIterations){
							                 x <- age
                               y <- cpgdata[,i]
                               fo <- as.formula(paste0("y~x+", paste0(colnames(covdata), collapse="+")))
							                 model1 <- lm(fo, data=covdata)
							                 if (method == "travel") {model2 <- lm((model1$residuals)^2 ~ x + I(x^2))
                                 # Calculate White test heteroskedasticity
                                 r2  <- summary(model2)$r.squared
                                 # White test statistic
                                 stat <- nSamples*r2
                                 # Calculate p value
                                 whiteP[i] <- pchisq(q=stat, df=2, lower.tail=F)
                                 names(whiteP)[i] <- names(cpgdata)[i]
                                # capture output from model summary
                                 coef <- summary(model2)$coef
                                 estVec_x[i] <- coef['x', 'Estimate']
                                 seVec_x[i] <- coef['x', 'Std. Error']
                                 pVec_x[i] <- coef['x', 'Pr(>|t|)']
                                 estVec_x2[i] <- coef['I(x^2)', 'Estimate']
                                 seVec_x2[i] <- coef['I(x^2)', 'Std. Error']
                                 pVec_x2[i] <- coef['I(x^2)', 'Pr(>|t|)']
                               }
                               }
 power <- length(whiteP[whiteP<alpha])/length(whiteP)
 return(list(power=power,whiteP=whiteP,betaX=estVec_x,seX=seVec_x,pX=pVec_x,betaX2=estVec_x2,se2=seVec_x2,p2=pVec_x2))
}
