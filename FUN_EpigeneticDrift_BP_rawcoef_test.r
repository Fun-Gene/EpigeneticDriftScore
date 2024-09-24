getBPdriftCpG_rawcoef <- function(method, cpg, nIterations, alpha){
                    pVec <- estVec <- seVec <- rep(NA,nIterations)
                    set.seed(111)
                    for(i in 1:nIterations){
							   x <- cpg[,1]
                               y <- cpg[, c(1+i)]
                               model1 <- lm(y ~ x)
							   if (method == "BP") { model2 <- lm((model1$residuals)^2 ~ x)}
                              # capture output from model summary
                               estVec[i] <- coef(summary(model2))['x', 'Estimate']
                               seVec[i] <- coef(summary(model2))['x', 'Std. Error']
                               pVec[i] <- coef(summary(model2))['x', 'Pr(>|t|)']
                               }
 power <- length(pVec[pVec<alpha])/length(pVec)
 return(list(power=power,p=pVec,beta=estVec,se=seVec))
}