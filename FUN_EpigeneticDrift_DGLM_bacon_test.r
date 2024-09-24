#Load all necessary libraries.

library(dglm)
library(tidyverse)
library(bacon)


getBPdriftCpG_runDGLM <- function(x,cpgnames)
{
  result_dglm <- matrix(nrow=length(cpgnames), ncol=9)
  colnames(result_dglm)<-c("cgID","effectsize.linear","std.error.linear","t.linear","p.linear","effectsize.disp","std.error.disp","t.disp","p.disp")
  for(CpG in 1:length(cpgnames))
  {
    print(CpG) #to check progress
    # 1. Set CpG
    probe <- cpgnames[CpG]
    # 2. set data.frame
    tempdata <- data.frame(Age=x[,"Age"], data = x[,probe])
    # 3. fit
    fit <- tryCatch({ dglm(formula = data ~ Age , dformula = ~ Age, data=tempdata, family = gaussian(link = identity))},error=identity)
    # 4.statistics
    if(is.null(fit$message)){ #good
    # 4.1 clock  #extract info about the linear model
      effectsize.linear <- summary(fit)$coefficients[2,1]
      std.error.linear<-summary(fit)$coefficients[2,2]
      t.value.linear <-summary(fit)$coefficients[2,3]
      p.value.linear <- summary(fit)$coefficients[2,4]

     # 4.2 drift  #extract info about the dispersion model
      effectsize.disp <- summary(fit$dispersion.fit)$coefficients[2]
      std.error.disp<-summary(fit$dispersion.fit)$coefficients[4]
      t.value.disp <- summary(fit$dispersion.fit)$coefficients[6]
      p.value.disp <- summary(fit$dispersion.fit)$coefficients[8]

      # 4.3 Summary statistics Output
      out <- matrix(c(probe,effectsize.linear,std.error.linear,t.value.linear,p.value.linear,effectsize.disp,std.error.disp,t.value.disp,p.value.disp), ncol=9)
      result_dglm[CpG,] <- out
    }else{                  #bad
      print("Didn't converge")
      out <- matrix(c(probe, fit$message, rep(NA, times = 7)), ncol=9)
      result_dglm[CpG,] <- out
    }

  }
  return(result_dglm)
}

DGLM_converge <- runDGLM(data)
set.seed(1)
# bacon 
bc_DGLM<- bacon(DGLM_converge[,3])

# Calculate bias and inflation
bias(bc_DGLM)       

inflation(bc_DGLM)    


#Extract the BACON-adjusted t-statistics and p-values
tstats_DGLM<-tstat(bc_DGLM)
pvals_DGLM <-pval(bc_DGLM)
padjs_DGLM <- as.matrix(p.adjust(pvals_DGLM, method="bonf"))

# Extract BACON-adjusted effectsize.
set.seed(1)
bc_DGLM_es<-bacon(NULL,DGLM_converge[,1],DGLM_converge[,2])
es_DGLM<-es(bc_DGLM_es)

# Extract BACON-adjusted standard error.
set.seed(1)
bc_DGLM_se<-bacon(NULL,DGLM_converge[,1],DGLM_converge[,2])
se_DGLM<-se(bc_DGLM_se)

