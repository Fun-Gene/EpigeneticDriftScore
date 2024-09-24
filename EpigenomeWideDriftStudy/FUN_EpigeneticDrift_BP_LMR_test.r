library(tidyverse)
library(parallel)
library(gamlss)
library(splines)

Cited from this link: https://github.com/JacobBergstedt/MIMETH/blob/main/Scripts/R_Scripts/infer_age_dispersion_include_cells_in_variance_884_samples.R

getBPdriftCpG_fit_age_dispersion <- function(y,fm, mf) {
						  mf$y <- y
						  mf <- na.omit(mf)
				# 1. null
				  null <- gamlss(fm,
                  data = mf,
                  control = gamlss.control(trace = FALSE))
				# 2. add age
				  alt <- update(null, formula. = ~ Age, what = "sigma")
				# 3. LRT
				  LRT <- deviance(null) - deviance(alt)
				  p_lrt <- pchisq(LRT, df = 1, lower.tail = FALSE)
				  est <- coef(alt, what = "sigma")["Age"]
				  tibble(Term = "Age",
						 Estimate = est,
						 P_value = p_lrt)
}

