library(tidyverse)
library(flexsurv)

## Generate data
exposure <- sample(0:1, 100, replace = TRUE)
times <- rexp(100, 10*exposure + 5 * (1-exposure))
mydf <- tibble(event = 1, times = times, exposure)

## For no covariates
## Run model
a <- flexsurvreg(Surv(times, event) ~ 1, data = mydf, 
            dist = "gengamma")
## For no covariates get mean using summary function
b <- summary(a, type = "mean")
## For no covariates get mean using mean_gengamma function
mean_gengamma(mu = a$coefficients[1], sigma = exp(a$coefficients[2]), Q = a$coefficients[3])
# note that the central estimate is exactly the same as this

## For one covariate
## run model
a2 <- flexsurvreg(Surv(times, event) ~ exposure, data = mydf, 
                  dist = "gengamma")
a2
## For no covariates get mean using summary function
b2 <- summary(a2, type = "mean", newdata = tibble(exposure = c(0,1)))
b2
## For no covariates get mean using mean_gengamma function
## Get mean for unexposed
mean_gengamma(mu = a2$coefficients[1], sigma = exp(a2$coefficients[2]), Q = a2$coefficients[3])
## get mean for exposed
mean_gengamma(mu = a2$coefficients[1] + a2$coefficients[4], sigma = exp(a2$coefficients[2]), Q = a2$coefficients[3])
## note is same answer


## get confidence interval by sampling from multivariate normal

## Initially for simple model wihtout covariates
smpls_coefs <- mvtnorm::rmvnorm(20000, 
                                mean = c(mu = a$res.t["mu", "est"], 
                                         sigma =  a$res.t ["sigma", "est"],
                                         Q = a$res.t["Q", "est"]),
                                sigma = vcov(a))
# exponentiate as per non-confidence interval analysis
## Compare to summary result from package
b

MakeCIMean <- function(mymatrix){
  res <- map_dbl(1:nrow(mymatrix), function(i) mean_gengamma(mu = mymatrix[i,"mu"], sigma = exp(mymatrix[i,"sigma"]), Q = mymatrix[i,"Q"]))
  res_smry = quantile(res, probs = c(0.025, 0.975))
  names(res_smry) <- c("lci", "uci")
  res_smry
}
MakeCIMean(smpls_coefs)

## Next for model with covariates (eg male and female)
smpls_coefs_withcov <- mvtnorm::rmvnorm(20000, 
                                mean = a2$res.t[,"est"],
                                sigma = vcov(a2))
head(smpls_coefs_withcov)
## Create exposures
smpls1 <- smpls_coefs_withcov
smpls1[, "mu"] <- smpls1[, "mu"] + smpls1[, "exposure"]
smpls1 <- smpls1[ , c("mu", "sigma", "Q")]
smpls2 <- smpls_coefs_withcov[ , c("mu", "sigma", "Q")]
b2
MakeCIMean(smpls1)
MakeCIMean(smpls2)
