library(tidyverse)
library(rstan)

## Used copilot to translate BUGS to Stan. Appears to recover values

## generate some aggregate data based on Gompertz ----
myvar = sample(0:1, 10, replace = TRUE)
myrate <- 0.2
mybeta <- 2
myrates <- if_else(myvar ==1, myrate*mybeta, myrate) 
fu <- rexp(10, 1)
mya <- 0.5
ps <- flexsurv::pgompertz(fu, shape = mya, rate = myrates)
ps_exp <- pexp(fu, rate = myrate)
ps
ps_exp
Ns <- runif(10, 100, 1000) %>% round()
r <- rbinom(10, size = Ns, prob = ps)

mdl <- rstan::stan_model("gomp_stan.stan")

fit <- sampling(mdl, data = list(r = r, 
                                 dt = fu,
                                 n = Ns,
                                 N = length(Ns),
                                 var1 = myvar))
fit

