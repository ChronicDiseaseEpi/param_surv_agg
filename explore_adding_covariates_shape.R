library(tidyverse)
library(rjags)
library(runjags)

# seems to underestimate zeta. Need to check which is best scale for modelling zeta

## generate some aggregate data based on Gompertz ----
myvar = sample(0:1, 10, replace = TRUE)
myrate <- 0.2
mybeta <- 2
myzeta <- 2
myrates  <- if_else(myvar ==1, myrate*mybeta, myrate) 
myshapes <- if_else(myvar ==1, mya*myzeta, mya) 
fu <- rexp(10, 1)
mya <- 0.5
ps <- flexsurv::pgompertz(fu, shape = myshapes, rate = myrates)
ps_exp <- pexp(fu, rate = myrate)
ps
ps_exp
Ns <- runif(10, 100, 1000) %>% round()
r <- rbinom(10, size = Ns, prob = ps)
r
## Run JAGS Gompertz model -----
jagscode_gomp <- "model {
  for (i in 1:N) { # N=number of data points in dataset
    
    # likelihood
    r[i] ~ dbin(p[i], n[i])
    p[i] <- 1 - exp( -b[i]/a[i] * (exp(a[i]*dt[i]) - 1) )
    # fixed effects model
    log(b[i]) <-  mu1 + var1[i]*beta
    a_lp[i] <- mu2 + var1[i]*zeta
  }
  # transformation does not work if a is zero so replace with a number
  # close to zero; assuming ifelse is vectorised
  a <- ifelse(a_lp == 0, 0.0001, a_lp)
  # priors
  mu1 ~ dnorm(0, 1)
  beta ~ dnorm(0,1)
  mu2 ~ dnorm(0, 1)
  zeta ~ dnorm(0, 1)
  
  for (i in 1:N) { # N=number of data points in dataset
    r_new[i] ~ dbin(p[i], n[i])
    p_new[i] <-  1 - exp( -b[i]/a[i] * (exp(a[i]*dt[i]) - 1) )
  }
  rate <- exp(mu1)

}"

test1 <- run.jags(model = jagscode_gomp, 
                  monitor = c("mu1", "rate", "beta", "p_new", "zeta"),
                  data = list(r = r, 
                              dt = fu,
                              n = Ns,
                              N = length(Ns),
                              var1 = myvar), 
                  n.chains = 2, sample = 2000)


## Pull and compare results -----
as_set <- c(log(myrate), myrate, log(mybeta), myzeta, ps)
as_set <- tibble(
  param = c("mu1",
            "rate",
            "beta",
            "zeta",
            paste0("p_new[", 1:10, "]")),
  est = as_set,
  mdl = "known")

tst1 <- test1$summaries %>% 
  as_tibble(rownames = "param") %>% 
  select(param, est = Mean) %>% 
  mutate(mdl = "gompertz")

cmpr <- bind_rows(as_set,
                  tst1)
cmpr %>% 
  spread(mdl, est)

