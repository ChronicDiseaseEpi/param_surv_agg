library(tidyverse)
library(rjags)
library(runjags)

## Explore how Gompertz is specified in flexsurv ----
MyExp <- function(rate, fu) {
  1 - exp(-rate*fu)
}
MyGmp <- function(rate, shape, fu) {
  1 - exp( (-rate/shape) * (exp(shape*fu) - 1) )
}
MyGmp2 <- function(rate, shape, fu) {
  if_else(shape == 0,
          1 - exp(-rate*fu),
          1 - exp( (-rate/shape) * (exp(shape*fu) - 1) ))
}

## specify shape and rate on non log scales
x <- expand.grid(rate = c(0.1, 1),
                 shape = c(1, 1.001, 1.5)) %>% 
  as_tibble() %>% 
  mutate(fu = 1.5)
x <- x %>% 
  mutate(
    r_exp = stats::pexp(fu, rate),
    my_exp = MyExp(rate, fu),
    fs_gmp = flexsurv::pgompertz(q = fu,
                                 shape = log(shape), 
                                 rate = rate),
    my_gmp = MyGmp(rate = rate, shape = log(shape), fu = fu),
    my_gmp2 = MyGmp2(rate, shape = log(shape), fu))
x

## generate some aggregate data based on Gompertz ----
myrate <- 1
fu <- rexp(10, 1)
mya <- 0.5
ps <- flexsurv::pgompertz(fu, shape = mya, rate = myrate)
ps_exp <- pexp(fu, rate = myrate)
ps
ps_exp
Ns <- runif(10, 100, 1000) %>% round()
r <- rbinom(10, size = Ns, prob = ps)

## Run JAGS exponential model ----
jagscode_exp <- "model {
  for (i in 1:N) { # N=number of data points in dataset

    # likelihood
    r[i] ~ dbin(p[i], n[i])
    p[i] <- 1 - exp(-b[i]*dt[i])
    # fixed effects model
    log(b[i]) <-  mu1 
  }
  # priors
  mu1 ~ dnorm(0, 1)
  for (i in 1:N) { # N=number of data points in dataset
      r_new[i] ~ dbin(p[i], n[i])
      p_new[i] <- 1 - exp(-b[i]*dt[i])
  }
  rate <- exp(mu1)
}"
test1 <- run.jags(model = jagscode_exp, 
                  monitor = c("mu1", "rate", "p_new"),
                  data = list(r = r, 
                              dt = fu,
                              n = Ns,
                              N = length(Ns)), 
                  n.chains = 2, sample = 2000)

## Run JAGS Gompertz model -----
jagscode_gomp <- "model {
  for (i in 1:N) { # N=number of data points in dataset
    
    # likelihood
    r[i] ~ dbin(p[i], n[i])
    p[i] <- 1 - exp( -b[i]/a * (exp(a*dt[i]) - 1) )
    # fixed effects model
    log(b[i]) <-  mu1 
  }
  # priors
  mu1 ~ dnorm(0, 1)
  a <- ifelse(a_raw == 0, 0.0001, a_raw)
  a_raw ~ dnorm(0, 1)
  
  for (i in 1:N) { # N=number of data points in dataset
    r_new[i] ~ dbin(p[i], n[i])
    p_new[i] <-  1 - exp( -b[i]/a * (exp(a*dt[i]) - 1) )
  }
  rate <- exp(mu1)

}"

test2 <- run.jags(model = jagscode_gomp, 
                  monitor = c("mu1", "rate",  "a", "p_new"),
                  data = list(r = r, 
                              dt = fu,
                              n = Ns,
                              N = length(Ns)), 
                  n.chains = 2, sample = 2000)

## Pull and compare results -----
as_set <- c(log(myrate), myrate, mya, ps)
as_set <- tibble(
  param = c("mu1",
            "rate",
            "a",
            paste0("p_new[", 1:10, "]")),
  est = as_set,
  mdl = "known")

tst1 <- test1$summaries %>% 
  as_tibble(rownames = "param") %>% 
  select(param, est = Mean) %>% 
  mutate(mdl = "exponential")

tst2 <- test2$summaries %>% 
  as_tibble(rownames = "param") %>% 
  select(param, est = Mean) %>% 
  mutate(mdl = "gompertz")

cmpr <- bind_rows(as_set,
                  tst1,
                  tst2)
cmpr %>% 
  spread(mdl, est)

