data {
  int<lower=0> N; // number of data points
  int<lower=0> r[N]; // observed counts
  int<lower=0> n[N]; // number of trials
  vector[N] dt; // time intervals
  vector[N] var1; // covariate
}

parameters {
  real mu1; // intercept
  real beta; // coefficient
  real a_raw; // raw parameter for 'a'
}

transformed parameters {
  real a;
  vector[N] b;
  vector[N] p;
  
  a = a_raw == 0 ? 0.0001 : a_raw; // ensure 'a' is not exactly zero
  for (i in 1:N) {
    b[i] = exp(mu1 + var1[i] * beta);
    p[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
  }
}

model {
  // priors
  mu1 ~ normal(0, 1);
  beta ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  
  // likelihood
  for (i in 1:N) {
    r[i] ~ binomial(n[i], p[i]);
  }
}

generated quantities {
  vector[N] r_new;
  vector[N] p_new;
  
  for (i in 1:N) {
    p_new[i] = 1 - exp(-b[i] / a * (exp(a * dt[i]) - 1));
    r_new[i] = binomial_rng(n[i], p_new[i]);
  }
  real rate = exp(mu1);
}
