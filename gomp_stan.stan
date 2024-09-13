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
  b = exp(mu1 + var1 * beta);
  p = 1 - exp(-b / a .* (exp(a * dt) - 1));
}

model {
  // priors
  mu1 ~ normal(0, 1);
  beta ~ normal(0, 1);
  a_raw ~ normal(0, 1);
  
  // likelihood
  r ~ binomial(n, p);
}

generated quantities {
  int r_new[N];
  vector[N] p_new;
  
  p_new = 1 - exp(-b / a .* (exp(a * dt) - 1));
  for (i in 1:N) {
    r_new[i] = binomial_rng(n[i], p_new[i]);
  }
  real rate = exp(mu1);
}
