
data {
  int<lower=1> T;  // Time points
  vector[T-1] tau;   // Width between time-points
  int y[T];        // Events
  vector[T] n;     // At risk
}

parameters {
  real beta_01;              // Initial coeff1
  real<lower=0> Z1;          // Variance coeff1
  real<lower=0> Z2;          // Variance coeff2
  real beta_02;              // Initial coeff2
  vector[T] zeta_tilde1;      // Tranformation of zeta2 (as in 8-schools example)
  vector[T-1] zeta_tilde2;      // Tranformation of zeta2 (as in 8-schools example)
  real<lower=0.7, upper=0.999> phi;  // Whether or not we know phi values (don't let = 1 as messes up extrap calcs)
}

transformed parameters {
    vector[T] beta1;         // State 1
    vector[T-1] beta2;         // State 2
    { // Don't want to save this
vector[T] zeta1;         // Innovations
vector[T-1] zeta2;         // Innovations

zeta1 = sqrt(Z1) * zeta_tilde1;
zeta2 = sqrt(Z2) * zeta_tilde2;
beta1[1] = beta_01 + zeta1[1];;
beta2[1] = beta_02 + zeta2[1];
for (t in 2:T-1) {
  beta1[t] = beta1[t-1] + beta2[t-1] * phi * tau[t-1] + zeta1[t];
  beta2[t] = beta2[t-1] * phi + zeta2[t];
}
beta1[T] = beta1[T-1] + beta2[T-1] * phi * tau[T-1] + zeta1[T];
}
}

model {
  Z1 ~ inv_gamma(1, 0.005);
  Z2 ~ inv_gamma(1, 0.005);
  zeta_tilde1 ~ normal(0, 1);
  zeta_tilde2 ~ normal(0, 1);
  y ~ poisson(exp(beta1) .* n);
}

generated quantities{
  real level;
  real trend;
  
  level = beta1[T];
  trend = beta2[T-1];
}
