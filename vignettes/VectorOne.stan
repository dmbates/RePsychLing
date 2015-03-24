data {
  int<lower=0>  N; // num observations
  int<lower=1>  K; // length of fixed-effects vector
  int<lower=1>  M; // num subjects
  int<lower=1>  J; // length of individual vector-valued random effects
  int<lower=1,upper=M> subj[N]; // subject indicator
  row_vector[K] X[N]; // model matrix for fixed-effects parameters
  row_vector[J] Z[N]; // generator model matrix for random-effects
  vector[N]     y; // response vector (reaction time)
}

parameters {
  cholesky_factor_corr[J] L; // Cholesky factor of corr in uncond r.e. dist
  vector<lower=0>[J] tau; // standard deviations of unconditional r.e. dist
  vector[K] beta;      // fixed-effects
  real<lower=0> sigma; // standard deviation of response given random effects
  vector[J] u[M];      // spherical random effects
}

transformed parameters {
  matrix[J,J] corr;
  corr <- tcrossprod(L);  // for monitoring the correlations
}

model {
  matrix[J,J] Lambda; 
  vector[J] b[M];
  tau ~ cauchy(0,2.5);
  L ~ lkj_corr_cholesky(2);
  Lambda <- diag_pre_multiply(tau,L);
  for (m in 1:M) {
    u[m] ~ normal(0,1);
    b[m] <- Lambda * u[m];
  }
  for (n in 1:N)
    y[n] ~ normal(X[n] * beta + Z[n] * b[subj[n]], sigma);
}
