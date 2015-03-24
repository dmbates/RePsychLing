data {
  int<lower=1> N;
  real lrt[N];					   //outcome
  real<lower=-0.50,upper=0.50> sze[N];		   //predictor
  real<lower=-0.75,upper=0.25> spt[N];		   //predictor
  real<lower=-0.50,upper=0.50> obj[N];		   //predictor
  real<lower=-0.75,upper=0.25> grv[N];		   //predictor
  real<lower=-0.50,upper=0.50> orn[N];		   //predictor
  real<lower=-0.375,upper=0.375> spt_orn[N];	   //interaction
  real<lower=-0.250,upper=0.250> obj_orn[N];	   //interaction
  real<lower=-0.375,upper=0.375> grv_orn[N];	   //interaction
  real<lower=-0.375,upper=0.375> sze_spt[N];	   //interaction
  real<lower=-0.250,upper=0.250> sze_obj[N];	   //interaction
  real<lower=-0.375,upper=0.375> sze_grv[N];	   //interaction
  real<lower=-0.250,upper=0.250> sze_orn[N];	   //interaction
  real<lower=-0.1875,upper=0.1875> sze_spt_orn[N]; //interaction
  real<lower=-0.1250,upper=0.1250> sze_obj_orn[N]; //interaction
  real<lower=-0.1875,upper=0.1875> sze_grv_orn[N]; //interaction
  int<lower=1> I;		 //number of subjects
  int<lower=1, upper=I> subj[N]; //subject id
}

parameters {
  vector[16] beta;		// intercept and slopes
  real<lower=0> sigma_e;	// residual sd
  vector<lower=0>[8] sigma_u;	// subj sd
  cholesky_factor_corr[8] L_u;
  matrix[8,I] z_u;
}

model {
  real mu[N];			// mu for likelihood
  matrix[I,8] u;		// random intercept and slopes subj
	
  # priors:
  beta ~ normal(0,100);
  sigma_e ~ normal(0,50);
  sigma_u ~ normal(0,50);
  L_u ~ lkj_corr_cholesky(4.0);
  to_vector(z_u) ~ normal(0,1);
	
  u <- (diag_pre_multiply(sigma_u,L_u) * z_u)';	// subj random effects
	
  for (n in 1:N)
    mu[n] <-	beta[1] + u[subj[n],1] +
      beta[2] * sze[n] +						    
      (beta[3] + u[subj[n],2]) * spt[n] +
      (beta[4] + u[subj[n],3]) * obj[n] +
      (beta[5] + u[subj[n],4]) * grv[n] +
      (beta[6] + u[subj[n],5]) * orn[n] +
      (beta[7]) * sze_spt[n] +
      (beta[8]) * sze_obj[n] +		
      (beta[9]) * sze_grv[n] +
      (beta[10]) * sze_orn[n] +
      (beta[11]+ u[subj[n],6]) * spt_orn[n] +
      (beta[12]+ u[subj[n],7]) * obj_orn[n] +
      (beta[13]+ u[subj[n],8]) * grv_orn[n] +
      beta[14] * sze_spt_orn[n] +
      beta[15] * sze_obj_orn[n] +
      beta[16] * sze_grv_orn[n];
  lrt ~ normal(mu,sigma_e);
}
