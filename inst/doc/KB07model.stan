data {
	int<lower=1> N;
	real RTtrunc[N];                     //outcome
    real<lower=-1,upper=1> S[N];	  //predictor
    real<lower=-1,upper=1> P[N];	  //predictor
    real<lower=-1,upper=1> C[N];	  //predictor
    real<lower=-1,upper=1> SP[N];	  //predictor
    real<lower=-1,upper=1> SC[N];	  //predictor
    real<lower=-1,upper=1> PC[N];	  //predictor
    real<lower=-1,upper=1> SPC[N];	  //predictor
	int<lower=1> I;                  //number of subjects
	int<lower=1> K;                  //number of items
	int<lower=1, upper=I> subj[N];   //subject id
	int<lower=1, upper=K> item[N];   //item id
}

parameters {
	vector[8] beta;			// intercept and slope
	real<lower=0> sigma_e;		// residual sd
	vector<lower=0>[8] sigma_u;	// subj sd
	vector<lower=0>[8] sigma_w;	// item sd
	cholesky_factor_corr[8] L_u;
	cholesky_factor_corr[8] L_w;
	matrix[8,I] z_u;
	matrix[8,K] z_w;
}

model {
	real mu[N]; 	// mu for likelihood
	matrix[I,8] u;	// random intercept and slopes subj
	matrix[K,8] w;
	
	# priors:
	beta ~ normal(0,2000);
	sigma_e ~ normal(0,500); // truncated normal, see parameters above
	sigma_u ~ normal(0,500);
	sigma_w ~ normal(0,500);
	L_u ~ lkj_corr_cholesky(4.0);
	L_w ~ lkj_corr_cholesky(4.0);
	to_vector(z_u) ~ normal(0,1);
	to_vector(z_w) ~ normal(0,1);
	
	u <- (diag_pre_multiply(sigma_u,L_u) * z_u)';	// subj random effects
	w <- (diag_pre_multiply(sigma_w,L_w) * z_w)';	// item random effects

	for (n in 1:N)
		mu[n] <- beta[1] + u[subj[n],1] + w[item[n],1] + 
		(beta[2] + u[subj[n],2] + w[item[n],2])*S[n] +
		(beta[3] + u[subj[n],3] + w[item[n],3])*P[n] +
		(beta[4] + u[subj[n],4] + w[item[n],4])*C[n] +
     	(beta[5] + u[subj[n],5] + w[item[n],5])*SP[n] +
        (beta[6] + u[subj[n],6] + w[item[n],6])*SC[n] +
        (beta[7] + u[subj[n],7] + w[item[n],7])*PC[n]+
        (beta[8] + u[subj[n],8] + w[item[n],8])*SPC[n];
	RTtrunc ~ normal(mu,sigma_e);        // likelihood
}

