## ----preliminaries,echo=FALSE,include=FALSE,cache=FALSE-----------------------------------
library(RePsychLing)
library(knitr)
library(rstan)
library(parallel)
library(xtable)
opts_chunk$set(comment=NA)
options(width=92,
        show.signif.stars = FALSE)

## ----kb07str------------------------------------------------------------------------------
str(kb07)

## ----X------------------------------------------------------------------------------------
X <- unname(model.matrix(~ 1+S+P+C+SP+SC+PC+SPC, kb07))
attr(X,"assign") <- NULL
str(X)

## ----XtX----------------------------------------------------------------------------------
crossprod(X)   # X'X

## ----standat------------------------------------------------------------------------------
standat <- '
data {
  int<lower=0>  N; // num observations
  int<lower=1>  K; // length of fixed-effects vector
  int<lower=0>  M; // num subjects
  int<lower=1>  J; // length of subj vector-valued random effects
  int<lower=0>  L; // num items
  int<lower=1>  I; // length of item vector-values random effects
  int<lower=1,upper=M> subj[N]; // subject indicator
  int<lower=1,upper=L> item[N]; // item indicator
  row_vector[K] X[N]; // model matrix for fixed-effects parameters
  row_vector[J] Zs[N]; // generator model matrix for subj random effects
  row_vector[I] Zi[N]; // generator model matrix for item random effects 
  vector[N]     y; // response vector (reaction time)
}

'

## ----stanpars-----------------------------------------------------------------------------
stanpars <- '
parameters {
  cholesky_factor_corr[J] Ls; // Cholesky factor of subj r.e. correlations
  cholesky_factor_corr[I] Li; // Cholesky factor of item r.e. correlations
  vector<lower=0>[J] taus; // standard deviations of unconditional subj r.e. dist
  vector<lower=0>[I] taui; // standard deviations of unconditional item r.e. dist
  vector[J] us[M];     // spherical subj random effects
  vector[I] ui[L];     // spherical item random effects
  vector[K] beta;      // fixed-effects
  real<lower=0> sigma; // standard deviation of response given random effects
}

'

## ----stantrans----------------------------------------------------------------------------
stantrans <- '
transformed parameters {
  matrix[J,J] corrs;
  matrix[I,I] corri;
  corrs <- tcrossprod(Ls);  // for monitoring subj correlations
  corri <- tcrossprod(Li);  // for monitoring item correlations
}

'

## ----model--------------------------------------------------------------------------------
stanmod <- '
model {
  matrix[J,J] Lambdas; 
  vector[J] bs[M];
  matrix[I,I] Lambdai; 
  vector[I] bi[L];
  taus ~ cauchy(0,2.5);
  taui ~ cauchy(0,2.5);
  Ls ~ lkj_corr_cholesky(2);
  Li ~ lkj_corr_cholesky(2);
  Lambdas <- diag_pre_multiply(taus,Ls);
  Lambdai <- diag_pre_multiply(taui,Li);
  for (m in 1:M) {
    us[m] ~ normal(0,1);
    bs[m] <- Lambdas * us[m];
  }
  for (l in 1:L) {
    ui[l] ~ normal(0,1);
    bi[l] <- Lambdai * ui[l];
  }
  for (n in 1:N)
    y[n] ~ normal(X[n] * beta + Zs[n] * bs[subj[n]] + Zi[n] * bi[item[n]], sigma);
}

'

## -----------------------------------------------------------------------------------------
model <- paste(standat, stanpars, stantrans, stanmod)

## ----maxdat-------------------------------------------------------------------------------
maxdat <- 
  within(list(), {
    N <- nrow(X)
    K <- J <- I <- ncol(X)
    M <- length(levels(kb07$subj))
    L <- length(levels(kb07$item))
    X <- Zs <- Zi <- unname(X)
    y <- kb07$RTtrunc
    subj <- as.integer(kb07$subj)
    item <- as.integer(kb07$item)
    }
    )
str(maxdat)

## ----maxmodel-----------------------------------------------------------------------------
maxmodel <- stan(model_name="maxmodel", model_code=model, data=maxdat, chains=0)

## ----KB07_stan,cache=TRUE-----------------------------------------------------------------
system.time(KB07_stan <-
  sflist2stanfit(
    mclapply(1:4, mc.cores = 4,    # adjust number of cores to suit 
      function(i) stan(fit = maxmodel, 
                       data = maxdat,
                       iter=2000,
                       chains = 1, 
                       chain_id = i, 
                       refresh = -1))
    )
  )

## ----KB07_results,cache=FALSE-------------------------------------------------------------
KB07_results<- summary(KB07_stan,
                       pars=c("beta", "sigma",
                              "taus","taui",
                              "corrs","corri"),
                       probs = c(0.025,  0.975), digits_summary = 3)
rownames(KB07_results$summary)

## ----printmaxmodel,echo=FALSE,eval=TRUE,cache=FALSE,results="asis"------------------------
print(xtable(KB07_results$summary[c(1:25,27:33,36:41,45:49,55:57,62:64,81,91:97,100:105),c(1,4,5)]), type="html")

## ----datreduced---------------------------------------------------------------------------
finaldat <- 
  within(list(), {
    N <- nrow(X)
    K <- ncol(X)
    J <- 1L
    I <- 2L
    M <- length(levels(kb07$subj))
    L <- length(levels(kb07$item))
    X <- X
    Zs <- X[, 1, drop=FALSE]
    Zi <- X[, 1:2]
    y <- kb07$RTtrunc
    subj <- as.integer(kb07$subj)
    item <- as.integer(kb07$item)
    }
    )
str(finaldat)

## -----------------------------------------------------------------------------------------
system.time(KB07_finalstan <-
  sflist2stanfit(
    mclapply(1:4, mc.cores = 4,    # adjust number of cores to suit 
      function(i) stan(fit = maxmodel, 
                       data = finaldat,
                       iter=2000,
                       chains = 1, 
                       chain_id = i, 
                       refresh = -1))
    )
  )

