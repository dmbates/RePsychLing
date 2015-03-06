##' Upper Cholesky factor of a 2 by 2 covariance matrix
##' 
##' @title Upper Cholesky factor of covariance
##' @param v1 positive numeric, variance of first coordinate
##' @param v2 positive numeric, variance of second coordinate
##' @param r numeric in the range [-1,1], correlation
##' @return The upper triangular Cholesky factor of the covariance matrix
mkCovR <- function(v1, v2, r)
    matrix(c(sqrt(v1),0,sqrt(v2)*c(r,sqrt(1-r*r))), ncol = 2L)

##' Create simulated data given a row of parameters and the no. of subjects and items
##' @title Return a dataframe with simulated data given a set of population parameters
##' @param nsubj number of subjects (default is 24)
##' @param nitem number of items (default is 24)
##' @param wsbi logical, is design between items (TRUE) or within items (default FALSE).
##' @param mcr.params vector of parameters for generation of data (see createParamMx)
##' @param missMeth method of generating missing data (default "random")
##' @param rigen logical, use a random-intercepts-only generative model (default FALSE)
##' @return a data frame of simulated data
##' @export
mkDf <- function(nsubj=24L, nitem=24L, wsbi=FALSE,
                 mcr.params=createParamMx(1L,firstseed=sample.int(1000000L,1L))[1L,],
                 missMeth=c("random", "none", "randomBig", "bycond", "bysubj", "bysubjcond"),
                 rigen=FALSE) {
    stopifnot(nitem %% 2L == 0L, nitem > 3L, # both nitem and nsubj must be even and > 3
              nsubj %% 2L == 0L, nsubj > 3L,
              !is.null(names(mcr.params)))
                                        # create subj and item factors with meaningful levels
    subjf <- factor(sprintf(paste("S%0",ceiling(log10(nsubj)),"d",sep=''),1:nsubj))
    itemf <- factor(sprintf(paste("I%0",ceiling(log10(nitem)),"d",sep=''),1:nitem))
    ans <- expand.grid(item=itemf,subj=subjf)
                                        # create the covariate according to wsbi
    cond <- rep.int(LETTERS[1:2], nitem %/% 2L) # condition factor levels within first subj
    if (wsbi) {
        cond <- c(cond,cond)            # levels of cond repeat for all subjects
    } else {
        cond <- c(cond,rev(cond))       # levels of cond reverse for succesive subjects
    }
    ans$cond <- factor(rep.int(cond, nsubj %/% 2L))
    contrasts(ans$cond) <- contr.sum    # use sum contrasts for (-1,+1) encoding
    attr(ans, "out.attrs") <- NULL      # remove some attributes added by expand.grid
    mm <- model.matrix(~ cond, ans)     # fixed-effects model matrix

    pars <- as.list(mcr.params)         # a list is easier to access by name than a named vector
    set.seed(pars$seed)
                                        #  matrices of random effects
    subjre <- matrix(rnorm(nsubj * 2L), ncol=2L, # subject random effects
                     dimnames=list(subjf,NULL)) %*% with(pars, mkCovR(t00, t11, rsub))
    itemre <- matrix(rnorm(nitem * 2L), ncol=2L, # item random effects
                     dimnames=list(subjf,NULL)) %*% with(pars, mkCovR(w00, w11, ritm))
    colnames(subjre) <- colnames(itemre) <- colnames(mm)
                                        # create the linear predictor
    linpred <- as.vector(mm %*% with(pars, c(int,eff))) # fixed-effects contribution
    if (rigen) {                        # random intercept only for subjects
        linpred <- linpred + subjre[ans$subj,1L]
    } else {
        linpred <- linpred + rowSums(subjre[ans$subj,] * mm)
    }
    if (wsbi) {                         # random intercept for items
        linpred <- linpred + itemre[ans$item,1L]
    } else {
        linpred <- linpred + rowSums(itemre[ans$item,] * mm)
    }
                                        # establish missing data values
    n <- nrow(mm)
    linpred[switch(
        match.arg(missMeth),
        random = sample.int(n,round(n * pars$miss)),
        none = integer(0),
        randomBig = which(runif(n) < runif(1L,pars$pMin,pars$pMax)),
        bycond = which(runif(n) < runif(2,pars$pMin,pars$pMax)[ans$cond]),
        bysubj = which(runif(n) < runif(nsubj,pars$pMin,pars$pMax)[ans$subj]),
        bysubjcond = {
            rates <- matrix(runif(2*nsubj, pars$pMin, pars$pMax), ncol=2L)
            which(runif(n) < rates[cbind(ans$subj,ans$cond)])
        })] <- NA
    ans$resp <- linpred + sqrt(pars$err) * rnorm(n)
    ans
}
