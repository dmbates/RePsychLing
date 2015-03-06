##' Generate default parameter ranges for the simulations.
##'
##' This generates default parameter ranges.  The environment returned
##' can subsequently be modified by the user.
##' @title Get the parameter ranges used in the Monte Carlo simulations.
##' @return An environment containing variables and default ranges shown in the function sources.
##' @export
getParamRanges <- function() {
    param.env <- new.env()
    z3 <- c(0,3)
    m88 <- c(-0.8,0.8)
    with(param.env, {
        pmissing.range <- c(0,.05)      # proportion of missing data
        pMin <- 0              # lower bound on missing data rate
        pMax <- 0.8            # upper bound on missing data rate
        icept.range <- c(-3,3) # range of intercept value, continuous simulations
        slope <- c(h0=0, h1=.8) # the treatment effect (H0 true; H0 false    )
        evar.range <- z3        # range for error variance
        t00.range <- z3         # subject variance for the intercept
        t11.range <- z3         # subject variance for the slope
        r01.subj.range <- m88   # by-subject intercept/slope correlation
        w00.range <- z3         # by-item intercept variance
        w11.range <- z3         # by-item slope variance
        r01.item.range <- m88   # by-item intercept/slope correlation
    })
    param.env
}
