##' Create a matrix of parameter values for simulations
##'
##' Generates a matrix of parameters.  Each row are the population parameters
##' used to generate data from a single hypothetical "experiment."
##' @title Create Population Parameter Matrix
##' @param nexp number of parameter-value sets to generate (default 100000)
##' @param simparam.env an environment containing ranges of population parameters
##' @param firstseed initial seed
##' @param h0 logical value indicating if h0 is true or false
##' @param outfile name of save file for parameter-value matrix
##' @return a matrix of generated population parameter values
##' @export
createParamMx <-
    function(nexp=100000L,
             simparam.env=getParamRanges(),
             firstseed=NULL,
             h0=TRUE,
             outfile=NULL) { 
        if (!is.null(firstseed)) set.seed(firstseed);
        stopifnot(is.environment(simparam.env))
        assign("nexp",nexp,envir=simparam.env)
        assign("h0",as.logical(h0[1]),envir=simparam.env)
        param.mx <-
            with(simparam.env,
                 cbind(
                     int = runif(nexp, min=icept.range[1], max=icept.range[2]),
                     eff = ifelse(h0, slope["h0"], slope["h1"]),
                     err = runif(nexp, min=evar.range[1], max=evar.range[2]),
                     miss = runif(nexp, min=pmissing.range[1], max=pmissing.range[2]),
                     pMin = pMin,
                     pMax = pMax,
                     t00 = runif(nexp, min=t00.range[1], max=t00.range[2]),
                     t11 = runif(nexp, min=t11.range[1], max=t11.range[2]),
                     rsub = runif(nexp, min=r01.subj.range[1], max=r01.subj.range[2]),
                     w00 = runif(nexp, min=w00.range[1], max=w00.range[2]),
                     w11 = runif(nexp, min=w11.range[1], max=w11.range[2]),
                     ritm = runif(nexp, min=r01.item.range[1], max=r01.item.range[2]),
                     seed = sample.int(.Machine$integer.max-1L, nexp)
                     )
                 )
        if (is.null(outfile)) return(param.mx)
        save(param.mx, file=outfile)
        invisible(param.mx)
    }
