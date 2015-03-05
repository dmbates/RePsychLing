#' @export
prcomp.merMod <- function(x) {
    chfs <- getME(x,"Tlist")  # list of lower Cholesky factors
    nms <- names(chfs)
    unms <- unique(nms)
    names(unms) <- unms
    svals <- function(m) {
        vv <- svd(m,nv=0L)
        names(vv) <- c("sdev","rotation")
        vv$center <- FALSE
        vv$scale <- FALSE
        class(vv) <- "prcomp"
        vv
    }
    structure(lapply(unms,function(m) svals(bdiag(chfs[which(nms == m)]))),
              class="prcomplist")
}
#' @export
summary.prcomplist <- function(l) {
    lapply(l,summary)
}
