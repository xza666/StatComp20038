# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A Metropolis sampler using Rcpp
#' @description A Metropolis sampler using Rcpp
#' @param sigma Variances  used for the proposal distribution
#' @param x0 Initial value
#' @param N the number of samples
#' @return a random sample of size \code{n} and the rejection rate
#' @examples
#' \dontrun{
#' rw=rwMC(2,25,1000)
#' plot(rw$x,type="l",lwd=2,xlab="sigma=2",ylab="X")
#' }
#' @export
rwMC <- function(sigma, x0, N) {
    .Call('_StatComp20038_rwMC', PACKAGE = 'StatComp20038', sigma, x0, N)
}

