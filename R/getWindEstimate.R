#' @importFrom stats optim
NULL
#' Estimate wind speed from a sample of ground speeds
#'
#' @param groundSpeeds matrix with two columns representing the ground speeds. 
#' @param phi numeric of length one giving the auto correlation. 
#' @param windStart numeric of length 2 giving the wind speed where to optimize from. 
#'
#' @return an list with parameter estimates
#' @export
#' @docType methods
#' @rdname getWindEstimate-methods
#'
#' @examples
#' s<-seq(0,2*pi, .1)
#' set.seed(34)
#' getWindEstimate(cbind(4*cos(s)+3+rnorm(length(s)), 4*sin(s)+2+rnorm(length(s))),0)
#' getWindEstimate(cbind(4*cos(s)+3+rnorm(length(s),sd=.2), 4*sin(s)+2+rnorm(length(s),sd=.2)),0)

setGeneric("getWindEstimate", function(groundSpeeds, phi, windStart=c(0,0)){
  standardGeneric("getWindEstimate")
  })

#' @rdname getWindEstimate-methods
#' @aliases getWindEstimate,matrix,numeric,ANY-method
setMethod("getWindEstimate",signature = signature(groundSpeeds="matrix", phi="numeric",windStart="ANY"), function(groundSpeeds, phi, windStart) {
  stopifnot(length(phi)==1)
  stopifnot(length(windStart)==2)
  optimResult = optim(windStart, calcVar, hessian = T, phi=phi, groundSpeeds=groundSpeeds)
  n = nrow(groundSpeeds)
  if (optimResult$convergence == 0) {
    if (det(optimResult$hessian) != 0)
      covar = solve(optimResult$hessian / optimResult$value * n / 2)
    else
      covar = NULL
    list(
      windEst = optimResult$par, residualVarAirSpeed = optimResult$value, covar =
        covar
    )
  }
  else
    NULL
}
)