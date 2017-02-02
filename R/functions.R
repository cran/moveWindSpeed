#' A function to generate an isThermallingFunction
#'
#' @param totalAngle the cumulative angle that is required to consider an trajectory thermalling
#' @param minMeanSpeed the minimal air speed that is required to decide of a track is thermalling
#'
#' @return a function is returned that based on a series of headings returns a logical value to indicate is a track is thermalling or not
#' @export
#'
#' @examples
#' fun<-getDefaultIsThermallingFunction(170)
#' fun(1:160)
#' fun(1:190, rep(2,190))
#' fun<-getDefaultIsThermallingFunction(170, 3)
#' fun(1:190, rep(2,190))
#' fun(1:190, rep(3.4,190))

getDefaultIsThermallingFunction <- function(totalAngle = 360, minMeanSpeed=NULL) {
  stopifnot(is.numeric(totalAngle))
  stopifnot(length(totalAngle)==1)
if(is.null(minMeanSpeed)){
  return(  function(headings, speeds) {
    diffHeading = diff(headings)
    diffHeading <- ((diffHeading + 180) %% 360) - 180
    max(abs(cumsum(diffHeading))) >= totalAngle
  })}else{
    stopifnot(is.numeric(minMeanSpeed))
    stopifnot(length(minMeanSpeed)==1)
    
    return(  function(headings, speeds) {
      diffHeading = diff(headings)
      diffHeading <- ((diffHeading + 180) %% 360) - 180
      max(abs(cumsum(diffHeading))) >= totalAngle & mean(speeds)>minMeanSpeed
    })
    }
}
#' Estimate the log likelihood 
#'
#' @param sigma the residual variance in airspeed
#' @param phi the autocorrelation used in the calculations
#'
#' @return the log likelihood
#' @export
#'
#' @examples
#' windEstimLogLik(c(1.3,.6,1.5,1.8),.3)
#' windEstimLogLik(c(1.3,.6,1.5,1.8),.5)

windEstimLogLik<-function(sigma, phi){
  sum(length(sigma) * log(sigma ^ 2) - log(1 - phi ^ 2))
  }
