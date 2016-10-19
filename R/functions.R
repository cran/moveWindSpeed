
#' getDefaultIsThermallingFunction
#'
#' @param totalAngle the cumilative angle that is required to consider an trajectory thermaling
#'
#' @return a function is returned that based on a series of headings returns a logical value to indicate is a track is thermaling or not
#' @export
#'
#' @examples
#' fun<-getDefaultIsThermallingFunction(170)
#' fun(1:160)
#' fun(1:190)
getDefaultIsThermallingFunction <- function(totalAngle = 360) {
  stopifnot(is.numeric(totalAngle))
  stopifnot(length(totalAngle)==1)
  function(headings) {
    diffHeading = diff(headings)
    diffHeading <- ((diffHeading + 180) %% 360) - 180
    max(abs(cumsum(diffHeading))) >= totalAngle
  }
}

#' getSamplingIsRegularFunction
#'
#' @param samplingIntervalSeconds the interval that is considered regular
#'
#' @return a function is returned that based on a series of timestamps decides if the segment is regular
#' @export
#'
#' @examples
#' fun<-getSamplingIsRegularFunction(10)
#' fun(Sys.time()+1:5)
#' fun(Sys.time()+c(0,10,20,30))
#' fun(Sys.time()+c(0,10,20,31))
getSamplingIsRegularFunction <- function(samplingIntervalSeconds) {
  stopifnot(is.numeric(samplingIntervalSeconds))
  stopifnot(length(samplingIntervalSeconds)==1)
  function(timestamps) {
    retval = all(diff(as.numeric(timestamps)) == samplingIntervalSeconds)
    if (is.na(retval))
      retval = F
    retval
  }
}
