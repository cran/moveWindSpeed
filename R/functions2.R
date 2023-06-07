#' Generate arguments for window size around focal point
#'
#' A function to translate an window size argument to a standardized argument.
#' @param windowSize a pair of positive integers determining the window size left and right of a focal point or an odd number determining the size of a symmetrical window
#'
#' @return windowSize a pair of positive integers determining the window size left and right of a focal point
#' @export
getWindowSizeLR = function(windowSize) {
  if(!is.numeric(windowSize) || !any(length(windowSize) == 1:2) || any(windowSize%%1!=0) || any(windowSize < 0)) stop("windowSize must be an integer or a pair of positive integers")
  if (length(windowSize) == 1) {
    if(!(windowSize %% 2 == 1)) stop("if a single value is given for windowSize, it must be an odd number")
    windowSize <- c(floor((windowSize - 1) / 2), ceiling((windowSize - 1) / 2))
  }
  windowSize
}

#' A function to generate functions used to check if a segment is regular
#'
#' @param isSamplingRegular a function which decides if a sequence of timestamps is regular or the interval which is considered regular
#'
#' @return a function which decides if a sequence of timestamps is regular
#' @export
#'
#' @examples
#' fun<-getIsSamplingRegularFunction(10)
#' fun(Sys.time()+1:5)
#' fun(Sys.time()+c(0,10,20,30))
#' fun(Sys.time()+c(0,10,20,31))
getIsSamplingRegularFunction <- function(isSamplingRegular) {
  stopifnot(length(isSamplingRegular) == 1)
  if (is.numeric(isSamplingRegular)) {
    stopifnot(length(isSamplingRegular)==1)
    isSamplingRegularFun <- function(timestamps) {
      retval = all(diff(as.numeric(timestamps)) == isSamplingRegular)
      if (is.na(retval))
        retval = F
      retval
    }
  } else{
    isSamplingRegularFun <- isSamplingRegular
  }
  if(!inherits(isSamplingRegularFun, "function")) stop("isSamplingRegular must be a number or a boolean function on a vector of timestamps")
  isSamplingRegularFun
}

#' A function to generate isFocalPoint functions
#'
#' @param isFocalPoint a function, a boolean array from which such a function can be built, or a list of indices
#'
#' @return a function which decides if wind estimation is performed for a point in the input data
#' @export
getIsFocalPointFunction <- function(isFocalPoint) {
  if (is.logical(isFocalPoint)) {
#    if(length(isFocalPoint) != nrow(data)) stop("If a boolean vector is passed for isFocalPoint, it must have the same size as the input data")
    isFocalPoint <- which(isFocalPoint)
  }
  if (is.numeric(isFocalPoint)) {
    u <- isFocalPoint
    isFocalPoint <- function(i, ts) {
      i %in% u
    }
  }
  if(!inherits(isFocalPoint, "function")) stop("isFocalPoint must be a vector of numbers, a boolean vector or a function")
  return(isFocalPoint)
}

