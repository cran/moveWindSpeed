#' Function to find good points for estimation of phi
#'
#' The function tries to find non overlapping windows for phi optimization.
#' @param data  An move object.
#' @param maxPointsToUseInEstimate The number of desired windows.
#' @param phiInitialEstimate  The initial value used for the autocorrelation when calculating the wind speed for finding suitable windows. 
#' @param windowSize  An odd number providing the window size
#' @param ...  passed on to getWindEstimates 
#'
#' @return a logical vector with the focal locations
#' @export
#'
#' @examples
#' data(storks)
#' which(findGoodPoints( storks[[2]],
#' windowSize = 29,  isSamplingRegular = 1,
#' isThermallingFunction = getDefaultIsThermallingFunction(360, 4),  maxPointsToUseInEstimate = 10,
#' phiInitialEstimate = 0  ))
findGoodPoints <- function(data,
                           maxPointsToUseInEstimate,
                           phiInitialEstimate,
                           windowSize, 
                           ...) {
  # rearrange points so that samples are spread out over the whole sample
  shuffledIndices = sample(seq(1, nrow(data), 1)) # windowSize
  # start with number of points requested
  numPoints = min(c(maxPointsToUseInEstimate*2, nrow(data)))
#  isGoodPoint = rep(F, nrow(data))
  w<-getWindowSizeLR(windowSize)
  goodPoints<-c()
  while (T) {
    # only use the first numPoints points from the shuffled list
    t<-min(c(numPoints, length(shuffledIndices)))
    isFocalPoint =  head(shuffledIndices,t)
    shuffledIndices<-tail(shuffledIndices,-t)
    windEst = getWindEstimates(
      data,
      phi = phiInitialEstimate,
      isFocalPoint =
        isFocalPoint,
      windowSize = windowSize,
      ...
    )
    # find points where the bird is thermalling
    pointsSelected<-isFocalPoint[( !is.na(windEst$estimationSuccessful) &
      windEst$estimationSuccessful)[isFocalPoint]]
    sOmit<-c()
    while(length(pointsSelected)){
      s<-pointsSelected[1]
      pointsSelected<-pointsSelected[-1]
      omit<-(s-sum(w)):(s+sum(w))
      sOmit<-c(omit, sOmit)
      pointsSelected<-pointsSelected[!(pointsSelected%in%omit)]
      goodPoints<-c(goodPoints,s)
    }
    # performance seems slightly better if this is taken out of loop
    shuffledIndices<- shuffledIndices[!(shuffledIndices%in%sOmit)]
    # stop if enough 'good' points have been found
    if (length(goodPoints) >= maxPointsToUseInEstimate)
      break
    # stop if there are no more input data
    if (length(shuffledIndices)<1)
      break
    numPoints = numPoints * 2
  }
  if (length(goodPoints) > maxPointsToUseInEstimate) {
    goodPoints = goodPoints[1:maxPointsToUseInEstimate]
  }
  
  if (length(goodPoints) < maxPointsToUseInEstimate)
  {
    warning('Desired number of locations could not be found')
  }
  isGoodPoint = rep(F, nrow(data))
  isGoodPoint[goodPoints]<-T
  return(isGoodPoint)
}
#' estimatePhi
#'
#' An function to estimate phi (the autocorrelation of speed) from data. This is done using iterative calls to the wind speed optimization on a selection of segments.
#' @param data An move object or stack.
#' @param isThermallingFunction The thermalling function to use.
#' @param maxPointsToUseInEstimate  Maximal number of desired windows for phi estimation
#' @param phiInitialEstimate  Initial phi estimate
#' @param isGoodPoint  The points to use for phi estimation as logical or numeric, if NULL then findGoodPoints is used.
#' @param returnPointsUsedInEstimate an logical value, if the segments used for phi estimation should also be returned.
#' @param windowSize  An window size, odd number or the start and end of the window relative to the focal point
#' @param ... extra arguments for getWindSpeedEstimates
#'
#' @return a list with phi and the log likelihood and the number of locations used
#' @export
#'
#' @examples
#' data(storks)
#' estimatePhi(
#'   storks[[2]],
#'   windowSize = 19,
#'   isSamplingRegular = 1,
#'   isThermallingFunction = getDefaultIsThermallingFunction(360, 4),
#'   maxPointsToUseInEstimate = 10
#' )
estimatePhi <-
  function(data,
           isThermallingFunction = getDefaultIsThermallingFunction(360, 4),
           maxPointsToUseInEstimate = 20,
           phiInitialEstimate = 0,
           isGoodPoint = NULL,
           returnPointsUsedInEstimate = F,
           windowSize = 29,
           ...) {
    if (is.null(isGoodPoint)) {
      isGoodPoint = findGoodPoints(
        data = data,
        isThermallingFunction = isThermallingFunction,
        windowSize = windowSize,
        phiInitialEstimate = phiInitialEstimate,
        maxPointsToUseInEstimate = maxPointsToUseInEstimate,
        ...
      )
    }
    if(is.logical(isGoodPoint)){
      stopifnot(length(isGoodPoint)==nrow(data))
    if(!any(isGoodPoint))
      stop("no suitable point found or provided for phi estimation")}
    segments<-getWindEstimates(
      data,
      isFocalPoint = isGoodPoint,
      isThermallingFunction = isThermallingFunction,
      windowSize = windowSize,
      returnSegmentList=T,
      ...
    )
    if(is.logical(isGoodPoint)){
    goodPts<-which(isGoodPoint)}
    f_lik = function(phi) {
      windEst <-
        getWindEstimates(
          segments,
          phi = phi,
          isThermallingFunction = isThermallingFunction,
          ...
        )
      sigma = windEst$residualVarAirspeed
      # remove good points if wind estimation fails for some phi:
      if (any(is.na(sigma))) {
        warning(
          paste(
            "Location(s) omited because it fails for phi=",
            phi,
            "first location:",
            goodPts[
                    is.na(sigma)][1]
          )
        )
        segments<<-segments[!is.na(sigma)]
        goodPts<<-goodPts[!is.na(sigma)]
      }
      sigma = sigma[!is.na(sigma)]
      # number of groups = number of good points
      windEstimLogLik(sigma, phi)
    }
    # optimize
    result = optim(0.5,
                   f_lik,
                   method = "Brent",
                   lower = 0,
                   upper = 1)
    if(is.numeric(isGoodPoint))
    {  nIsGoodPoint<-length(isGoodPoint)
    }else{
      nIsGoodPoint<-sum(isGoodPoint)
      }
    retval = list(
      phi = result$par,
      logLik = result$value,
      n = nIsGoodPoint
    )
    if (returnPointsUsedInEstimate){
      retval[["isPointUsedInEstimate"]] = isGoodPoint
      retval[["segments"]]<-segments
    }
    return(retval)
  }
