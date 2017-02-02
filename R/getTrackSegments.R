#' An helper function to extract trajectory segments for wind estimation from a track
#'
#' @param data A two column dataframe.
#' @param timestamps A series of POSIXct timestamps as long as the data.
#' @param windowSize The window size (odd number) or two numbers giving the start and end of a window around a focal point. 
#' @param isFocalPoint an function taking location numbers and timestamps that is used to see if a location should be considered as an focal point. It can for example be used to speed up calculations by only considering every second location. An numeric value can also be provided then only these locations are considered
#' @param isSamplingRegular Either an numeric or a function that is used to decide if a series of timestamps is regular. If numeric than it should correspond to the interval in seconds.
#' @param focalSampleBefore An argument to be used if data is not the start of the location count.
#'
#' @return A list of ground speeds
#' @export
#'
#' @examples
#' length(getTrackSegments(data.frame(1:40,1:40), Sys.time()+1:40))
#' length(getTrackSegments(data.frame(1:40,1:40), Sys.time()+c(1:25,36:50), windowSize=11))
#' str(getTrackSegments(data.frame(1:40,1:40), Sys.time()+1:40, windowSize=39))

getTrackSegments <- function(data,
                             timestamps,
                             windowSize = 29,
                             isFocalPoint = function(i, ts) {
                               TRUE
                             },
                             isSamplingRegular = 1,
                             focalSampleBefore = 0) {
  windowSize = getWindowSizeLR(windowSize)
  isFocalPoint = getIsFocalPointFunction(isFocalPoint)
  isSamplingRegular = getIsSamplingRegularFunction(isSamplingRegular)
  stopifnot(ncol(data) == 2)
  stopifnot(nrow(data)== length(timestamps))
  focalPoints <- 1:nrow(data)
  focalPoints <-
    focalPoints[focalPoints > windowSize[1] &
                  focalPoints <= (nrow(data) - windowSize[2]) &
                  isFocalPoint(focalSampleBefore + focalPoints, timestamps)]
  if (!(length(focalPoints) > 0)) {
    stop("No sample locations selected, either through the focal function or large windowSize")
  }
  selRange <- -windowSize[1]:windowSize[2]
  segments <- lapply(focalPoints, '+', selRange)
  times <- lapply(segments , function(x, i) {
    x[i]
  }, x = timestamps)
  isGoodSegment <- unlist(lapply(times, isSamplingRegular))
  if (!is.null(data$animalIds)) {
    animals <- lapply(segments , function(x, i) {
      x[i]
    }, x = data$animalIds)
    isGoodSegment = isGoodSegment &
      unlist(lapply(animals, function(l)
        all(l == l[1])))
  }
  if (!any(isGoodSegment))
    stop("No sample locations selected due to irregular sampling")
  focalPoints = focalPoints[isGoodSegment]
  segments = segments[isGoodSegment]
  
  groundSpeedsL <-
    lapply(segments, function(i, x) {
      x[i,]
    }, as.matrix(data[, 1:2], rownames.force = F))
  naGrnd <- !unlist(lapply(lapply(groundSpeedsL, is.na), any))
  focalPoints <- focalPoints[naGrnd]
  segments <- segments[naGrnd]
  groundSpeedsL <- groundSpeedsL[naGrnd]
  if (length(groundSpeedsL) == 0)
    stop("No sample locations selected due to missing ground speed values")
  names(groundSpeedsL)<-focalPoints
  groundSpeedsL<-mapply(`rownames<-`, groundSpeedsL, segments, SIMPLIFY = F)
  return(groundSpeedsL)
}
