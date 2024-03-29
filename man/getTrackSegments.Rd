% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getTrackSegments.R
\name{getTrackSegments}
\alias{getTrackSegments}
\title{An helper function to extract trajectory segments for wind estimation from a track}
\usage{
getTrackSegments(
  data,
  timestamps,
  windowSize = 29,
  isFocalPoint = function(i, ts) {
     TRUE
 },
  isSamplingRegular = 1,
  focalSampleBefore = 0
)
}
\arguments{
\item{data}{A two column dataframe.}

\item{timestamps}{A series of POSIXct timestamps as long as the data.}

\item{windowSize}{The window size (odd number) or two numbers giving the start and end of a window around a focal point.}

\item{isFocalPoint}{an function taking location numbers and timestamps that is used to see if a location should be considered as an focal point. It can for example be used to speed up calculations by only considering every second location. An numeric value can also be provided then only these locations are considered}

\item{isSamplingRegular}{Either an numeric or a function that is used to decide if a series of timestamps is regular. If numeric than it should correspond to the interval in seconds.}

\item{focalSampleBefore}{An argument to be used if data is not the start of the location count.}
}
\value{
A list of ground speeds
}
\description{
An helper function to extract trajectory segments for wind estimation from a track
}
\examples{
length(getTrackSegments(data.frame(1:40,1:40), Sys.time()+1:40))
length(getTrackSegments(data.frame(1:40,1:40), Sys.time()+c(1:25,36:50), windowSize=11))
str(getTrackSegments(data.frame(1:40,1:40), Sys.time()+1:40, windowSize=39))
}
