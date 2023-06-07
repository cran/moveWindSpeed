#' Generate wind estimates for a trajectories or data frame with wind speeds
#'
#' @param data Move object, MoveStack or data.frame containing wind speeds
#' @param timestamps timestamps of the speed observations
#' @param windowSize a numeric vector of length 1 or 2, if length 1 it is the size of the focal window data will be assigned to the central location. If length 2 the window size is \code{sum(windowSize)+1)} and the first element is the number of location before the focal locations, the second is the number of locations after the focal location.
#' @param phi the auto correlation of air speed.
#' @param isFocalPoint an function that based on location number and timestamps returns a logical vector if location should be included. Or a numeric/logical vector indicating the location numbers.
#' @param isSamplingRegular either a function that determines based on a vector of timestamps if the sampling interval is regular or a numeric value that corresponds to the time interval between observations in the dataset that is regular
#' @param isThermallingFunction An function that based on a series of headings and speeds (wind corrected) decides if an segment should be considered thermalling.
#' @param columnNamesWind The column names used for storing the data in the returned objected after it has been calculated.
#' @param groundSpeedXY an character of length 2 containing column names from the move object that need to be used as the x and y component of the ground speed vector
#' @param focalSampleBefore The number of locations that occurred before the move object fed in the getWindEstimates function, used in case stacks are provided for example. This is most cases not useful for users.
#' @param returnSegmentList a logical value indicating if the list of segments to estimate wind over should be returned instead of the estimates
#' @param referenceGroundSpeed a number indicating which of the grounds speed vectors to take as a reference for air speed, by default the 0th/middle location of the window if that is specified by one number.
#' @param ... other possible arguments currently nothing else is implemented
#'
#' @return a Move object, dataframe or a MoveStack depending on input
#' @export
#' @docType methods
#' @rdname getWindEstimates-methods
#'
#' @examples
#' data("storks")
#' # run example for reduced dataset
#' windEst<-getWindEstimates(storks[format(timestamps(storks),"%H")=="12",][[2:3]])
#' # Use evolution status 2 to avoid using rgdal (set using sp)
#' set_evolution_status(2L)
#' windEst<-spTransform(windEst, center=TRUE)
#' plot(windEst)
#' # only plot few arrows of estimates
#' s<-windEst$estimationSuccessful & format(timestamps(windEst), "%S")=='00'
#' # enlarge arrows 30 times
#' arrows(coordinates(windEst)[s,1],coordinates(windEst)[s,2],
#'    coordinates(windEst)[s,1]+ windEst$windX[s]*30,
#'    coordinates(windEst)[s,2]+windEst$windY[s]*30)
setGeneric("getWindEstimates", function(data, timestamps, ...) {
  standardGeneric("getWindEstimates")
})

#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,MoveStack,missing-method
setMethod("getWindEstimates", signature = signature(data = "MoveStack", timestamps =
                                                      "missing") ,
          function(data, timestamps, ...) {
            l <- list(...)
            if ('isFocalPoint' %in% names(l))
              if (is.logical(l[['isFocalPoint']])) {
                stopifnot(length(l[['isFocalPoint']]) == sum(n.locs(data)))
                l[['isFocalPoint']] <- which(l[['isFocalPoint']])
              }
            slist <-         move::split(data)
            x <-
              mapply(
                function(...) {
                  tryCatch(
                    getWindEstimates(...),
                    error = function(e) {
                      if (e$message == 'No sample locations selected, either through the focal function or large windowSize') {
                        NULL
                      } else{
                        stop(e)
                      }
                    }
                  )
                },
                slist
                ,
                focalSampleBefore = head(c(0, cumsum(n.locs(
                  data
                ))),-1),
                MoreArgs = l,
                SIMPLIFY = F
              )
            if (all(s <- unlist(lapply(x, is.null))))
            {
              stop('No sample locations selected, either through the focal function or large windowSize')
            }
            
            if (unique(unlist(lapply(x[!s], class))) != "Move")
            {
              return(unlist(recursive = F, x))
            }

            res <- moveStack(ifelse(s, slist, x))
            return(res)
          })
#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,Move,missing-method
setMethod("getWindEstimates", signature = signature(data = "Move", timestamps =
                                                      "missing") ,
          function(data, timestamps, groundSpeedXY = NULL, ...) {
            if (!is.null(groundSpeedXY))
            {
              stopifnot(is.character(groundSpeedXY))
              stopifnot(length(groundSpeedXY) == 2)
              stopifnot(all(groundSpeedXY %in% names(data)))
              dataDf <- data.frame(data[, groundSpeedXY])[, groundSpeedXY]
            } else{
              spd <- speed(data)
              ang <- angle(data)
              dataDf <-
                data.frame(vx = spd * sin(ang / 180 * pi),
                           vy = spd * cos(ang / 180 * pi))
              dataDf <- rbind(dataDf, NA)
            }
            res <-
              callGeneric(data = dataDf,
                          timestamps = move::timestamps(data),
                          ...)
            if(!is.data.frame(res))
            {
              return(res)
            }
            slot(data, "data") <- cbind(slot(data, 'data'), res)
            return(data)
          })

#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,data.frame,POSIXct-method
setMethod("getWindEstimates", signature = signature(data = "data.frame", timestamps =
                                                      'POSIXct') ,
          function(data,
                   timestamps,
                   windowSize = 29,
                   isFocalPoint = function(i, ts) {
                     TRUE
                   },
                   isSamplingRegular = 1,
                   focalSampleBefore = 0,
                   returnSegmentList=F,                   referenceGroundSpeed=NULL,
                   ...)
          {
            l <- nrow(data)
            segList <-
              getTrackSegments(
                data,
                windowSize = windowSize,
                isFocalPoint = isFocalPoint,
                isSamplingRegular = isSamplingRegular, 
                focalSampleBefore=focalSampleBefore, timestamps=timestamps
              )
            data <- segList
            if(returnSegmentList){
              return(data)
            }
            w <- getWindowSizeLR(windowSize)
            if (is.null(referenceGroundSpeed)) {
              referenceGroundSpeed <- which(((-w[1]):w[2]) == 0)
            }
            res <-
              callGeneric(data = data, referenceGroundSpeed = referenceGroundSpeed, ...)
            resFull<-res[F,][1:l,]
            resFull[as.numeric(names(segList)),]<-res
            return(resFull)
          })
#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,list,ANY-method
setMethod("getWindEstimates", signature = signature(data = "list", timestamps =
                                                      'ANY') ,
          function(data,
                   timestamps,
                   phi = 0,
                   isThermallingFunction = getDefaultIsThermallingFunction(360, 4),
                   columnNamesWind = c(
                     "estimationSuccessful",
                     "residualVarAirspeed",
                     "windX",
                     "windY",
                     "windVarX",
                     "windVarY",
                     "windCovarXY",
                     "windVarMax",
                     'airX',
                     'airY'
                   ),
                   referenceGroundSpeed=NULL,
                   ...)
          {
            n<-length(data)
            windEst <- mclapply(data, getWindEstimate, phi = phi)
            res <- data.frame(matrix(NA, ncol = length(columnNamesWind), nrow = n))
            names(res) <- columnNamesWind
            estimationSuccessful <- !unlist(lapply(windEst, is.null))
            airVectors <-
              mapply('-',
                     lapply(data[estimationSuccessful], t),
                     lapply(windEst[estimationSuccessful], '[[', 'windEst'),
                     SIMPLIFY =
                       F)  
            airHeadings <- lapply(airVectors
                                  , function(x) {
                                    atan2(x[1, ], x[2, ]) / pi * 180
                                  })
            airSpeeds <- lapply(airVectors
                                , function(x) {
                                  sqrt(colSums(x ^ 2))
                                })
            isThermalling = rep(F, n)
            isThermalling[estimationSuccessful] <-
              unlist(mcmapply(
                headings = airHeadings,
                speeds = airSpeeds
                ,
                isThermallingFunction
              ))

            estimationSuccessful = estimationSuccessful & isThermalling
            res[, columnNamesWind[1]] <- estimationSuccessful
            if(all(!estimationSuccessful))
            {  return(res)}
            res[estimationSuccessful, columnNamesWind[2:4]] <-
              matrix(unlist(lapply(
                windEst[estimationSuccessful], '[', c("residualVarAirSpeed", "windEst")
              )), ncol = 3, byrow = T)
            covarList <- lapply(windEst[estimationSuccessful], '[[', 'covar')
            nonNullCovar = rep(F, n)
            nonNullCovar[estimationSuccessful] <- !unlist(lapply(covarList, is.null))
            res[estimationSuccessful & nonNullCovar, columnNamesWind[5:7]] <-
              matrix(unlist(covarList[nonNullCovar[estimationSuccessful]]), ncol = 4, byrow = T)[, c(1, 4, 2)]
            nonNaNCovar = rep(F, n)
            nonNaNCovar[estimationSuccessful & nonNullCovar] <-
              !unlist(lapply(lapply(covarList[nonNullCovar], is.nan), any))
            res[estimationSuccessful & nonNullCovar & nonNaNCovar, columnNamesWind[8]] <-
              unlist(lapply(covarList[nonNullCovar[estimationSuccessful] & nonNaNCovar[estimationSuccessful]], eigen, only.values = T))[c(T, F)]
            if(!is.null(referenceGroundSpeed))
            {
              grndS<-do.call("rbind",lapply(data,'[', referenceGroundSpeed, T))
              res[, columnNamesWind[9:10]]<-
                grndS-res[,c("windX","windY")]

            }
            rownames(res)<-names(data)
            return(res)
})
