#' @importFrom parallel mclapply
#' @importFrom move angle
#' @importFrom move timestamps
#' @importFrom move speed
#' @importFrom move split
#' @importFrom move moveStack
#' @importClassesFrom move Move
#' @importClassesFrom move MoveStack
NULL
#' Title
#'
#' @param data Move object, MoveStack or data.frame containing wind speeds
#' @param timestamps timestamps of the speed observations
#' @param windowSize a numeric vector of length 1 or 2, if length 1 it is the size of the focal window data will be assigned to the central location. If length 2 the window size is \code{sum(windowSize)+1)} and the first element is the number of location before the focal locations, the second is the number of locations after the focal location.
#' @param phi todo
#' @param isFocalPoint todo
#' @param isSamplingRegular either a function that determines based on a vector of timestamps if the sampling interval is regular or a numeric value that corresponds to the time interval between observations in the dataset that is regular
#' @param hasVariationInHeadingFunction todo
#' @param columnNamesWind todo
#' @param ... other possible arguments currently nothing else is implented
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
#' windEst<-spTransform(windEst, center=TRUE)
#' plot(windEst)
#' # only plot few arrows of estimates
#' s<-windEst$estimationSuccessful & format(timestamps(windEst), "%S")=='00'
#' # enlarge arrows 30 times
#' arrows(coordinates(windEst)[s,1],coordinates(windEst)[s,2], 
#'    coordinates(windEst)[s,1]+ windEst$windX[s]*30, 
#'    coordinates(windEst)[s,2]+windEst$windY[s]*30)
setGeneric("getWindEstimates", function(data,timestamps,...){
  standardGeneric("getWindEstimates")
  })

#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,MoveStack,missing-method
setMethod("getWindEstimates",signature = signature(data="MoveStack", timestamps="missing") ,
          function(data,timestamps,...){
x<-lapply(move::split(data), getWindEstimates,...)
            
            return(moveStack(x))
            }
)
#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,Move,missing-method
setMethod("getWindEstimates",signature = signature(data="Move", timestamps="missing") ,
          function(data,timestamps,...){
            
                      dataDf<-data.frame(vx=speed(data)*sin(angle(data)/180*pi),vy=speed(data)*cos(angle(data)/180*pi))
                             res<-callGeneric(data=rbind(dataDf,NA),timestamps=move::timestamps(data), ...)
          slot(data,"data")<-cbind(slot(data,'data'), res)
            return(data)
            }
)
                   
#' @rdname getWindEstimates-methods
#' @aliases getWindEstimates,data.frame,POSIXct-method
setMethod("getWindEstimates",signature = signature(data="data.frame", timestamps='POSIXct') ,
 function(data,timestamps,
                                  windowSize=29,
                                  phi = 0.5,
                                  isFocalPoint = function(i, ts){TRUE},
                                  isSamplingRegular = 1,
                                  hasVariationInHeadingFunction = getDefaultIsThermallingFunction(360),
                                  columnNamesWind = c(
                                    "estimationSuccessful", "residualVarAirspeed", "windX", "windY", "windVarX", "windVarY", "windCovarXY", "windVarMax", "isThermalling"
                                  ),...)
{
   # first check some conditions
   stopifnot(nrow(data)==length(timestamps))
   stopifnot(any(length(windowSize)==1:2))
   stopifnot(is.numeric(windowSize))
   stopifnot(all(windowSize>=0))
   stopifnot(length(isSamplingRegular)==1)
   if(is.numeric(isSamplingRegular)){
     isSamplingRegularFun<-getSamplingIsRegularFunction(isSamplingRegular)}else{
       isSamplingRegularFun<-isSamplingRegular
     }
   stopifnot(class(isSamplingRegularFun)=="function")
   res<-data.frame(matrix(NA, ncol=9, nrow=nrow(data)))
   names(res)<-columnNamesWind
   if(length(windowSize)==1){
windowSize<-c(floor((windowSize-1)/2),ceiling((windowSize-1)/2))
}
  selRange <- -windowSize[1]:windowSize[2]
  samples <- 1:nrow(data)
  samples <-
    samples[samples > windowSize[1] &
              samples <= (nrow(data) - windowSize[2]) &
              isFocalPoint(samples, timestamps)]
  if(!(length(samples)>0)) {stop("No sample locations selected, either throught the focal function or large windowSize")}

  stopifnot(length(phi)==1)
  stopifnot(is.numeric(phi))
  
  oo <- lapply(samples,'+', selRange)
    
  times <- lapply(oo ,function(x,i) {x[i]}, x = timestamps)
  reg <- unlist(lapply(times, isSamplingRegularFun))
  
  groundSpeedsL <-
    lapply(oo[reg],function(i,x) {
      x[i,]
    }, as.matrix(data, rownames.force = F))
  
  naGrnd <- !unlist(lapply(lapply(groundSpeedsL, is.na),any))
  samples <- samples[reg][naGrnd]
  
  groundSpeedsL <- groundSpeedsL[naGrnd]
  if(length(groundSpeedsL)==0)
    return(res)
  windEst <- mclapply(groundSpeedsL,getWindEstimate, phi = phi)

  res[samples,columnNamesWind[1]] <- (s <-
                                         !unlist(lapply(windEst, is.null)))
  samples <- samples[s]
  windEst <- windEst[s]
  res[samples,columnNamesWind[2:4]] <-
    matrix(unlist(lapply(
      windEst, '[',c("residualVarAirSpeed","windEst")
    )), ncol = 3, byrow = T)
  covarList <- lapply(windEst,'[[','covar')
  nonNullCovar <- !unlist(lapply(covarList, is.null))
  res[samples[nonNullCovar], columnNamesWind[5:7]] <-
    matrix(unlist(covarList[nonNullCovar]), ncol = 4, byrow = T)[,c(1,4,2)]
  nonNaNCovar <- !unlist(lapply(lapply(covarList[nonNullCovar], is.nan), any))
  res[samples[nonNullCovar][nonNaNCovar], columnNamesWind[8]] <-
    unlist(lapply(covarList[nonNullCovar][nonNaNCovar], eigen,only.values = T))[c(T,F)]
  res[samples, columnNamesWind[9]] <-
    unlist(mclapply(
      lapply(mapply(
        '-',lapply(groundSpeedsL[s],t), lapply(windEst,'[[','windEst'), SIMPLIFY =
          F
      ),function(x) {
        atan2(x[1,],x[2,]) / pi * 180
      }), hasVariationInHeadingFunction
    ))
  
  res[!(res[,columnNamesWind[9]] |
           is.na(res[,columnNamesWind[9]])), columnNamesWind[2:8]] <- NA
  return(res)
})
