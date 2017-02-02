#' Example stork data.
#'
#' A dataset containing location data of 6 juvenile storks (Ciconia ciconia) on the 18th of august when migration just started. On several occasion the birds use thermals.
#'
#' @format A MoveStack consisting of 22333 locations
#'
#' @source \url{http://www.movebank.org/}
#' @examples
#' data("storks")
"storks"


if (F) {
  require(move)
  l <- movebankLogin()
  s <-
    grep("2015",
         value = T,
         grep('Life Tra', getMovebankStudies(l), value = T))
  a <- getMovebankAnimals(s, l)
  
  ids <-
    unique(a$individual_id[sub(').*', '', trim(sub('.*eobs', '', a$animalName_deployment))) %in% c("3666", "3929", "3930", "3916", "3907", "3909")])
  trks <-
    getMovebankData(s,
                    ids,
                    login = l,
                    removeDuplicatedTimestamps = T)
  trksSel <- trks[as.Date(timestamps(trks)) == "2014-08-18", ]
  for (i in c(
    "eobs_start_timestamp",
    "eobs_activity",
    "eobs_accelerations_raw",
    "sensor_type_id",
    "eobs_acceleration_axes",
    "eobs_acceleration_sampling_frequency_per_axis",
    "eobs_battery_voltage",
    "eobs_fix_battery_voltage",
    "eobs_activity_samples",
    "eobs_acceleration_axis",
    "eobs_key_bin_checksum",
    "timestamp",
    "location_long",
    "location_lat",
    "deployment_id"
  )) {
    trksSel@data[, i] <- NULL
    trksSel@dataUnUsedRecords[, i] <- NULL
  }
  
  validObject(trksSel)
  storks <- trksSel
  require(devtools)
  use_data(storks)
}