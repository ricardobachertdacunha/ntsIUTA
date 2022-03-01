

### patFragmentSettingsDefault -----

#' @title patFragmentSettingsDefault-class
#'
#' @return An \linkS4class{ntsSettings} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object
#' based on the functionality from \pkg{patRoon} through \pkg{mzR}.
#'
#' @export
#'
patFragmentSettingsDefault <- function() {

  return(
    new(
      "ntsSettings",
      algorithm = "mzr",
      settings = list(
        maxMSRtWindow = 10,
        precursorMzWindow = 4,
        topMost = NULL,
        avgFeatParams = list(
          clusterMzWindow = 0.005,
          topMost = 500,
          minIntensityPre = 20,
          minIntensityPost = 20,
          avgFun = mean,
          method = "distance",
          pruneMissingPrecursorMS = FALSE,
          retainPrecursorMSMS = TRUE
        ),
        avgFGroupParams = list(
          clusterMzWindow = 0.008,
          topMost = 500,
          minIntensityPre = 20,
          minIntensityPost = 20,
          avgFun = mean,
          method = "distance",
          pruneMissingPrecursorMS = FALSE,
          retainPrecursorMSMS = TRUE
        )
      )
    )
  )
}




### fragmentSettingsDefault -----

#' @title fragmentSettingsDefault-class
#'
#' @return An \linkS4class{ntsSettings} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object.
#'
#' @export
#'
fragmentSettingsDefault <- function() {

  return(
    new(
      "ntsSettings",
      algorithm = "ntsiuta",
      settings = list(
        isolationTimeWindow = 5,
        isolationMassWindow = 1.3,
        clusteringMethod = "distance",
        clusteringUnit = "ppm",
        clusteringWindow = 25,
        minIntensityPre = 100,
        minIntensityPost = 100,
        asPatRoon = TRUE,
        mergeCEs = TRUE,
        mergeBy = NULL
      )
    )
  )
}




#' @title generateMS2
#'
#' @param object An \linkS4class{ntsData} or a \linkS4class{featureGroups} object to
#' extract MS2 data according to percursor ions.
#' @param algorithm A character vector with the algorithm to generate MS2 data.
#' @param settings A list of settings to use for extraction of MS2 data.
#'
#' @return A \linkS4class{MSPeakLists} when \code{object} is a \linkS4class{featureGroups}
#' or an \linkS4class{ntsData} with MS2 data added to the features data.table
#' when \code{object} is an \linkS4class{ntsData} object.
#'
#' @export
#'
#' @importFrom patRoon getDefAvgPListParams
#' @importMethodsFrom patRoon generateMSPeakLists
#'
#'
generateMS2 <- function(object,
                        algorithm = NULL,
                        settings = NULL) {

  pat <- object
  if (checkmate::testClass(object, "ntsData")) {
    pat <- pat@pat
    if (is.null(algorithm)) algorithm <- fragmentsParameters(object)@algorithm
    if (is.null(settings)) settings <- fragmentsParameters(object)@settings
  }

  if (algorithm != "ntsiuta") {
    ag <- list(fGroups = pat, algorithm = algorithm)

    MS2 <- suppressWarnings(
      do.call(generateMSPeakLists,
        c(ag, settings)
      )
    )
  } else {
    MS2 <- suppressWarnings(
      do.call(extractMSn,
        c(list(object = object), settings)
      )
    )
  }
  #TODO add option to add MS2 to the ntsData object

  return(MS2)
}
