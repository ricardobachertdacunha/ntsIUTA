

### patFragmentSettingsDefault -----

#' @title patFragmentSettingsDefault
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

#' @title fragmentSettingsDefault
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
        clusteringWindow = 5,
        minIntensityPre = 10,
        minIntensityPost = 10,
        asPatRoon = TRUE,
        mergeVoltages = TRUE,
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
                        algorithm = NA_character_,
                        settings = NULL) {

  pat <- object
  if (checkmate::testClass(object, "ntsData")) {
    pat <- pat@pat
    if (is.na(algorithm)) algorithm <- fragmentsParameters(object)@algorithm
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

  return(MS2)
}




#' @title loadMS2
#'
#' @description Function to load MS2 data from features in the feature table
#' of a \linkS4class{ntsData} object. The MS2 data is added to the \code{ms2}
#' slot of the \linkS4class{ntsData} object. If the \code{algorithm} and
#' \code{settings} arguments are not given, the fucntion uses the parameters in the
#' \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character vector with the algorithm to generate MS2 data.
#' @param settings A list of settings to use for extraction of MS2 data.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} with MS2 data added to the \code{ms2} slot.
#'
#' @export
#'
loadMS2 <- function(object,
                    algorithm = NA_character_,
                    settings = NULL,
                    save = FALSE) {

  assertClass(object, "ntsData")

  if (is.na(algorithm)) {
    algorithm <- fragmentsParameters(object)@algorithm
  }

  if (is.na(algorithm)) {
    warning("Fragments algorihtm not defined!")
    return(object)
  }

  if (is.null(settings)) {
    settings <- fragmentsParameters(object)@settings
  }

  MS2 <- generateMS2(
    object,
    algorithm = algorithm,
    settings = settings
  )

  object@ms2 <- MS2

  if (is.logical(save)) if (save) saveObject(object = object)

  if (is.character(save)) saveObject(object = object, filename = save)

  return(object)
}
