

### fillingSettingsDefaultXCMS -----

#' @title fillingSettingsDefaultXCMS
#'
#' @return An \linkS4class{ntsSettings} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object.
#'
#' @export
#'
fillingSettingsDefaultXCMS <- function() {

  return(
    new(
      "ntsSettings",
      algorithm = "xcms",
      settings = list(
        xcms::ChromPeakAreaParam()
      )
    )
  )
}




#' @title peakFilling
#'
#' @description Recursive integration for filling missing peaks within each feature.
#' The function \code{\link[xcms]{fillChromPeaks}} from \pkg{xcms} package can be used,
#' giving the respectives parameters with the argument \code{settings}.
#'
#' @param object An \linkS4class{ntsData} object containing features.
#' @param algorithm A character vector defining the algorithm to be used for filling peaks.
#' @param settings A list with parameters according to the defined algorithm.
#' When the algorithm is set to \emph{xcms3}, the settings are the S4 class objects
#' \linkS4class{FillChromPeaksParam} or \linkS4class{ChromPeakAreaParam}.
#' See \code{?\link[xcms]{fillChromPeaks}} for more information.
#'
#' @return An \linkS4class{XCMSnExp} object including filled missing peaks.
#'
#' @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#'
#' @importMethodsFrom patRoon getXCMSnExp analysisInfo
#' @importFrom patRoon importFeatureGroupsXCMS3
#' @importMethodsFrom xcms fillChromPeaks
#' @importClassesFrom xcms ChromPeakAreaParam FillChromPeaksParam XCMSnExp
#' @importClassesFrom patRoon featuresSIRIUS featureGroupsSIRIUS featureGroups
#' @importFrom checkmate assertClass testClass
#'
peakFilling <- function(object, algorithm = NULL, settings = NULL, save = TRUE) {

  checkmate::assertClass(object, "ntsData")

  if (is.null(algorithm)) algorithm <- fillingParameters(object)@algorithm

  if (is.null(settings)) settings <- fillingParameters(object)@settings

  if (algorithm == "xcms") {

    if (checkmate::testClass(settings, "list")) settings <- settings[[1]]

    Exp <- patRoon::getXCMSnExp(object@pat, loadRawData = TRUE)

    if (hasAdjustedRetentionTime(object)) {
      adjRT <- lapply(object@scans, function(x) x$adjustedRetentionTime)
      adjRT <- lapply(seq_len(length(adjRT)), function(x, adjRT) {
        temp <- adjRT[[x]]
        lenSeq <- seq_len(length(temp))
        names(temp) <- paste0(
          "F",
          x,
          ".S",
          sprintf(paste0("%0.", nchar(max(lenSeq)), "d"), lenSeq)
        )
        return(temp)
      }, adjRT = adjRT)

      names(adjRT) <- as.character(seq_len(length(adjRT)))

      xFD <- new("MsFeatureData")
      xFD$adjustedRtime <- adjRT
      xFD$chromPeakData <- Exp@msFeatureData$chromPeakData
      xFD$chromPeaks <- Exp@msFeatureData$chromPeaks
      xFD$featureDefinitions <- Exp@msFeatureData$featureDefinitions

      Exp@msFeatureData <- xFD
    }

    Exp <- xcms::fillChromPeaks(Exp, param = settings)

    pat <- patRoon::importFeatureGroupsXCMS3(Exp, patRoon::analysisInfo(object@pat))

    object@pat <- pat

    object <- buildPeaksTable(object)

    object <- buildFeatureTable(object)

    xFeats <- new("featuresSIRIUS",
      analysisInfo = pat@features@analysisInfo,
      features = pat@features@features
    )

    xFeatGroups <- new("featureGroupsSIRIUS",
      groups = pat@groups,
      groupInfo = pat@groupInfo,
      analysisInfo = pat@analysisInfo,
      features = xFeats,
      ftindex = pat@ftindex
    )

    object@pat <- xFeatGroups
  }

  object <- fillingParameters(object, algorithm = algorithm, settings = settings)

  if (is.logical(save)) if (save) saveObject(object = object)

  if (is.character(save)) saveObject(object = object, filename = save)

  return(object)
}
