

#' peakPicking
#'
#' @description Finds chromatographic peaks from centroided MS data in the
#' samples of a given \linkS4class{ntsData} object.
#' The peak picking uses the \pkg{patRoon} package and
#' the following algorithms are available:
#' "xcms3", "xcms", "openms", "envipick", "sirius", "kpic2", "safd".
#' The parameters depend on the algorithm chosen.
#' See ?\pkg{patRoon} for further information.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character string with the algorithm to use for peak picking.
#' @param settings List of parameter settings for the specified algorithm.
#' See \link[patRoon]{findFeatures} for more information.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object with peaks per sample and an updated
#' \linkS4class{features} object added to the slot \code{pat}.
#'
#' @references
#' \insertRef{Helmus2021}{ntsIUTA}
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importClassesFrom patRoon features
#' @importFrom patRoon findFeatures featureTable
#'
peakPicking <- function(object = NULL,
                        algorithm = NULL,
                        settings = NULL,
                        save = TRUE) {

  assertClass(object, "ntsData")

  if (is.null(algorithm)) {
    algorithm <- pickingParameters(object)@algorithm
  }

  if (is.na(algorithm)) {
    warning("Peak picking algorihtm not defined!")
    return(object)
  }

  if (is.null(settings)) {
    settings <- pickingParameters(object)@settings
  }

  sinfo <- data.frame(path = dirname(object@samples$file),
                      analysis = object@samples$sample,
                      group = object@samples$replicate,
                      blank = object@samples$blank)

  sinfo$blank[is.na(sinfo$blank)] <- ""

  ag <- list(
    analysisInfo = sinfo,
    algorithm = algorithm
  )

  pat <- do.call(
    findFeatures,
    c(ag, settings,
    verbose = TRUE)
  )

  object@pat <- pat

  object <- buildPeaksTable(object)

  object <- pickingParameters(
    object,
    algorithm = algorithm,
    settings = settings
  )

  if (save) saveObject(object = object)

  if (is.character(save)) {
    saveObject(
      object = object,
      filename = save
    )
  }

  return(object)
}




#' @title buildPeaksTable
#'
#' @param object An \linkS4class{ntsData} object
#' containing a \linkS4class{features}
#' object in the slot \code{pat}.
#'
#' @return A \link[data.table]{data.table} containing
#' information for peaks in each sample.
#'
#' @importClassesFrom patRoon features
#' @importFrom dplyr select
#' @importMethodsFrom xcms chromPeaks
#' @importFrom data.table rbindlist as.data.table setnames
#' @importFrom checkmate testClass
#'
buildPeaksTable <- function(object) {

  cat("Building peaks table... \n")

  pat <- object@pat

  peaks <- pat@features

  if (checkmate::testClass(peaks, "features")) peaks <- peaks@features

  peaks <- rbindlist(peaks, idcol = "sample")

  peaks <- setnames(
    peaks,
    c("ID", "ret", "retmin", "retmax"),
    c("id", "rt", "rtmin", "rtmax"),
  )

  if (checkmate::testClass(pat, "featuresXCMS3") | checkmate::testClass(pat, "featureGroupsXCMS3")) {
    extra <- as.data.table(chromPeaks(pat@xdata, isFilledColumn = TRUE))
    extra <- setnames(extra, c("maxo", "into"), c("intensity", "area"))
    extra[, sample := samples(object)[extra$sample]]

    if (nrow(peaks) == nrow(extra) & all(peaks$mz == extra$mz)) {
      newCols <- colnames(extra)[!colnames(extra) %in% colnames(peaks)]
      peaks <- cbind(peaks, extra[, newCols, with = FALSE])
    }
  }

  peaks <- setnames(peaks, "group", "feature", skip_absent = TRUE)

  rpl <- replicates(object)
  names(rpl) <- samples(object)
  peaks[, replicate := rpl[peaks$sample]][]

  peaks$id <- seq_len(nrow(peaks))
  peaks$id <- paste0(
    "P",
    peaks$id,
    "_M",
    round(peaks$mz, digits = 0),
    "_R",
    round(peaks$rt, digits = 0)
  )

  peaks <- select(peaks,
    id,
    sample,
    replicate,
    mz, rt,
    intensity, area,
    mzmin, mzmax,
    rtmin, rtmax,
    everything()
  )

  object@peaks <- peaks

  return(object)
}
