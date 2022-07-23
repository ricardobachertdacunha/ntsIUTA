

#### pickingParameters -----

#' @describeIn ntsData Getter for peak picking parameter settings.
#'
#' @export
#'
setMethod("pickingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@picking)


#' @describeIn ntsData Setter for peak picking parameter settings.
#'
#' @param algorithm A character vector with the name of the algorithm to be applied.
#' @param settings A list with parameter settings.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("pickingParameters", c("ntsData", "character", "ANY"), function(object, algorithm, settings) {

  if (!testChoice(algorithm, c("xcms3", "xcms", "openms", "envipick", "sirius", "kpic2", "safd"))) {
    warning("Peak picking algorithm not recognized. See ?peakPicking for more information.")
    return(object)
  }

  if (!checkmate::testClass(settings, "list")) settings <- list(settings)

  object@parameters@picking@algorithm <- algorithm

  object@parameters@picking@settings <- settings

  return(object)
})




#### groupingParameters -----

#' @describeIn ntsData Getter for peak grouping parameters.
#'
#' @export
#'
setMethod("groupingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@grouping)


#' @describeIn ntsData Setter for peak grouping parameters.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("groupingParameters", c("ntsData", "character", "ANY"), function(object, algorithm, settings) {

  if (!testChoice(algorithm, c("xcms3", "xcms", "openms"))) {
    warning("Algorithm not recognized. See ?makeFeatures for more information.")
    return(object)
  }

  if (!checkmate::testClass(settings, "list")) settings <- list(settings)

  object@parameters@grouping@algorithm <- algorithm

  object@parameters@grouping@settings <- settings

  return(object)
})




#### fillingParameters -----

#' @describeIn ntsData Getter for fill missing parameters.
#'
#' @export
#'
setMethod("fillingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@filling)


#' @describeIn ntsData Setter for fill missing parameters.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("fillingParameters", c("ntsData", "character", "ANY"), function(object, algorithm, settings) {

  if (missing(settings) & algorithm == "xcms3") {
    pfill <- fillingSettingsDefaultXCMS()
    settings <- pfill@settings
  }

  if (!checkmate::testClass(settings, "list")) settings <- list(settings)

  object@parameters@filling@algorithm <- algorithm

  object@parameters@filling@settings <- settings

  return(object)
})




#### fragmentsParameters -----

#' @describeIn ntsData Getter for fragments extraction parameters.
#'
#' @export
#'
setMethod("fragmentsParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@fragments)


#' @describeIn ntsData Setter for fragments extraction parameters.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("fragmentsParameters", c("ntsData", "character", "ANY"), function(object, algorithm, settings) {

  if (algorithm == "patroon") {
    ms2p <- patFragmentSettingsDefault()
    algorithm <- ms2p@algorithm
    settings <- ms2p@settings
  }

  if (algorithm == "ntsiuta" & missing(settings)) {
    ms2p <- fragmentSettingsDefault()
    settings <- ms2p@settings
  }

  if (!checkmate::testClass(settings, "list")) settings <- list(settings)

  object@parameters@fragments@algorithm <- algorithm

  object@parameters@fragments@settings <- settings

  return(object)
})




#### annotationParameters -----

#' @describeIn ntsData Getter for annotation parameters.
#'
#' @export
#'
setMethod("annotationParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@annotation)


#' @describeIn ntsData Setter for annotation parameters.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("annotationParameters", c("ntsData", "character", "ANY"), function(object, algorithm, settings) {

  if (!checkmate::testClass(settings, "list")) settings <- list(settings)

  object@parameters@annotation@algorithm <- algorithm

  object@parameters@annotation@settings <- settings

  return(object)
})




### Save and Load Parameters -----

#' @describeIn ntsData Method to save the \linkS4class{ntsParameters} in the project path
#' or other defined by the argument \code{path}.
#' A name for the \emph{rds} file can be specified by the argument \code{filename}, without format.
#'
#' @param path A character vector with a folder location.
#' @param filename A character string with a file name.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("saveParameters", "ntsData", function(object, path = NULL, filename = NULL) {

  if (missing(path)) {
    path <- path(object)
  } else if (is.null(path)) {
    path <- path(object)
  }

  if (missing(filename)) {
    filename <- "ntsParameters"
  } else if (is.null(path)) {
    filename <- "ntsParameters"
  }

  saveObject(parameters = object@parameters, path = path, filename = filename)
})


#' @describeIn ntsData Method to load an existing \linkS4class{ntsParameters} from disk.
#' The complete location of the file should be given by the argument \code{filepath}.
#' The default is the \link{choose.files} function.
#'
#' @param filepath A character string with the complete file path.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("loadParameters", "ntsData", function(object, filepath = NULL) {

  if (missing(filepath)) {
    filepath <- choose.files()
  } else if (is.null(filepath)) {
    filepath <- choose.files()
  }

  object@parameters <- readRDS(filepath)

  return(object)
})
