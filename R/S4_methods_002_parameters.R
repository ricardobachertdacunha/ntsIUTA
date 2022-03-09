

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
#' @param object An \linkS4class{ntsData} object.
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
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("fillingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@filling)


#' @describeIn ntsData Setter for fill missing parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character vector with the name of the algorithm to be applied.
#' @param settings A list with parameters matching the defined algorithm.
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
#' @param object An \linkS4class{ntsData} object.
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
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("annotationParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@annotation)


#' @describeIn ntsData Setter for annotation parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character vector with the name of the algorithm to be applied.
#' @param settings A list with parameters matching the defined algorithm.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("annotationParameters", c("ntsData", "character", "ANY"), function(object, algorithm, settings) {

  # if (!checkmate::testChoice(algorithm, c("xcms3", "xcms", "openms"))) {
  #   warning("Algorithm not recognized. See ?makeFeatures for more information.")
  #   return(object)
  # }

  if (!checkmate::testClass(settings, "list")) settings <- list(settings)

  object@parameters@annotation@algorithm <- algorithm

  object@parameters@annotation@settings <- settings

  return(object)
})
































#' @title AlteredCameraParam
#'
#' @param sigma The multiplier of the standard deviation for grouping features
#' by retention time.
#' @param perfwhm Percentage of the overlapping FWHM to group features.
#' @param cor_eic_th Threshold for feature EIC correlation in each sample.
#' @param cor_exp_th Threshold for intensity correlation across samples.
#' @param pval p-value threshold for testing correlation of significance.
#' @param validateIsotopePatterns Logical, set to \code{TRUE} for validating
#' the annotated isotopes with the \emph{kegg} database.
#' @param ppmIsotopes The expected mass deviation (in ppm) to find isotopes.
#' @param noise numeric.
#' @param searchAdducts Logical, set to \code{TRUE} to screen for adducts after finding isotopes.
#' @param ppmAdducts The expected mass deviation (in ppm) to find adducts.
#' @param extendedList Logical, set to \code{TRUE} to use the extended list of
#' adducts. The default is \code{FALSE}.
#'
#' @return An \linkS4class{AlteredCameraParam} class object for annotation.
#'
#' @export
#'
AlteredCameraParam <- function(
  sigma = 6,
  perfwhm = 0.4,
  cor_eic_th = 0.75,
  cor_exp_th = 0.75,
  pval = 0.1,
  validateIsotopePatterns = TRUE,
  ppmIsotopes = 40,
  noise = 300,
  searchAdducts = TRUE,
  ppmAdducts = 5,
  extendedList = TRUE) {

  paramobj <- do.call(new, c("AlteredCameraParam", as.list(environment())))

  paramobj

  return(paramobj)

}
