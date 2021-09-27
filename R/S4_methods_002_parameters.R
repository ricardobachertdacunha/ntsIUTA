

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




#' MS2param
#'
#' @param maxMSRtWindow .
#' @param precursorMzWindow .
#' @param clusterMzWindow .
#' @param topMost .
#' @param minIntensityPre .
#' @param minIntensityPost .
#'
#' @return An \linkS4class{MS2param} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object.
#'
#' @export
#'
MS2param <- function(maxMSRtWindow = 10,
                     precursorMzWindow = 1.3,
                     clusterMzWindow = 0.003,
                     topMost = 50,
                     minIntensityPre = 10,
                     minIntensityPost = 10) {

  paramobj <- do.call(new, c("MS2param", as.list(environment())))

  paramobj

  return(paramobj)

}




### add parameters methods -----

#### peak picking -----

#' @describeIn ntsData Getter for peak picking parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("peakPickingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@peakPicking)


#' @describeIn ntsData Setter for peak picking parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character vector with the name of the algorithm to be applied.
#' @param param A list with parameters matching the defined algorithm.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("peakPickingParameters", c("ntsData", "character", "ANY"), function(object, algorithm, param) {

  if (!testChoice(algorithm, c("xcms3", "xcms", "openms", "envipick", "sirius", "kpic2", "safd"))) {
    warning("Peak picking algorithm not recognized. See ?peakPicking for more information.")
    return(object)
  }

  if (class(param) != "list") param <- list(param)

  object@parameters@peakPicking@algorithm <- algorithm
  
  object@parameters@peakPicking@param <- param

  return(object)

})

#### peak grouping -----

#' @describeIn ntsData Getter for peak grouping parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("peakGroupingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@peakGrouping)


#' @describeIn ntsData Setter for peak grouping parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character vector with the name of the algorithm to be applied.
#' @param param A list with parameters matching the defined algorithm.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("peakGroupingParameters", c("ntsData", "character", "ANY"), function(object, algorithm, param) {
  
  if (!testChoice(algorithm, c("xcms3", "xcms", "openms"))) {
    warning("Algorithm not recognized. See ?makeFeatures for more information.")
    return(object)
  }
  
  if (class(param) != "list") param <- list(param)
  
  object@parameters@peakGrouping@algorithm <- algorithm
  
  object@parameters@peakGrouping@param <- param
  
  return(object)
  
})

#### fill missing -----

#' @describeIn ntsData Getter for fill missing parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("fillMissingParameters", c("ntsData", "missing", "missing"), function(object) object@parameters@fillMissing)


#' @describeIn ntsData Setter for fill missing parameters.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param algorithm A character vector with the name of the algorithm to be applied.
#' @param param A list with parameters matching the defined algorithm.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("fillMissingParameters", c("ntsData", "character", "ANY"), function(object, algorithm, param) {
  
  if (!testChoice(algorithm, c("xcms3"))) {
    warning("Algorithm not recognized for filling features with missing peaks.
            See ?makeFeatures for more information.")
    return(object)
  }
  
  if (class(param) != "list") param <- list(param)
  
  object@parameters@fillMissing@algorithm <- algorithm
  
  object@parameters@fillMissing@param <- param
  
  return(object)
  
})

#### annotation -----

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
#' @param param A list with parameters matching the defined algorithm.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#'
setMethod("annotationParameters", c("ntsData", "character", "ANY"), function(object, algorithm, param) {
  
  if (!testChoice(algorithm, c("alteredcamera"))) {
    warning("Algorithm not recognized for annotation of features.
            See ?annotateFeatures for more information.")
    return(object)
  }
  
  if (class(param) != "list") param <- list(param)
  
  object@parameters@annotation@algorithm <- algorithm
  
  object@parameters@annotation@param <- param
  
  return(object)
  
})






