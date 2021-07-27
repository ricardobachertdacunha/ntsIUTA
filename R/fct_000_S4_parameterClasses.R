

### AlteredCameraParam -----

#' @title AlteredCameraParam
#'
#' @slot sigma The multiplier of the standard deviation for grouping features
#' by retention time.
#' @slot perfwhm Percentage of the overlapping FWHM to group features.
#' @slot cor_eic_th Threshold for feature EIC correlation in each sample.
#' @slot calcCaS Calculate correlation across samples.
#' @slot cor_exp_th Threshold for intensity correlation across samples.
#' @slot calcIso Logical, calculates based on isotopic information.
#' @slot pval p-value threshold for testing correlation of significance.
#' @slot validateIsotopePatterns Logical, set to \code{TRUE} for validating
#' the annotated isotopes with the \emph{kegg} database.
#' @slot ppmIsotopes The expected mass deviation (in ppm) to find isotopes.
#' @slot mzabs numeric.
#' @slot noise numeric.
#' @slot searchAdducts Logical, set to \code{TRUE} to screen for adducts after finding isotopes.
#' @slot ppmAdducts The expected mass deviation (in ppm) to find adducts.
#' @slot extendedList Logical, set to \code{TRUE} to use the extended list of
#' adducts. The default is \code{FALSE}.
#'
#' @return An \linkS4class{AlteredCameraParam} object containing parameters for
#' annotation of features in an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("AlteredCameraParam",
  slots = c(
    sigma = "numeric",
    perfwhm = "numeric",
    cor_eic_th = "numeric",
    calcCaS = "logical",
    cor_exp_th = "numeric",
    calcIso = "logical",
    pval = "numeric",
    validateIsotopePatterns = "logical",
    ppmIsotopes = "numeric",
    mzabs = "numeric",
    noise = "numeric",
    searchAdducts = "logical",
    ppmAdducts = "numeric",
    extendedList = "logical"
))



#' @title AlteredCameraParam
#'
#' @param sigma The multiplier of the standard deviation for grouping features
#' by retention time.
#' @param perfwhm Percentage of the overlapping FWHM to group features.
#' @param cor_eic_th Threshold for feature EIC correlation in each sample.
#' @param calcCaS Calculate correlation across samples.
#' @param cor_exp_th Threshold for intensity correlation across samples.
#' @param pval p-value threshold for testing correlation of significance.
#' @param validateIsotopePatterns Logical, set to \code{TRUE} for validating
#' the annotated isotopes with the \emph{kegg} database.
#' @param ppmIsotopes The expected mass deviation (in ppm) to find isotopes.
#' @param mzabs numeric.
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
#' @examples
#'
AlteredCameraParam <- function(
  sigma = 5,
  perfwhm = 0.45,
  cor_eic_th = 0.85,
  calcCaS = TRUE,
  cor_exp_th = 0.85,
  pval = 0.05,
  validateIsotopePatterns = TRUE,
  ppmIsotopes = 50,
  mzabs = 0.01,
  noise = 350,
  searchAdducts = TRUE,
  ppmAdducts = 5,
  extendedList = FALSE) {

  paramobj <- do.call(new, c("AlteredCameraParam", as.list(environment())))

  paramobj

  return(paramobj)

}



### MS2param -----

#' @title MS2param
#'
#' @slot maxMSRtWindow .
#' @slot precursorMzWindow .
#' @slot clusterMzWindow .
#' @slot topMost .
#' @slot minIntensityPre .
#' @slot minIntensityPost .
#'
#' @return An \linkS4class{MS2param} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("MS2param",
  slots = c(
    maxMSRtWindow = "numeric",
    precursorMzWindow = "numeric",
    clusterMzWindow = "numeric",
    topMost = "numeric",
    minIntensityPre = "numeric",
    minIntensityPost = "numeric"
  ),
  prototype = list(
    maxMSRtWindow = 10,
    precursorMzWindow = 1.3,
    clusterMzWindow = 0.003,
    topMost = 50,
    minIntensityPre = 10,
    minIntensityPost = 10
  )
)



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
