

#' @title filterPeaks
#'
#' @description Filter peaks from a \linkS4class{ntsData} object based on
#' file, retention time and emph{m/z}.
#' See \code{\link[patRoon]{filter}} for more information about the arguments.
#'
#' @param obj A \linkS4class{ntsData} object containing peaks.
#' @param mz The target emph{mz} or emph{mz} range, which should be given as a numeric vector with length 2.
#' @param ppm The mass deviation in ppm from the target emph{mz} to get peaks.
#' @param rt The target retention time to get peaks.
#' @param rtWindow The retention time deviation or the range (given as a numeric vector with length 2) to get peaks.
#' @param rtUnit The unit of the time arguments. Possible values are "sec" (the default) and "min".
#' @param absMinIntensity The absolute minimum intensity to obtain peaks.
#' @param chromWidthRange The min and max limits for chromatographic peak width given as a numeric vector with length 2.
#' Max can be \code{Inf} to specify no maximum limit.
#' @param negate Logical, set to \code{TRUE} to apply the opposite of the defined filters. The default is \code{FALSE}.
#' @param dataframe Logical, set to \code{TRUE} to return a data.frame instead of a \code{features} object. Usefull for inspection.
#'
#' @return A filtered \linkS4class{features} object or a \code{data.frame} when dataframe argument is set to \code{TRUE}.
#'
#' @references
#' \insertRef{Helmus2021}{ntsIUTA}
#'
#' @export
#'
#' @importClassesFrom patRoon features
#'
#' @examples
#'
filterPeaks <- function(peaks = peaks,
                        fileIndex = NULL,
                        mz = NULL,
                        ppm = NULL,
                        rt = NULL,
                        rtWindow = NULL,
                        rtUnit = "sec",
                        absMinIntensity = NULL,
                        chromWidthRange = NULL,
                        negate = FALSE,
                        dataframe = TRUE) {
  
  x <- peaks
  
  x <- x[fileIndex]
  
  mzr <- mzrBuilder(mz = mz, ppm = ppm)
  
  rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
  
  x <- patRoon::filter(obj = x,
                       absMinIntensity = absMinIntensity,
                       relMinIntensity = NULL,
                       retentionRange = rtr,
                       mzRange = mzr,
                       mzDefectRange = NULL,
                       chromWidthRange = chromWidthRange,
                       negate = negate)
  
  if (dataframe) x <- patRoon::as.data.frame(x)
  
  return(x)
  
}
