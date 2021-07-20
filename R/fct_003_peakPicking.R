

#' @title peakPicking
#' @description Finds chromatographic peaks from centroided MS data
#' in an \linkS4class{ntsData} object.
#' The \code{peakPicking} function can use both \pkg{xcms}
#' and \pkg{patRoon} packages, depending on the specified algorithm.
#' Possible algorithms are: "xcms3", "xcms", "openms" and "envipick".
#' The function arguments depend on the algorithm chosen.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param algorithm The algorithm to use for peak picking.
#' One of "xcms3" (the default), "xcms", "openms" or "envipick".
#' @param param List of arguments for the specified algorithm.
#' See \code{\link[patRoon]{findFeatures}} for more information.
#' For instance, for "xcms3" as \code{algorithm}, \code{param} should be given
#' with the method parameters for the peak finding.
#' See \href{https://rdrr.io/bioc/xcms/man/chromatographic-peak-detection.html}{\code{xcms::findChromPeaks}}
#' for more information.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object containing a \linkS4class{features}
#' object in the slot \code{patdata}.
#'
#' @references
#' \insertRef{Helmus2021}{ntsIUTA}
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#'
#' @export
#'
#' @importClassesFrom patRoon features
#' @importFrom patRoon findFeatures featureTable
#' @importMethodsFrom MSnbase filterFile fileNames
#' @importMethodsFrom xcms chromPeaks
#' @importFrom dplyr select everything rename
#' @importFrom data.table rbindlist
#'
#' @examples
#'
peakPicking <- function(obj = NULL,
                        algorithm = NULL,
                        param = NULL,
                        save = TRUE) {
  
  if (is.null(algorithm)) algorithm <- obj@algorithms$peakPicking
  
  checkmate::checkChoice(algorithm, c("xcms3", "xcms", "openms", "envipick"))
  
  if (is.null(param)) param <- obj@parameters$peakPicking
  
  if (length(param) == 0) return(print("Peak picking parameters not found or not given."))
  
  sinfo <- data.frame(path = dirname(obj@samples$file),
                      analysis = obj@samples$sample,
                      group = obj@samples$group,
                      blank = obj@samples$blank)

  ag <- list(analysisInfo = sinfo, algorithm = algorithm)

  pat <- do.call(findFeatures, c(ag, param, verbose = TRUE))
  
  obj@patdata <- pat
  
  if (class(pat) == "featuresXCMS3") {
    peaks <- as.data.frame(chromPeaks(pat@xdata, isFilledColumn = TRUE))
    peaks$group <- sapply(peaks$sample, FUN = function(x) x <- obj@samples$group[x])
    peaks$sample <- sapply(peaks$sample, FUN = function(x) x <- obj@samples$sample[x])
    peaks <- split(peaks, f = peaks$sample)
    peaks <- as.data.frame(rbindlist(peaks))
    peaks$ID <- 1:nrow(peaks)
    peaks <- select(peaks, ID, sample, group, everything())
    peaks <- rename(peaks, intensity = maxo, area = into)
  } else {
    peaks <- featureTable(pat)
    for (i in seq_len(nrow(obj@samples))) peaks[[i]]$sample <- obj@samples$sample[i]
    for (i in seq_len(nrow(obj@samples))) peaks[[i]]$group <- obj@samples$group[i]
    peaks <- as.data.frame(rbindlist(peaks))
    peaks$ID <- 1:nrow(peaks)
    peaks <- select(peaks, ID, sample, group, everything())
    peaks <- rename(peaks, rt = ret, rtmin = retmin, rtmax = retmax)
  }
  
  obj@peaks <- peaks
  
  obj@algorithms$peakPicking <- algorithm
  
  obj@parameters$peakPicking <- param
  
  if (save) saveObject(obj = obj)
  
  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}
