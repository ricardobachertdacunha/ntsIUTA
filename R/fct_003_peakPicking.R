

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
#' @importFrom checkmate assertClass testChoice
#' @importClassesFrom patRoon features
#' @importFrom patRoon findFeatures featureTable
#' @importMethodsFrom MSnbase filterFile fileNames
#' @importMethodsFrom xcms chromPeaks
#' @importFrom dplyr select everything rename
#' @importFrom data.table rbindlist
#'
peakPicking <- function(obj = NULL,
                        algorithm = NULL,
                        param = NULL,
                        save = TRUE) {

  assertClass(obj, "ntsData")

  if (length(obj@MSnExp) == 0) {
    print("Importing data to MSnExp slot...")
    obj <- importRawData(obj, removeEmptySpectra = TRUE, save = FALSE)
  }
  
  if (is.null(algorithm)) algorithm <- obj@algorithms$peakPicking

  if (!testChoice(algorithm, c("xcms3", "xcms", "openms", "envipick", "sirius", "kpic2", "safd"))) {
    warning("Peak picking algorithm not recognized. See ?peakPicking for more information.")
    return(obj)
  }

  if (is.null(param)) param <- obj@parameters$peakPicking

  # if (length(param) == 0) {
  #   warning("Peak picking parameters not found or not given.")
  #   return(obj)
  # }

  sinfo <- data.frame(path = dirname(obj@samples$file),
                      analysis = obj@samples$sample,
                      group = obj@samples$group,
                      blank = obj@samples$blank)

  ag <- list(analysisInfo = sinfo, algorithm = algorithm)

  pat <- do.call(findFeatures, c(ag, param, verbose = TRUE))

  obj@patdata <- pat

  obj <- buildPeakList(obj)

  obj@algorithms$peakPicking <- algorithm

  obj@parameters$peakPicking <- param

  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}


#' @title buildPeakList
#'
#' @param obj An \linkS4class{ntsData} object containing a \linkS4class{features}
#' object in the slot \code{patdata}.
#'
#' @return A data.frame containing information for peaks in each sample.
#' 
#' @importClassesFrom patRoon features
#' @importMethodsFrom patRoon featureTable
#' @importFrom dplyr rename select
#' @importMethodsFrom xcms chromPeaks
#' @importFrom data.table rbindlist
#'
buildPeakList <- function(obj) {
  
  pat <- obj@patdata
  
  peaks <- base::as.data.frame(data.table::rbindlist(featureTable(pat), idcol = "sample"))
  
  if ("group" %in% colnames(peaks)) peaks <- rename(peaks, feature = group)
  
  peaks <- rename(peaks, rt = ret, rtmin = retmin, rtmax = retmax)
  
  peaks$group <- sapply(peaks$sample, FUN = function(x) x <- obj@samples$group[which(obj@samples$sample == x)])
  
  if (class(pat) == "featuresXCMS3" | class(pat) == "featureGroupsXCMS3") {
    extra <- base::as.data.frame(chromPeaks(pat@xdata, isFilledColumn = TRUE))
    extra <- rename(extra, intensity = maxo, area = into)
    peaks <- cbind(peaks, extra[,!(colnames(extra) %in% colnames(peaks))])
  }
  
  peaks$ID <- seq_len(nrow(peaks))
  
  rownames(peaks) <- peaks$ID
  
  peaks <- select(peaks, ID, sample, group, everything())
  
  obj@peaks <- peaks
  
  return(obj)
  
}
