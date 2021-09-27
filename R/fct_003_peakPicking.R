

#' peakPicking
#'
#' @description Finds chromatographic peaks from centroided MS in the
#' samples of a given \linkS4class{ntsData} object.
#' The peak picking uses the \pkg{patRoon} package and
#' the following algorithms are possible:
#' "xcms3", "xcms", "openms", "envipick", "sirius", "kpic2", "safd".
#' The parameters depend on the algorithm chosen.
#' See ?\pkg{patRoon} for further information.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param algorithm The algorithm to use for peak picking.
#' @param param List of arguments for the specified algorithm.
#' See \code{\link[patRoon]{findFeatures}} for more information.
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

  if (is.null(algorithm)) algorithm <- peakPickingParameters(obj)@algorithm

  if (is.na(algorithm)) {
    warning("Peak picking algorihtm not defined!")
    return(obj)
  }
  
  if (is.null(param)) param <- peakPickingParameters(obj)@param

  sinfo <- data.frame(path = dirname(obj@samples$file),
                      analysis = obj@samples$sample,
                      group = obj@samples$group,
                      blank = obj@samples$blank)

  sinfo$blank[is.na(sinfo$blank)] <- ""
  
  ag <- list(analysisInfo = sinfo, algorithm = algorithm)

  pat <- do.call(findFeatures, c(ag, param, verbose = TRUE))

  obj@patdata <- pat

  obj <- buildPeakList(obj)

  obj <- peakPickingParameters(obj, algorithm = algorithm, param = param)
  
  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}




#' buildPeakList
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
