

#' @title extractEIC
#' @description Extracts an ion chromatogram (EIC) from raw data of a specified \emph{m/z}.
#'
#' @param obj An \linkS4class{ntsData} object with one or more files.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param mz Target \emph{m/z} to the EIC.
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' With Zero or \code{NULL} the complete retention time will be used.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{sec}.
#' @param msLevel The MS level to extract the data. For the momment, only 1 is possible.
#' @param normIntensity Logical, set to \code{TRUE} for normalizing the intensity values between 1 and 0 for each file in the given \code{raw} object.
#'
#' @return A \code{data.frame} with the columns \code{file}, \code{rt}, \code{mz} and \code{i}
#' representing the mzML file index, the retention time, the \emph{m/z} and the intensity, respectively.
#'
#' @export
#'
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importClassesFrom ProtGenerics ProcessingStep
#' @importFrom ProtGenerics ProcessingStep executeProcessingStep
#' @importMethodsFrom ProtGenerics filterMz rtime
#' @importMethodsFrom MSnbase filterFile filterRt filterMz filterMsLevel rtime
#' @importFrom methods as
#' @importFrom data.table rbindlist
#'
#' @examples
#'
extractEIC <- function(obj = NULL, fileIndex = NULL,
                       mz = NULL, ppm = NULL,
                       rt = NULL, rtWindow = NULL,
                       rtUnit = "sec", msLevel = 1,
                       normIntensity = FALSE) {

  # raw <- rawData
  # fileIndex <- 4:5
  # mz <- 233.0243
  # rt <- NULL
  # rtUnit <- "min"
  # ppm <- 20
  # rtWindow <- NULL
  # msLevel <- 1
  

  if (!is.null(fileIndex)) obj <- obj[fileIndex]

  mzr <- NULL

  if (!is.null(mz)) mzr <- mzrBuilder(mz = mz, ppm = ppm)

  rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)

  if (is.null(rtr)) {
    tt <- rtime(obj@MSnExp)
    rtr <- c(min(tt), max(tt))
  }

  raw <- filterRt(object = obj@MSnExp, rt = rtr)

  raw <- filterMsLevel(raw, msLevel. = msLevel)

  if (!is.null(mzr)) raw <- suppressWarnings(filterMz(raw, mz = mzr))

  raw <- suppressWarnings(methods::as(raw, "data.frame"))

  if (normIntensity) {
    raw <- split(raw, raw$file)
    for (j in seq_len(length(raw))) {
      raw[[j]]$i <- (raw[[j]]$i - min(raw[[j]]$i)) / (max(raw[[j]]$i) - min(raw[[j]]$i))
    }
    raw <- rbindlist(lapply(raw, function(x) as.data.frame.list(x)), fill = TRUE)
  }

  return(raw)

}
