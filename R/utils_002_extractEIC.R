

#' @title extractEIC
#' @description Extracts the extracted ion chromatogram (EIC) from \code{rawData} of a specified \emph{m/z}.
#'
#' @param raw An \linkS4class{OnDiskMSnExp} object with one or more files.
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
#' @importMethodsFrom ProtGenerics filterMz
#' @importMethodsFrom MSnbase filterFile filterRt filterMz filterMsLevel rtime
#' @importFrom methods as
#' @importFrom data.table rbindlist
#'
#' @examples
#'
extractEIC <- function(raw = rawData, fileIndex = NULL,
                       mz = NULL, ppm = NULL,
                       rt = NULL, rtWindow = NULL,
                       rtUnit = "sec", msLevel = 1,
                       normIntensity = FALSE) {

  #Examples
  # extractEIC(ntsIUTA::rawDataExample, fileIndex = 1, mz = 233.0243, ppm = 20)

  # raw <- rawData
  # fileIndex <- 4:5
  # mz <- 233.0243
  # rt <- NULL
  # rtUnit <- "min"
  # ppm <- 20
  # rtWindow <- NULL
  # msLevel <- 1

  if (rtUnit == "min") if (!is.null(rt)) rt <- rt * 60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow * 60

  if (!is.null(fileIndex)) {
    raw <- filterFile(raw, fileIndex)
  }

  mzr <- NULL

  if (!is.null(mz)) {
    if (length(mz) == 1) {
      if (is.null(ppm)) ppm <- 20
      mzr <- c(mz - ((ppm / 1E6) * mz), mz + ((ppm / 1E6) * mz))
    }
    if (length(mz) == 2) {
      mzr <- c(mz[1], mz[2])
    }
  }

  rtr <- c(min(rtime(raw)), max(rtime(raw)))

  if (!is.null(rt)) {
    rtr <- c((rt) - ifelse(!is.null(rtWindow), rtWindow, 60),
             (rt) + ifelse(!is.null(rtWindow), rtWindow, 60))
  }

  if (is.null(rt)) if (unique(!is.null(rtWindow))) if (length(rtWindow) == 2) {
    rtr <- c(rtWindow[1], rtWindow[2])
  }

  raw <- filterRt(object = raw, rt = rtr, msLevel. = c(1, 2))

  raw <- filterMsLevel(raw, msLevel. = msLevel)

  if (!is.null(mzr)) {
    raw <- base::suppressWarnings(filterMz(raw, mz = mzr, msLevel. = msLevel))
  }

  raw <- base::suppressWarnings(methods::as(raw, "data.frame"))

  if (normIntensity) {
    raw <- split(raw, raw$file)
    for (j in seq_len(length(raw))) {
      raw[[j]]$i <- (raw[[j]]$i - min(raw[[j]]$i)) / (max(raw[[j]]$i) - min(raw[[j]]$i))
    }
    raw <- data.table::rbindlist(lapply(raw, function(x) as.data.frame.list(x)), fill = TRUE)
  }

  return(raw)

}
