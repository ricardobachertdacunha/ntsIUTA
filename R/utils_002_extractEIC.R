

#' @title extractEIC
#' @description Extracts the extracted ion chromatogram (EIC) from \code{rawData} of a specified \emph{m/z}.
#'
#' @param x An \linkS4class{OnDiskMSnExp} object with one or more files.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param mz Target \emph{m/z} to the EIC.
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' With Zero or \code{NULL} the complete retention time will be used.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{min}.
#' @param msLevel The MS level to extract the data. Possible values are 1 or 2.
#' @param normalize Logical, set to \code{TRUE} for normalizing the intensity values between 1 and 0 for each file in the given \code{x}.
#'
#' @return
#' 
#' @export
#' 
#' @import magrittr
#' @importFrom MSnbase filterFile rtime filterRt filterMz filterMsLevel
#' @importFrom methods as
#' @importFrom data.table rbindlist
#'
#' @examples
#' 
#' 
#' 
extractEIC <- function(x = rawData1, fileIndex = NULL,
                       mz = NULL, ppm = NULL,
                       rt = NULL, rtWindow = NULL,
                       rtUnit = "sec", msLevel = 1, normalize = FALSE) {
  
  #Examples
  # extractEIC(ntsIUTA::rawDataExample, fileIndex = 1, mz = 233.0243, ppm = 20)
  
  
  # require(MSnbase)
  # require(magrittr)
  # require(data.table)
  
  # x = rawData
  # fileIndex = 1:2
  # mz <- 233.0243
  # rt <- 0
  # rtUnit = "min"
  # ppm <- 20
  # rtWindow = 0
  # msLevel = 1
  
  if (rtUnit == "min") {rt <- rt*60}
  if (rtUnit == "min") {rtWindow <- rtWindow*60}
  
  if (!base::is.null(fileIndex)) {x <- MSnbase::filterFile(x, fileIndex)}
  
  rtr <- c(base::min(MSnbase::rtime(x)), base::max(MSnbase::rtime(x)))
  mzr <- NULL

  if (TRUE %in% (mz != 0 & !base::is.null(mz))) {
    
    if (base::length(mz) == 1) { mzr <- c(mz - ((ppm/1E6)*mz), mz + ((ppm/1E6)*mz)) }
    if (base::length(mz) == 2) { mzr <- c(mz[1], mz[2]) }
    
  } else {stop("m/z must be specified and different than 0.")}
  
  
  if (!base::is.null(rt)) { rtr <- c((rt) - base::ifelse(!base::is.null(rtWindow), rtWindow, 1),
                                     (rt) + base::ifelse(!base::is.null(rtWindow), rtWindow, 1)) }
  if (base::is.null(rt)) if (base::unique(!base::is.null(rtWindow))) if (base::length(rtWindow) == 2) { rtr <- c(rtWindow[1], rtWindow[2]) }
  
  x <- x %>% MSnbase::filterRt(rtr, msLevel. = c(1,2)) %>% MSnbase::filterMz(mzr, msLevel. = 1)
  
  x <- MSnbase::filterMsLevel(x, msLevel. = msLevel)
  if (base::length(x) == 0) stop("No MS data available")
  
  x <- base::suppressWarnings(methods::as(x, "data.frame"))
  x <- base::split(x, x$file)
  
  if (normalize)
  {
    #Make Standard Normalization Betwenn 1 and 0
    for (j in 1:length(x))
    {
      x[[j]]$i <- (x[[j]]$i - base::min(x[[j]]$i)) / (base::max(x[[j]]$i) - base::min(x[[j]]$i))
    }
  }
  
  #Export as csv
  x <- data.table::rbindlist(base::lapply(x, function(x) base::as.data.frame.list(x)), fill = TRUE)
  
  return(x)
  
}















