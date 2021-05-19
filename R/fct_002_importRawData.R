

#' @title importRawData
#' @description Function to import MS data from the mzML files listed in \code{sampleInfo} from the \code{setupProject} object.
#' Files should contain centroided spectra.
#' The function uses the function \code{\link[MSnbase]{readMSData}} from the \code{MSnbase} package to read the mzML files.
#'
#' @param sampleInfo The sample \code{data.table} obtained by the \code{\link{setupProject}} function.
#' Note that the \code{sampleInfo} can be filtered beforehand to exclude unwanted mzML files or to redefine the name of replicate sample groups.
#' @param rtFilter A numeric vector with length 2 defining the minimum and maximum chromatographic retention time for the listed files respectively.
#' If the files have different time dimensions, \code{rtFilter} should be \code{NULL}.
#' @param timeUnit If \code{rtFilter} is given in minutes please keep the default value \code{min} otherwise change to \code{sec} for seconds.
#' @param msLevel The MS dimensions for the rtFilter to be applied. Default is for both MS1 and MS2 with \code{c(1,2)}.
#' @param centroidedData Logical, set to \code{TRUE} for mzML files with centroided data or \code{FALSE} for profile data.
#' \code{NA} will collect all the data from the mzML file.
#' @param removeEmptySpectra Logical, set to TRUE if empty spectra should be removed.
#' It is recommended to remove empty spectra as it may cause issues during creation of features.
#' @param save Logical, set to \code{TRUE} to save the object rawData in the \strong{rData} folder.
#' If \strong{rData} folder does not exist in given \code{projPath} from \code{setup}, the folder will be created.
#' @param projPath If save is \code{TRUE}, the saving location or project path should be given.
#' The default is the \code{projPath} in the object \code{setup}.
#'
#' @return The output is a standard \linkS4class{OnDiskMSnExp} object,
#' which can be used within the workflow of \pkg{ntsIUTA} but also within \pkg{Bioconductor} packages, such as \pkg{xcms}.
#'
#' @references
#' \insertRef{MSnbase1}{ntsIUTA}
#' \insertRef{MSnbase2}{ntsIUTA}
#'
#' @export
#'
#' @importFrom MSnbase readMSData filterRt filterEmptySpectra
#'
#' @examples
#' sampleInfo <- ntsIUTA::setupProject(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)
#' rawData <- ntsIUTA::importRawData(sampleInfo[1,], save = FALSE, centroidedData = TRUE)
#' rawData
#'
importRawData <- function(sampleInfo = sampleInfo,
                          rtFilter = NULL,
                          timeUnit = "min",
                          msLevel = c(1, 2),
                          centroidedData = TRUE,
                          removeEmptySpectra = TRUE,
                          save = FALSE,
                          projPath = base::getwd()) {

  rawData <- base::suppressWarnings(
    MSnbase::readMSData(sampleInfo$filePath[drop = TRUE],
                        pdata = methods::new("NAnnotatedDataFrame",
                        base::data.frame(sample_name = sampleInfo$sample,
                                         sample_group = sampleInfo$group)),
                        mode = "onDisk",
                        centroided. = centroidedData,
                        smoothed. = FALSE)
  )

  if (!base::is.null(rtFilter)) {
    if (timeUnit == "min") rtFilter <- rtFilter * 60
    rawData <- MSnbase::filterRt(rawData, rt = rtFilter, msLevel. = msLevel)
  }

  if (removeEmptySpectra) rawData <- base::suppressWarnings(MSnbase::filterEmptySpectra(rawData))

  if (save) {
    rData <- base::paste0(projPath, "\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(rawData, file = base::paste0(rData, "\\rawData.rds"))
  }

  return(rawData)

}



#' @title centroidProfileData
#' @description Centroiding of profile data with additional possibility for data smoothing before centroiding and \emph{m/z} refinement.
#' The \code{centroidProfileData} function combines functions \code{smooth} and \code{pickPeaks}
#' from the \code{MSnbase} package, see references.
#'
#' @param x A \linkS4class{OnDiskMSnExp} object with profile data for centroiding.
#' @param halfwindow Sets the window size for centroiding as \code{2 * halfwindow + 1}.
#' The \code{halfwindow} should be slightly larger than the full width at half maximum of the profile peak.
#' @param SNR The signal-to-noise ratio to consider a local maximum as peak.
#' @param noiseMethod The method to estimate the noise level. Possible methods are "MAD" (the default) and "SuperSmoother".
#' See \code{?MSnbase::pickPeaks} for more information.
#' @param smoothing Logical, set to \code{TRUE} for applying smothing to the profile data before centroiding. The default is FALSE.
#' @param methodSmoothing Method for data smoothing. The possible methods are "SavitzkyGolay" (the default) and "MovingAverage".
#' See \code{?MSnbase::smooth} for more information and arguments, which are passed by \code{...}.
#' @param ... Arguments for selected smoothing method. See \code{?MSnbase::smooth} for possible arguments for each method.
#' @param methodRefineMz Method for refinement. Possible methods are "none" (the default, for not applying \emph{m/z} refinement), "kNeighbors" and "descendPeak".
#' See \code{?MSnbase::pickPeaks} for more information.
#' @param k When refine method is "kNeighbors", \code{k} is number of closest signals to the centroid.
#' @param signalPercentage When refine method is "descendPeak", \code{signalPercentage} is the minimum signal percentage of centroids to refine \emph{m/z}.
#' @param stopAtTwo Logical, when refine method is "descendPeak", set to \code{TRUE} for allowing two consecutive equal or higher signals.
#' \code{FALSE} will stop when one equal or higher centroid is found.
#' @param save Logical, set to \code{TRUE} to replace the original files by the centroided files in disk.
#' The location is taken from the originbal file paths.
#'
#' @return Centroided \linkS4class{OnDiskMSnExp} object.
#' When \code{save} is set to TRUE, the profile data in the original mzML or mzXML files are replaced by the centroided data.
#'
#' @references
#' \insertRef{MSnbase2}{ntsIUTA}
#'
#' @export
#'
#' @importFrom MSnbase smooth pickPeaks fileNames writeMSData
#'
#' @examples
#'
centroidProfileData <- function(x,
                                halfwindow = 2,
                                SNR = 0,
                                noiseMethod = "MAD",
                                smoothing = FALSE,
                                methodSmoothing = "SavitzkyGolay",
                                methodRefineMz = "kNeighbors",
                                k = 1,
                                signalPercentage = 10, stopAtTwo = TRUE,
                                save = FALSE, ...) {

  if (smoothing) {
    x <- x %>% MSnbase::smooth(method = methodSmoothing, ...)
  }

  if (methodRefineMz == "kNeighbors") {
    x <- MSnbase::pickPeaks(x,
                            halfWindowSize = halfwindow,
                            SNR = SNR,
                            noiseMethod = noiseMethod,
                            refineMz = methodRefineMz,
                            k = k)
  } else {
    if (methodRefineMz == "descendPeak") {
      x <- MSnbase::pickPeaks(x,
                              halfWindowSize = halfwindow,
                              SNR = SNR,
                              noiseMethod = noiseMethod,
                              refineMz = methodRefineMz,
                              signalPercentage = 0.1,
                              stopAtTwo = TRUE)
    } else {
      x <- MSnbase::pickPeaks(x,
                              halfWindowSize = halfwindow,
                              SNR = SNR,
                              noiseMethod = noiseMethod,
                              refineMz = "none")
    }
  }

  if (save) {
    fls_new <- MSnbase::fileNames(x)
    MSnbase::writeMSData(x, file = fls_new)
  }

  return(x)

}
