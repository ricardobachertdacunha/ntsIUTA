

#' @title importRawData
#' @description Function to import MS data from the MS files listed in
#' the \linkS4class{ntsData} object. Files should contain centroided spectra.
#' The function \code{\link[MSnbase]{readMSData}} from the
#' \code{MSnbase} package is used to read the MS files.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param rtFilter A numeric vector with length 2 defining the minimum
#' and maximum chromatographic retention time for the listed MS files.
#' @param rtUnit The unit of the \code{rtFilter}.
#' Possible values are \code{min} (the default) and \code{sec}.
#' @param msLevel The MS dimensions for the rtFilter to be applied.
#' The default is both MS1 and MS2 using \code{c(1,2)}.
#' @param centroidedData Logical, set to \code{TRUE} for MS files
#' with centroided data or \code{FALSE} for profile data.
#' \code{NA} will collect all the data from the MS files.
#' @param removeEmptySpectra Logical, set to TRUE if empty spectra should be removed.
#' It is recommended to remove empty spectra as it may cause issues during creation of features.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return The \linkS4class{ntsData} object including a standard
#' \linkS4class{OnDiskMSnExp} object in the MSnExp slot. Note, that
#' the \linkS4class{OnDiskMSnExp} object can also be used within
#' the workflow of \pkg{Bioconductor} packages.
#'
#' @references
#' \insertRef{MSnbase1}{ntsIUTA}
#' \insertRef{MSnbase2}{ntsIUTA}
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom BiocParallel SnowParam register
#' @importFrom parallel detectCores
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importFrom MSnbase readMSData
#' @importMethodsFrom MSnbase filterRt filterEmptySpectra
#' @importFrom methods new
#'
#' @examples
#' path <- system.file(package = "ntsIUTA", dir = "extdata")
#' dt <- setupProject(path = path, save = FALSE)
#' dt <- importRawData(dt[1], save = FALSE, centroidedData = TRUE)
#'
importRawData <- function(obj = NULL,
                          rtFilter = NULL,
                          rtUnit = "min",
                          msLevel = c(1, 2),
                          centroidedData = TRUE,
                          removeEmptySpectra = TRUE,
                          save = FALSE) {

  assertClass(obj, "ntsData")

  snow <- SnowParam(workers = detectCores() - 1,
                    type = "SOCK",
                    exportglobals = FALSE,
                    progressbar = TRUE)

  register(snow, default = TRUE)

  msFiles <- obj@samples$file[drop = TRUE]
  sample_name <- obj@samples$sample
  sample_group <- obj@samples$group

  if (length(sample_name) == 0) {
    warning("There are not samples in the ntsData object.")
    return(obj)
  }

  raw <- suppressWarnings(
    readMSData(msFiles,
               pdata = new("NAnnotatedDataFrame",
                data.frame(sample_name = sample_name,
                           sample_group = sample_group)),
               msLevel. = NULL,
               mode = "onDisk",
               centroided. = centroidedData,
               smoothed. = FALSE)
  )

  if (!is.null(rtFilter)) {
    if (rtUnit == "min") rtFilter <- rtFilter * 60
    raw <- filterRt(raw, rt = rtFilter, msLevel. = msLevel)
  }

  if (removeEmptySpectra) raw <- filterEmptySpectra(raw)

  obj@MSnExp <- raw

  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}




#' @title centroidProfileData
#' @description Centroiding of profile data with additional possibility
#' for data smoothing before centroiding and \emph{m/z} refinement.
#' The \code{centroidProfileData} function combines functions \code{smooth}
#' and \code{pickPeaks} from the \code{MSnbase} package, see references.
#'
#' @param obj A \linkS4class{ntsData} object with profile data for centroiding.
#' @param halfwindow Sets the window size for centroiding as \code{2 * halfwindow + 1}.
#' The \code{halfwindow} should be slightly larger than the full width
#' at half maximum of the profile peak.
#' @param SNR The signal-to-noise ratio to consider a local maximum as peak.
#' @param noiseMethod The method to estimate the noise level.
#' Possible methods are "MAD" (the default) and "SuperSmoother".
#' See \code{?MSnbase::pickPeaks} for more information.
#' @param smoothing Logical, set to \code{TRUE} for applying smothing
#' to the profile data before centroiding. The default is FALSE.
#' @param methodSmoothing Method for data smoothing.
#' The possible methods are "SavitzkyGolay" (the default) and "MovingAverage".
#' See \code{?MSnbase::smooth} for more information and arguments,
#' which are passed by \code{...}.
#' @param ... Arguments for selected smoothing method.
#' See \code{?MSnbase::smooth} for possible arguments for each method.
#' @param methodRefineMz Method for refinement.
#' Possible methods are "none" (the default, for not applying \emph{m/z} refinement),
#' "kNeighbors" and "descendPeak". See \code{?MSnbase::pickPeaks} for more information.
#' @param k When refine method is "kNeighbors",
#' \code{k} is number of closest signals to the centroid.
#' @param signalPercentage When refine method is "descendPeak",
#' \code{signalPercentage} is the minimum signal percentage of centroids to refine \emph{m/z}.
#' @param stopAtTwo Logical, when refine method is "descendPeak",
#' set to \code{TRUE} for allowing two consecutive equal or higher signals.
#' \code{FALSE} will stop when one equal or higher centroid is found.
#' @param save Logical, set to \code{TRUE} to replace
#' the original files by the centroided files in disk.
#' The location is taken from the originbal file paths.
#'
#' @return Centroided \linkS4class{ntsData} object.
#' When \code{save} is set to TRUE, the profile data in the original
#' mzML or mzXML files is replaced by the centroided data.
#'
#' @references
#' \insertRef{MSnbase2}{ntsIUTA}
#'
#' @export
#'
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importMethodsFrom MSnbase fileNames smooth pickPeaks writeMSData
#'
centroidProfileData <- function(obj,
                                halfwindow = 2,
                                SNR = 0,
                                noiseMethod = "MAD",
                                smoothing = FALSE,
                                methodSmoothing = "SavitzkyGolay",
                                methodRefineMz = "kNeighbors",
                                k = 1,
                                signalPercentage = 10, stopAtTwo = TRUE,
                                save = FALSE, ...) {

  raw <- obj@MSnExp

  if (smoothing) {
    raw <- raw %>% MSnbase::smooth(method = methodSmoothing, ...)
  }

  if (methodRefineMz == "kNeighbors") {
    raw <- pickPeaks(raw,
                     halfWindowSize = halfwindow,
                     SNR = SNR,
                     noiseMethod = noiseMethod,
                     refineMz = methodRefineMz,
                     k = k)
  } else {
    if (methodRefineMz == "descendPeak") {
      raw <- pickPeaks(raw,
                       halfWindowSize = halfwindow,
                       SNR = SNR,
                       noiseMethod = noiseMethod,
                       refineMz = methodRefineMz,
                       signalPercentage = signalPercentage,
                       stopAtTwo = TRUE)
    } else {
      raw <- pickPeaks(raw,
                       halfWindowSize = halfwindow,
                       SNR = SNR,
                       noiseMethod = noiseMethod,
                       refineMz = "none")
    }
  }

  obj@MSnExp <- raw

  if (save) {
    fls_new <- fileNames(raw)
    writeMSData(raw, file = fls_new)
  }

  return(obj)

}




#' @title extractEIC
#' @description Extracts an ion chromatogram (EIC) from raw data of a specified \emph{m/z}.
#'
#' @param obj An \linkS4class{ntsData} object with one or more files.
#' @param samples The index or names of the sample/s to extract the centroids or profile data.
#' @param mz Target \emph{m/z} to the EIC.
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' With Zero or \code{NULL} the complete retention time will be used.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{sec}.
#' @param msLevel The MS level to extract the data. For the moment, only 1 is possible.
#' @param normIntensity Logical, set to \code{TRUE} for normalizing the intensity values between 1 and 0 for each file in the given \code{raw} object.
#'
#' @return A \code{data.frame} with the columns \code{file}, \code{rt}, \code{mz} and \code{i}
#' representing the mzML file index, the retention time, the \emph{m/z} and the intensity, respectively.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importClassesFrom ProtGenerics ProcessingStep
#' @importFrom ProtGenerics ProcessingStep executeProcessingStep
#' @importMethodsFrom ProtGenerics filterMz rtime
#' @importMethodsFrom MSnbase filterFile filterRt filterMz filterMsLevel rtime
#' @importFrom methods as
#' @importFrom data.table rbindlist
#'
extractEIC <- function(obj = NULL, samples = NULL,
                       mz = NULL, ppm = NULL,
                       rt = NULL, rtWindow = NULL,
                       rtUnit = "sec", msLevel = 1,
                       normIntensity = FALSE) {

  assertClass(obj, "ntsData")

  if (length(obj@MSnExp) == 0) {
    warning("No raw data found in the MSnExp slot of the ntsData object.
            Use the function importRawData to import data from raw files")
  }
  
  if (is.character(samples)) {
    if (FALSE %in% (samples %in% obj@samples$sample)) {
      warning("Given sample names not found in the ntsData object!")
    }
    samples <- which(obj@samples$sample %in% samples)
  }
  
  if (!is.null(samples)) obj@MSnExp <- filterFile(obj@MSnExp, file = samples)

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
