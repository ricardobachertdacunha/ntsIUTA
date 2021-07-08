

#' @title importRawData
#' @description Function to import MS data from the MS files listed in
#' the \linkS4class{ntsData} object. Files should contain centroided spectra.
#' The function \code{\link[MSnbase]{readMSData}} from the
#' \code{MSnbase} package is used to read the MS files.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param rtFilter A numeric vector with length 2 defining the minimum
#' and maximum chromatographic retention time for the listed MS files.
#' @param timeUnit The unit of the \code{rtFilter}.
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
#' dt <- importRawData(obj[1], save = FALSE, centroidedData = TRUE)
#'
importRawData <- function(obj = NULL,
                          rtFilter = NULL,
                          timeUnit = "min",
                          msLevel = c(1, 2),
                          centroidedData = TRUE,
                          removeEmptySpectra = TRUE,
                          save = FALSE) {
  
  snow <- SnowParam(workers = detectCores() - 1,
                    type = "SOCK",
                    exportglobals = FALSE,
                    progressbar = TRUE)
  
  register(snow, default = TRUE)
  
  msFiles <- obj@samples$file[drop = TRUE]
  sample_name <- obj@samples$sample
  sample_group <- obj@samples$group
  
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
    if (timeUnit == "min") rtFilter <- rtFilter * 60
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
#' @examples
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



#' concatenate_OnDiskMSnExp
#'
#' @param ... A list of OnDiskMSnExp objects to concatenate using \code{c()}.
#'
#' @return An OnDiskMSnExp object.
concatenate_OnDiskMSnExp <- function(...) {
  x <- list(...)
  if (length(x) == 0)
    return(NULL)
  if (length(x) == 1)
    return(x[[1]])
  ## Check that all are XCMSnExp objects.
  if (!all(unlist(lapply(x, function(z) is(z, "OnDiskMSnExp")))))
    stop("All passed objects should be 'OnDiskMSnExp' objects")
  ## Check processingQueue
  procQ <- lapply(x, function(z) z@spectraProcessingQueue)
  new_procQ <- procQ[[1]]
  is_ok <- unlist(lapply(procQ, function(z)
    !is.character(all.equal(new_procQ, z))
  ))
  if (any(!is_ok)) {
    warning("Processing queues from the submitted objects differ! ",
            "Dropping the processing queue.")
    new_procQ <- list()
  }
  ## processingData
  fls <- lapply(x, function(z) z@processingData@files)
  startidx <- cumsum(lengths(fls))
  ## featureData
  featd <- lapply(x, fData)
  ## Have to update the file index and the spectrum names.
  for (i in 2:length(featd)) {
    featd[[i]]$fileIdx <- featd[[i]]$fileIdx + startidx[i - 1]
    rownames(featd[[i]]) <- MSnbase:::formatFileSpectrumNames(
      fileIds = featd[[i]]$fileIdx,
      spectrumIds = featd[[i]]$spIdx,
      nSpectra = nrow(featd[[i]]),
      nFiles = length(unlist(fls))
    )
  }
  featd <- do.call(rbind, featd)
  featd$spectrum <- 1:nrow(featd)
  ## experimentData
  expdata <- lapply(x, function(z) {
    ed <- z@experimentData
    data.frame(instrumentManufacturer = ed@instrumentManufacturer,
               instrumentModel = ed@instrumentModel,
               ionSource = ed@ionSource,
               analyser = ed@analyser,
               detectorType = ed@detectorType,
               stringsAsFactors = FALSE)
  })
  expdata <- do.call(rbind, expdata)
  expdata <- new("MIAPE",
                 instrumentManufacturer = expdata$instrumentManufacturer,
                 instrumentModel = expdata$instrumentModel,
                 ionSource = expdata$ionSource,
                 analyser = expdata$analyser,
                 detectorType = expdata$detectorType)
  
  ## protocolData
  protodata <- lapply(x, function(z) z@protocolData)
  if (any(unlist(lapply(protodata, nrow)) > 0))
    warning("Found non-empty protocol data, but merging protocol data is",
            " currently not supported. Skipped.")
  ## phenoData
  pdata <- do.call(rbind, lapply(x, pData))
  res <- new(
    "OnDiskMSnExp",
    phenoData = new("NAnnotatedDataFrame", data = pdata),
    featureData = new("AnnotatedDataFrame", featd),
    processingData = new("MSnProcess",
                         processing = paste0("Concatenated [", date(), "]"),
                         files = unlist(fls), smoothed = NA),
    experimentData = expdata,
    spectraProcessingQueue = new_procQ)
  if (validObject(res))
    res
}
