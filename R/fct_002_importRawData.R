


#' @title importRawData
#' @description Function to import MS data from the mzML files listed in \code{sampleInfo} from the \code{makeSetup} object.
#' Files should contain centroided spectra.
#' The function uses the function \code{\link[MSnbase]{readMSData}} from the \code{MSnbase} package to read the mzML files.
#' 
#' @param sampleInfo The sample \code{data.table} obtained by the \code{\link{makeSetup}} function.
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
#'
#' @examples
#' 
#' 
#' 
importRawData <- function(sampleInfo = setup$sampleInfo, rtFilter = NULL, timeUnit = "min",
                          msLevel = c(1,2), centroidedData = TRUE, removeEmptySpectra = TRUE, save = TRUE, projPath = setup$projPath) {
  
  #Examples
  # setup <- ntsIUTA::makeSetup(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)
  # rawDataExample <- ntsIUTA::importRawData(setup$sampleInfo[1:3,], save = FALSE, centroidedData = TRUE)
  # rawDataExample
  
  rawData <- base::suppressWarnings(MSnbase::readMSData(sampleInfo$filePath[drop = TRUE],
                                                        pdata = methods::new("NAnnotatedDataFrame",
                                                                             base::data.frame(sample_name = sampleInfo$sample,
                                                                                              sample_group = sampleInfo$group)),
                                                        mode = "onDisk",
                                                        centroided. = centroidedData,
                                                        smoothed. = FALSE))
  
  if(!base::is.null(rtFilter))
  {
    if(timeUnit == "min") rtFilter = rtFilter*60
    rawData <- MSnbase::filterRt(rawData, rt = rtFilter, msLevel. = msLevel)
  }
  
  if(removeEmptySpectra) rawData <- base::suppressWarnings(MSnbase::filterEmptySpectra(rawData))
  
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(rawData, file = base::paste0(rData,"\\rawData.rds"))
  }
  
  return(rawData)
  
  
  #Options to save information from raw data
  #raw[["rawchrom"]] <- chromatogram(rawData, aggregationFun = "max")
  #raw[["emptyspec"]] <- length(which(MSnbase::peaksCount(rawData) == 0))
  
  #TODO For examples afterwards
  # file <- dir(system.file(package = "MSnbase", dir = "extdata"),
  #             full.name = TRUE,
  #             pattern = "mzXML$")
  
  #TODO Para verificar se of ficheiros tem centroided MS data ou nao. Pode sugerir para fazer centroiding
  #TestDataTypeCent <- table(MSnbase::isCentroided(MSnbase::filterFile(rawData_cent, 1)), useNA = "always")
  #TestDataTypeCent <- table(MSnbase::isCentroidedFromFile(MSnbase::filterFile(rawData_cent, 1)), useNA = "always")
  
  
}



#TODO To add better descriptions

#' @title centroidProfileData
#' @description Centroid profile MS data with additional smoothing and refinement of the m/z.
#' 
#' @param rawData The \linkS4class{OnDiskMSnExp} object with profile data for centroiding.
#' @param smoothing Logical, set to \code{TRUE} for applying smothing to the profile data.
#' @param refineMZ Logical, set to \code{TRUE} for applying \emph{m/z} refinement before centroiding.
#' @param methodSmoothing Method for smoothing. 
#' @param halfWindowSize The window size of the filter.
#' @param polynomialOrder Polynomial order for smothing, default is 4.
#' @param methodRefineMz Method for refinement.
#' @param k Number of closest signals to the centroid.
#' @param signalPercentage Minimum signal percentage of centroids to refine \emph{m/z}.
#' @param stopAtTwo Logical, set to \code{TRUE} for allowed two consecutive equal or higher signals.
#' \code{FALSE} will stop when one equal or higher centroid is found.
#' @param save Logical, set to \code{TRUE} to save the object rawData in the \strong{rData} folder.
#' Note, that the profile files will be replaced by the produced centroided files.
#' @param projPath If save is \code{TRUE}, the saving location or project path should be given.
#' The default is the \code{projPath} in the object \code{setup}.
#'
#' @return Centroiding of profile data from mzML and mzXML files in an \linkS4class{OnDiskMSnExp} object.
#' 
#' @export
#'
#' @import magrittr
#' @importFrom MSnbase smooth pickPeaks fileNames writeMSData
#'
#' @examples
#' 
#' 
#' 
centroidProfileData <- function(rawData, smoothing = TRUE, refineMZ = TRUE,
                                methodSmoothing = "SavitzkyGolay", halfWindowSize = 5, polynomialOrder = 4,
                                methodRefineMz = "kNeighbors", k = 2, signalPercentage = 10, stopAtTwo = TRUE,
                                save = TRUE, projPath = setup$projPath) {
  
  #require(magrittr)
  
  if (smoothing)
  {
    rawData <- rawData %>% MSnbase::smooth(method = methodSmoothing,
                                           halfWindowSize = halfWindowSize,
                                           polynomialOrder = polynomialOrder)
  }
  
  if (refineMZ)
  {
    if (methodRefineMz == "kNeighbors") {
      rawData <- MSnbase::pickPeaks(rawData, refineMz = methodRefineMz, k = k) }
    
    if (methodRefineMz == "descendPeak") {
      rawData <- MSnbase::pickPeaks(rawData, refineMz = methodRefineMz, signalPercentage = 0.1, stopAtTwo = TRUE) }
  }
  
  #Save centroided files in disk
  if (save) {
    fls_new <- MSnbase::fileNames(rawData)
    MSnbase::writeMSData(data_cent, file = fls_new)
    
  }
  
  #TODO Finish function to centroid data
  return(rawData)
  
}
