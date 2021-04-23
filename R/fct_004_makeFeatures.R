



#' @title makeFeatures
#' @description The \code{makeFeatures} consists of alignment and grouping of chromatographic peaks across samples.
#' The function uses methods from the package \pkg{xcms}, specifically \code{\link[xcms]{adjustRtime}} and \code{\link[xcms]{groupChromPeaks}}
#' methods. The input is a \linkS4class{XCMSnExp} object or a \code{list} of \linkS4class{XCMSnExp} objects
#' named with the given sample replicate group during the \code{setup}. If the input is a \code{list}, \linkS4class{XCMSnExp} objects
#' are concatenated to produce a merged \linkS4class{XCMSnExp} object,
#' containing features (\emph{i.e.} grouped chromatographic peaks across samples).
#'
#' @param peaksData A \linkS4class{XCMSnExp} object or a \code{list} of \linkS4class{XCMSnExp} objects named according to the 
#' sample replicate group represented by each \linkS4class{XCMSnExp} object as obtained by \code{\link{peakPicking}}.
#' @param paramPreGrouping If \code{\link[xcms]{adjustRtime}} is preformed with the method \code{PeakGroups},
#' a pre-grouping of peaks is required for alignment. The \code{paramPreGrouping} is the parameters obtained by the selected grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} for more information.
#' @param paramAlignment The parameters for the chosen alignment method. See documentation of \code{\link[xcms]{adjustRtime}} for more information.
#' @param paramGrouping The parameters for the chosen grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} for more information.
#' @param ppmForFillingGroups The mass (in ppm) to expand the \emph{mz} for filling missing peaks in incomplete features.
#' See \code{\link[xcms]{fillChromPeaks}} for more information.
#' @param save Logical, set to \code{TRUE} to save the generated \code{XCMSnExp} object in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#' @param maxMultiProcess Logical, set to \code{TRUE} to enable max parallel processing. Changes the number of workers to the maximum available.
#'
#' @return An \linkS4class{XCMSnExp} containing features which are grouped and aligned peaks across samples.
#' 
#'  @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#' 
#' @importFrom BiocParallel registered SnowParam register bpparam
#' @importFrom parallel detectCores
#' @importFrom xcms sampleGroups groupChromPeaks adjustRtime fillChromPeaks ChromPeakAreaParam
#'
#' @examples
#' 
#' 
#' 
makeFeatures <- function(peaksData = peaksData,
                         paramPreGrouping = instParam$preGrouping,
                         paramAlignment = instParam$alignment,
                         paramGrouping = instParam$grouping,
                         ppmForFillingGroups = 5,
                         save = TRUE, projPath = setup$projPath, maxMultiProcess = TRUE) {
  
  #Examples
  #setup <- ntsIUTA::makeSetup(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)
  # setup$sampleInfo[1:3,"group"] <- "Sample"
  # rawDataExample <- ntsIUTA::importRawData(setup$sampleInfo[1:3,], save = FALSE, centroidedData = TRUE)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = param, save = FALSE)
  # paramList <- paramListExample
  # instParam <- getInstParam(paramList)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = instParam$PP, save = FALSE)
  # featDataExample <- ntsIUTA::makeFeatures(peaksDataExample, paramPreGrouping = instParam$preGrouping, paramAlignment = instParam$alignment, paramGrouping = instParam$grouping, save = FALSE)
  # featDataExample
  
  
  
  # require(xcms)
  # require(BiocParallel)
  # require(parallel)
  # require(grDevices)
  # require(RColorBrewer)
  # require(dplyr)
  # require(graphics)
  
  #Enable full parallel processing
  if (maxMultiProcess)
  {
    snow <- BiocParallel::registered("SnowParam")
    if (snow$workers < parallel::detectCores())
    {
      snow <- BiocParallel::SnowParam(workers = parallel::detectCores(), type = "SOCK", exportglobals = FALSE)
      BiocParallel::register(snow, default = TRUE)
    }
  }
  
  
  #Concatenation for peaksData objects before grouping
  if (base::class(peaksData) == "XCMSnExp")
  {
    peaksDataUnified <- peaksData
  } else {
    
    if (base::length(peaksData) == 1)
    {peaksDataUnified <- peaksData[[1]]}
    else {peaksDataUnified <- base::do.call(c, base::unlist(peaksData, recursive = FALSE))}
    
  }
  
  if (base::length(peaksDataUnified$sample_name) > 1) {
    
    #Pregrouping necessary for alignment with PeakGroups method from xcms
    if(base::class(paramAlignment) == "PeakGroupsParam")
    {
      xcms::sampleGroups(paramPreGrouping) <- peaksDataUnified$sample_group
      featData <- xcms::groupChromPeaks(peaksDataUnified, param = paramPreGrouping)
    }
    
    
    featData <- xcms::adjustRtime(featData, msLevel = 1, param = paramAlignment)
    
    xcms::sampleGroups(paramGrouping) <- peaksDataUnified$sample_group
    featData <- xcms::groupChromPeaks(featData, param = paramGrouping)
    
  } else {
    
    xcms::sampleGroups(paramGrouping) <- peaksDataUnified$sample_group
    featData <- xcms::groupChromPeaks(peaksDataUnified, param = paramGrouping)
    
  }
  
  #Fill missing peaks in feature groups
  #Old parameter used: FillChromPeaksParam(ppm = ppmForFillingGroups)
  featData <- base::suppressWarnings(xcms::fillChromPeaks(featData, param = xcms::ChromPeakAreaParam(),
                                                          BPPARAM = BiocParallel::bpparam("SnowParam")))
  
  #TODO add optinally to apply correction of retention time, avoiding another function from xcms
  # if (length(peaksDataUnified$sample_name) > 1) {
  #   featData <- xcms::applyAdjustedRtime(featData)
  # }
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(featData, file = base::paste0(rData,"\\featData.rds"))
  }
  
  return(featData)
}



