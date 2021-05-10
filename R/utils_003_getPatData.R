


#' @title getPatData
#'
#' @param featData A \linkS4class{XCMSnExp} object with feature data to be converted.
#' @param sampleInfo The \code{sampleInfo} obtained by the \code{\link{setupProject}} function, matching the \code{featData}.
#' @param save Logical, set to \code{TRUE} to save the generated \linkS4class{featureGroups} object in the disk.
#' @param projPath The \code{projPath} in the \code{setup} object, or the location where to save the generated R object.
#'
#' @return A \linkS4class{featureGroups} object corresponding matching the given \linkS4class{XCMSnExp} object.
#' 
#' @importFrom MSnbase fileNames
#' @importFrom patRoon importFeatureGroupsXCMS3
#' 
#' @export
#'
#' @examples
#' 
#' 
#' 
getPatData <- function(featData, sampleInfo, save = FALSE, projPath = setup$projPath) {
  
  #Examples
  # patDataExample <- getPatData(ntsIUTA::featDataExample, ntsIUTA::setupExample$sampleInfo)
  # patDataExample
  
  patSampleInfo <- sampleInfo[sampleInfo$sample %in% featData$sample_name,]
  patSampleInfo$filePath <- base::dirname(sampleInfo$filePath[sampleInfo$sample %in% featData$sample_name])
  patSampleInfo <- base::as.data.frame(patSampleInfo[,1:4])
  base::colnames(patSampleInfo) <- c("path", "analysis", "group", "blank")
  
  if (!base::all.equal(patSampleInfo$path, base::dirname(MSnbase::fileNames(featData))))
  {
    stop("File locations do not match between featData and sampleInfo.")
  }
  
  patData <- patRoon::importFeatureGroupsXCMS3(featData, patSampleInfo)
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(patData, file = base::paste0(rData,"\\patData.rds"))
  }
  
  return(patData)
  
}
