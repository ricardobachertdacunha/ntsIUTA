

#' @title peakFilling
#'
#' @description Recursive integration for samples with missing peaks
#' using the function \code{\link[xcms]{fillChromPeaks}}
#' from the \pkg{xcms} package.
#'
#' @param XCMSfeatures An \linkS4class{XCMSnExp} object containing features.
#' @param paramFill An object of class \code{FillChromPeaksParam} or \code{ChromPeakAreaParam}
#' containing the parameters to apply the recursive integration.
#' See \code{?\link[xcms]{fillChromPeaks}} for more information.
#'
#' @return An \linkS4class{XCMSnExp} object including filled missing peaks.
#'
#' @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#'
#' @importMethodsFrom xcms fillChromPeaks
#' @importClassesFrom xcms ChromPeakAreaParam FillChromPeaksParam
#'
peakFilling <- function(object, algorithm = NULL, settings = NULL, save = TRUE) {


  x <- getXCMSnExp(x)


  paramFill <- xcms::ChromPeakAreaParam()



  XCMSfeatures <- suppressWarnings(xcms::fillChromPeaks(XCMSfeatures, param = paramFill))

  return(XCMSfeatures)

 if (length(paramFill) > 0) {
    
    if (class(x) != "XCMSnExp")
    
    if (is.list(paramFill)) paramFill <- paramFill[[1]]
    
    x <- recursiveIntegration(XCMSfeatures = x, paramFill = paramFill)

    sinfo <- data.frame(path = dirname(fileNames(x)),
                        analysis = x$sample_name,
                        group = x$sample_group,
                        blank = obj@patdata@analysisInfo$blank)
    
    sinfo$blank[is.na(sinfo$blank)] <- ""
    
    x <- importFeatureGroupsXCMS3(x, sinfo)

  }











}
