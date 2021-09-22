

#' @title makeFeatures
#' @description The \code{makeFeatures} consists of alignment and grouping of
#' chromatographic peaks across samples. The function uses algorithms
#' from the package \pkg{xcms} or from \pkg{patRoon} based on "openms".
#' Additionally, a recursive integration for samples with missing peaks
#' can be applied using the function \code{\link[xcms]{fillChromPeaks}}
#' from the \pkg{xcms} package.
#'
#' @param obj A \linkS4class{ntsData} object containing peaks
#' for one or more samples.
#' @param algorithm The algorithm to use for peak alignment and grouping.
#' One of "xcms3" (the default) or "openms".
#' When "openms" the \pkg{patRoon} package is used.
#' @param paramGrouping The parameters for the chosen grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} or
#' \code{\link[patRoon]{groupFeatures}} for more information.
#' @param recursive Logical, set to \code{TRUE} for applying recursive
#' integration for samples with missing peaks.
#' @param paramFill An object of class \code{FillChromPeaksParam}
#' or \code{ChromPeakAreaParam} containing the parameters
#' to apply the recursive integration.
#' #' See \code{?\link[xcms]{fillChromPeaks}} for more information.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object
#' containing peaks aligned and grouped as features across samples accessible
#' in the slot \code{features}.
#'
#'  @references
#' \insertRef{Helmus2021}{ntsIUTA}
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#'
#' @export
#'
#' @importFrom checkmate assertClass testChoice
#' @importClassesFrom patRoon featureGroups
#' @importClassesFrom xcms XCMSnExp
#' @importFrom patRoon groupFeatures getXCMSnExp importFeaturesXCMS3 groupFeatIndex
#' @importMethodsFrom patRoon as.data.table as.data.frame
#' @importMethodsFrom xcms groupChromPeaks adjustRtime fillChromPeaks featureDefinitions
#' @importFrom dplyr select rename everything
#' @importFrom data.table rbindlist
#'
makeFeatures <- function(obj = NULL,
                         algorithm = NULL,
                         paramGrouping = NULL,
                         recursive = TRUE,
                         paramFill = NULL,
                         save = TRUE) {

  assertClass(obj, "ntsData")

  x <- obj@patdata

  if (is.null(algorithm)) algorithm <- obj@algorithms$makeFeatures

  if (!testChoice(algorithm, c("xcms3", "xcms", "openms"))) {
    warning("Algorithm not recognized. See ?makeFeatures for more information.")
    return(obj)
  }

  if (is.null(paramGrouping)) paramGrouping <- obj@parameters$peakGrouping

  if (is.null(paramFill)) paramFill <- obj@parameters$fillMissing

  if (algorithm == "xcms3") {
      paramGrouping$groupParam@sampleGroups <- x@analysisInfo$group
    if (paramGrouping$rtalign) {
      paramGrouping$preGroupParam@sampleGroups <- x@analysisInfo$group
    }
  }
  
  ag <- list(obj = x, algorithm = algorithm)
  
  x <- do.call(groupFeatures, c(ag, paramGrouping, verbose = TRUE))
  
  
  if (recursive) {
    
    if (class(x) != "XCMSnExp") x <- getXCMSnExp(x)
    
    x <- recursiveIntegration(XCMSfeatures = x, paramFill = paramFill)

    sinfo <- data.frame(path = dirname(fileNames(x)),
                        analysis = x$sample_name,
                        group = x$sample_group,
                        blank = obj@patdata@analysisInfo$blank)

    x <- importFeatureGroupsXCMS3(x, sinfo)

  }


  obj@algorithms$makeFeatures <- algorithm

  obj@parameters$peakGrouping <- paramGrouping

  obj@parameters$fillMissing <- paramFill

  obj@patdata <- x

  obj <- buildPeakList(obj)
  
  obj <- buildFeatureList(obj)

  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}



#' @title recursiveIntegration
#' @description Recursive integration for samples with missing peaks can be applied
#' using the function \code{\link[xcms]{fillChromPeaks}} from the \pkg{xcms} package.
#'
#' @param XCMSfeatures An \linkS4class{XCMSnExp} object containing features.
#' @param paramFill An object of class \code{FillChromPeaksParam} or \code{ChromPeakAreaParam}
#' containing the parameters to apply the recursive integration.
#' #' See \code{?\link[xcms]{fillChromPeaks}} for more information.
#'
#' @return An \linkS4class{XCMSnExp} object including filled missing peaks.
#'
#' @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#'
recursiveIntegration <- function(XCMSfeatures = XCMSfeatures,
                                 paramFill = xcms::ChromPeakAreaParam()) {

  XCMSfeatures <- suppressWarnings(fillChromPeaks(XCMSfeatures, param = paramFill))

  return(XCMSfeatures)

}


#' buildFeatureList
#'
#' @param obj An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with updated features slot.
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.frame
#' @importFrom dplyr select everything
#' @importFrom stats na.omit
#'
buildFeatureList <- function(obj) {
  
  pat <- obj@patdata
  
  feat <- patRoon::as.data.frame(pat, average = TRUE)
  
  feat <- rename(feat, rt = ret, ID = group)
  
  rgs <- obj@samples$group
  
  rg <- unique(rgs)
  
  sp <- obj@samples$sample
  
  names(rgs) <- sp
  
  for (i in seq_len(length(rg))) {
    feat[, paste0(rg[i], "_sd%")] <- apply(
      X = patRoon::as.data.table(pat, average = FALSE)[, .SD, .SDcols = sp[rgs == rg[i]]],
      MARGIN = 1, function(x) {
        round(ifelse(sd(x) != 0, (sd(x) / mean(x)) * 100, NA), digits = 0)})
  }
  
  pl <- obj@peaks
  
  index <- lapply(seq_len(nrow(feat)), function(h) pl$ID[pl$feature %in% feat$ID[h]])
  
  names(index) <- feat$ID
  
  feat$mzmin <- unlist(lapply(index, function(h) {min(pl[na.omit(h), "mzmin", drop = TRUE])}))
  feat$mzmax <- unlist(lapply(index, function(h) {max(pl[na.omit(h), "mzmax", drop = TRUE])}))
  feat$rtmin <- unlist(lapply(index, function(h) {min(pl[na.omit(h), "rtmin", drop = TRUE])}))
  feat$rtmax <- unlist(lapply(index, function(h) {max(pl[na.omit(h), "rtmax", drop = TRUE])}))
  
  feat$dppm <- round(((feat$mzmax - feat$mzmin) / feat$mz) * 1E6, digits = 1)
  
  feat$width <- round(feat$rtmax - feat$rtmin, digits = 0)
  
  feat$pIdx <- I(index)
  
  pl$group <- factor(pl$group, levels = rg)
  
  feat$npeaks <- I(lapply(index, function(h) {as.data.frame(table(pl$group[na.omit(h)]))$Freq}))
  
  feat$hasFilled <- unlist(lapply(index, function(x) {
    1 %in% pl[na.omit(x), "is_filled"]
  }))
  
  feat <- select(feat, ID, mz, rt, dppm, width, everything())
  
  obj@features <- feat
  
  return(obj)
  
}
