

#' @title makeFeatures_Old
#'
#' @description Grouping and alignment of peaks across samples.
#' The peak grouping uses the \pkg{patRoon} package and
#' the following algorithms are possible: "xcms3", "xcms", "openms".
#' See ?\pkg{patRoon} for further information.
#' When fill missing parameters (paramFill) are defined
#' a recursive integration for samples with missing peaks
#' can be applied through the \pkg{xcms} package.
#'
#' @param obj A \linkS4class{ntsData} object containing peaks.
#' @param algorithm The algorithm to use for peak alignment and grouping.
#' @param paramGrouping The parameters for the chosen grouping method.
#' See \code{\link[patRoon]{groupFeatures}} for more information.
#' @param paramFill A list of parameters for recursive integration of missing peaks.
#' See \code{\link[xcms]{fillChromPeaks}} for more information.
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
makeFeatures_Old <- function(obj = NULL,
                         algorithm = NULL,
                         paramGrouping = NULL,
                         paramFill = NULL,
                         save = TRUE) {

  assertClass(obj, "ntsData")

  x <- obj@patdata

  if (is.null(algorithm)) algorithm <- peakGroupingParameters(obj)@algorithm
  
  if (is.null(paramGrouping)) paramGrouping <- peakGroupingParameters(obj)@param
  
  if (is.null(paramFill)) paramFill <- fillMissingParameters(obj)@param
  
  if (is.na(algorithm)) {
    warning("Peak grouping algorihtm not defined!")
    return(obj)
  }
  
  if (algorithm == "xcms3") {
      paramGrouping$groupParam@sampleGroups <- x@analysisInfo$group
    if (paramGrouping$rtalign) {
      paramGrouping$preGroupParam@sampleGroups <- x@analysisInfo$group
    }
  }
  
  ag <- list(obj = x, algorithm = algorithm)
  
  x <- do.call(groupFeatures, c(ag, paramGrouping, verbose = TRUE))
  
  
  if (length(paramFill) > 0) {
    
    if (class(x) != "XCMSnExp") x <- getXCMSnExp(x)
    
    if (is.list(paramFill)) paramFill <- paramFill[[1]]
    
    x <- recursiveIntegration(XCMSfeatures = x, paramFill = paramFill)

    sinfo <- data.frame(path = dirname(fileNames(x)),
                        analysis = x$sample_name,
                        group = x$sample_group,
                        blank = obj@patdata@analysisInfo$blank)
    
    sinfo$blank[is.na(sinfo$blank)] <- ""
    
    x <- importFeatureGroupsXCMS3(x, sinfo)

  }
  
  obj <- peakGroupingParameters(obj, algorithm = algorithm, param = paramGrouping)
  
  obj <- fillMissingParameters(obj, algorithm = "xcms3", param = paramFill)

  obj@patdata <- x

  obj <- buildPeakList(obj)
  
  obj <- buildFeatureList(obj)

  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}




#' recursiveIntegration
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
recursiveIntegration <- function(XCMSfeatures = XCMSfeatures,
                                 paramFill = xcms::ChromPeakAreaParam()) {

  XCMSfeatures <- suppressWarnings(xcms::fillChromPeaks(XCMSfeatures, param = paramFill))

  return(XCMSfeatures)

}




#' buildFeatureList
#'
#' @param obj An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with updated features slot.
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

  feat$mzmin <- unlist(lapply(index, function(h) min(pl[na.omit(h), "mzmin", drop = TRUE])))
  feat$mzmax <- unlist(lapply(index, function(h) max(pl[na.omit(h), "mzmax", drop = TRUE])))
  feat$rtmin <- unlist(lapply(index, function(h) min(pl[na.omit(h), "rtmin", drop = TRUE])))
  feat$rtmax <- unlist(lapply(index, function(h) max(pl[na.omit(h), "rtmax", drop = TRUE])))

  feat$dppm <- round(((feat$mzmax - feat$mzmin) / feat$mz) * 1E6, digits = 1)

  feat$width <- round(feat$rtmax - feat$rtmin, digits = 0)

  feat$pIdx <- I(index)

  pl$group <- factor(pl$group, levels = rg)

  feat$npeaks <- I(lapply(index, function(h) {
    as.data.frame(table(pl$group[na.omit(h)]))$Freq
  }))

  feat$hasFilled <- unlist(lapply(index, function(x) {
    1 %in% pl[na.omit(x), "is_filled"]
  }))

  feat <- select(feat, ID, mz, rt, dppm, width, everything())

  obj@features <- feat

  return(obj)

}




#' updateFeatureList
#'
#' @param obj An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with updated features slot.
#'
#' @importMethodsFrom patRoon as.data.frame
#' @importFrom dplyr select everything
#' @importFrom stats na.omit
#'
updateFeatureList <- function(obj) {

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

  # TODO Issue, takes longer than buildFeatureList because the indices are not ordered
  feat$mzmin <- unlist(lapply(index, function(h) min(pl$mzmin[pl$ID %in% na.omit(h)])))
  feat$mzmax <- unlist(lapply(index, function(h) max(pl$mzmax[pl$ID %in% na.omit(h)])))
  feat$rtmin <- unlist(lapply(index, function(h) min(pl$rtmin[pl$ID %in% na.omit(h)])))
  feat$rtmax <- unlist(lapply(index, function(h) max(pl$rtmax[pl$ID %in% na.omit(h)])))

  feat$dppm <- round(((feat$mzmax - feat$mzmin) / feat$mz) * 1E6, digits = 1)

  feat$width <- round(feat$rtmax - feat$rtmin, digits = 0)

  feat$pIdx <- I(index)

  pl$group <- factor(pl$group, levels = rg)

  feat$npeaks <- I(lapply(index, function(h) {as.data.frame(table(pl$group[na.omit(h)]))$Freq}))

  feat$hasFilled <- unlist(lapply(index, function(x) {
    1 %in% pl[na.omit(x), "is_filled"]
  }))

  feat <- select(feat, ID, mz, rt, dppm, width, everything())

  featOld <- obj@features

  featOld <- featOld[, c("ID", colnames(featOld)[!colnames(featOld) %in% colnames(feat)])]

  feat <- left_join(feat, featOld, by = "ID")

  obj@features <- feat

  return(obj)

}
