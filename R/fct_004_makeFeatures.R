

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
#' @param rtAlignment Logical, set to \code{TRUE} (the default) to preform
#' retention time aligment besides peak grouping.
#' @param paramAlignment Applicable for algorithm "xcms3" only,
#' the parameters for the chosen retention time alignment method.
#' See documentation of \code{\link[xcms]{adjustRtime}} for more information.
#' For alignment with the method \code{PeakGroups}, a list of length two
#' should be given, where the first element
#' is the pre-grouping parameters using \code{\link[xcms]{groupChromPeaks}}
#' and the second element the actual alignment parameters.
#' @param paramGrouping The parameters for the chosen grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} or
#' \code{\link[patRoon]{groupFeatures}} for more information.
#' @param recurvive Logical, set to \code{TRUE} for applying recursive
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
#' @importClassesFrom patRoon featureGroups
#' @importClassesFrom xcms XCMSnExp
#' @importFrom patRoon groupFeatures getXCMSnExp importFeaturesXCMS3 as.data.table as.data.frame groupFeatIndex
#' @importMethodsFrom xcms groupChromPeaks adjustRtime fillChromPeaks featureDefinitions
#' @importFrom dplyr count select rename everything
#' @importFrom data.table rbindlist
#'
#' @examples
#'
makeFeatures <- function(obj = NULL,
                         algorithm = "xcms3",
                         rtAlignment = TRUE,
                         paramAlignment = NULL,
                         paramGrouping = NULL,
                         recursive = TRUE,
                         paramFill = NULL,
                         save = TRUE) {

  x <- obj@patdata
  
  if (is.null(paramAlignment)) paramAlignment <- obj@parameters$peakAlignment
  
  if (is.null(paramGrouping)) paramGrouping <- obj@parameters$peakGrouping
  
  if (is.null(paramFill)) paramFill <- obj@parameters$fillMissing
  
  if (algorithm == "xcms3") {
    x <- getXCMSnExp(x)
    if (rtAlignment) {
      if (length(x$sample_name) > 1) {
        if (class(paramAlignment[[2]]) == "PeakGroupsParam") {
          paramAlignment[[1]]@sampleGroups <- x$sample_group
          x <- groupChromPeaks(x, param = paramAlignment[[1]])
          x <- adjustRtime(x, msLevel = 1, param = paramAlignment[[2]])
        } else {
          x <- adjustRtime(x, msLevel = 1, param = paramAlignment)
        }
      }
    }
    paramGrouping@sampleGroups <- x$sample_group
    x <- groupChromPeaks(x, param = paramGrouping)

  } else {
    ag <- list(feat = x, rtalign = rtAlignment, algorithm = algorithm)
    x <- do.call(groupFeatures, c(ag, paramGrouping, verbose = TRUE))
  }

  if (recursive) {
    if (class(x) != "XCMSnExp") {
      x <- getXCMSnExp(x)
    }
    x <- recursiveIntegration(XCMSfeatures = x, paramFill = paramFill)
  }

  obj@parameters$peakAlignment <- paramAlignment
  
  obj@parameters$peakGrouping <- paramGrouping
  
  obj@parameters$fillMissing <- paramFill
  
  if (class(x) == "XCMSnExp") {
    sinfo <- data.frame(path = dirname(fileNames(x)),
                        analysis = x$sample_name,
                        group = x$sample_group,
                        blank = rep(NA, length(x$sample_name)))

    x <- importFeatureGroupsXCMS3(x, sinfo)
  }

  obj@patdata <- x
  
  #updates slot peaks with filled peaks, possibly added to patRoon new version
  if (class(x) == "featureGroupsXCMS3") {
    peaks <- as.data.frame(chromPeaks(x@xdata, isFilledColumn = TRUE))
    peaks$group <- sapply(peaks$sample, FUN = function(x) x <- obj@samples$group[x])
    peaks$sample <- sapply(peaks$sample, FUN = function(x) x <- obj@samples$sample[x])
    peaks <- split(peaks, f = peaks$sample)
    peaks <- as.data.frame(rbindlist(peaks))
    peaks$ID <- 1:nrow(peaks)
    peaks <- select(peaks, ID, sample, group, everything())
    peaks <- rename(peaks, intensity = maxo, area = into)
    obj@peaks <- peaks
  }
  
  #add maker for feature list
  obj <- buildFeatureList2(obj)
  
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
#' @export
#'
recursiveIntegration <- function(XCMSfeatures = XCMSfeatures,
                                 paramFill = xcms::ChromPeakAreaParam()) {

  XCMSfeatures <- suppressWarnings(fillChromPeaks(XCMSfeatures, param = paramFill))

  return(XCMSfeatures)

}




buildFeatureList2 <- function(obj) {
  
  x <- obj@patdata
  
  feat <- patRoon::as.data.frame(x, average = TRUE)
  
  feat <- rename(feat, rt = ret, ID = group)
  
  rgs <- obj@samples$group
  
  rg <- unique(rgs)
  
  sp <- obj@samples$sample
  
  names(rgs) <- sp
  
  for (i in seq_len(length(rg))) {
    feat[, paste0(rg[i], "_sd%")] <- apply(
      X = patRoon::as.data.table(x, average = FALSE)[, .SD, .SDcols = sp[rgs == rg[i]]],
      MARGIN = 1, function(x) {
        round(ifelse(sd(x) != 0, (sd(x)/mean(x))*100, NA), digits = 0)})
  }
  
  pl <- obj@peaks
  
  fidx <- patRoon::groupFeatIndex(x)
  fidx[fidx == 0] <- NA 
  total <- nrow(pl[pl$sample %in% sp[1],])
  for (i in 2:length(sp)) {
    newtotal <- nrow(pl[pl$sample %in% sp[i],])
    fidx[i, ] <- fidx[i, ] + total
    total <- total + newtotal
  }
  fidx <-  as.list(fidx)
  
  feat <- feat[!sapply(feat[, rg] == 0, function(x) sum(x) == length(rg)), ]
  
  fidx <- fidx[names(fidx) %in% feat$ID]
  
  feat$mzmin <- unlist(lapply(fidx, function(x) {min(pl[na.omit(x), "mzmin", drop = T])}))
  feat$mzmax <- unlist(lapply(fidx, function(x) {max(pl[na.omit(x), "mzmax", drop = T])}))
  feat$rtmin <- unlist(lapply(fidx, function(x) {min(pl[na.omit(x), "rtmin", drop = T])}))
  feat$rtmax <- unlist(lapply(fidx, function(x) {max(pl[na.omit(x), "rtmax", drop = T])}))
  
  feat$dppm <- round(((feat$mzmax - feat$mzmin) / feat$mz) * 1E6, digits = 1)
  
  feat$width <- round(feat$rtmax - feat$rtmin, digits = 0)
  
  feat$pIdx <- I(fidx)
  
  fidx2 <- fidx
  fidx2 <- do.call(rbind, fidx2)
  colnames(fidx2) <- rgs
  fidx2[is.na(fidx2)] <- 0
  fidx2[fidx2 > 0] <- 1
  fidx2 <- sapply(rg, function(x) rowSums(fidx2[, grepl(x, rgs)]))
  feat$npeaks <- I(lapply(rownames(fidx2), function(x) {fidx2[x, , drop = TRUE]}))
  
  feat$hasFilled <- unlist(lapply(fidx, function(x) {1 %in% pl[na.omit(x), "is_filled"]}))
  
  feat <- select(feat, ID, mz, rt, dppm, width, everything())
  
  obj@features <- feat
  
  return(obj)
  
}
