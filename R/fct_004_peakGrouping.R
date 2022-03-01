

#' @title peakGrouping
#'
#' @description Grouping and alignment of peaks across samples.
#' The peak grouping uses the \pkg{patRoon} package and
#' the following algorithms are possible: "xcms3", "xcms", "openms".
#' See ?\pkg{patRoon} for further information.
#'
#' @param object A \linkS4class{ntsData} object containing peaks.
#' @param algorithm The algorithm to use for peak alignment and grouping.
#' @param settings The respective parameter settings.
#' See \link[patRoon]{groupFeatures} for more information.
#' @param simplify Logical, set to \code{TRUE} to convert the patRoon object
#' to a \linkS4class{featureGroupsSIRIUS} object which is simpler and faster for subsetting.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object
#' containing peaks aligned and grouped as features across samples.
#'
#'  @references
#' \insertRef{Helmus2021}{ntsIUTA}
#'
#' @export
#'
#' @importFrom checkmate assertClass testClass
#' @importClassesFrom patRoon featureGroups
#' @importFrom patRoon groupFeatures
#' @importMethodsFrom patRoon as.data.table
#'
peakGrouping <- function(object = NULL,
                         algorithm = NULL,
                         settings = NULL,
                         simplify = TRUE,
                         save = TRUE) {

  assertClass(object, "ntsData")

  x <- object@pat

  if (checkmate::testClass(x, "featureGroups")) x <- x@features

  if (is.null(algorithm)) algorithm <- groupingParameters(object)@algorithm

  if (is.null(settings)) settings <- groupingParameters(object)@settings

  if (is.na(algorithm)) {
    warning("Peak grouping algorihtm not defined!")
    return(object)
  }

  if (algorithm == "xcms3") {
      settings$groupParam@sampleGroups <- x@analysisInfo$group
    if (settings$rtalign) {
      settings$preGroupParam@sampleGroups <- x@analysisInfo$group
    }
  }

  ag <- list(obj = x, algorithm = algorithm)

  x <- do.call(groupFeatures, c(ag, settings, verbose = TRUE))

  object <- groupingParameters(object, algorithm = algorithm, settings = settings)

  object@pat <- x

  object <- buildPeaksTable(object)

  object <- buildFeatureTable(object)

  if (simplify) {
    xFeats <- new("featuresSIRIUS",
      analysisInfo = x@features@analysisInfo,
      features = x@features@features
    )

    xFeatGroups <- new("featureGroupsSIRIUS",
      groups = x@groups,
      groupInfo = x@groupInfo,
      analysisInfo = x@analysisInfo,
      features = xFeats,
      ftindex = x@ftindex
    )

    object@pat <- xFeatGroups
  }

  if (save) saveObject(object = object)

  if (is.character(save)) saveObject(object = object, filename = save)

  return(object)
}




#' buildFeatureTable
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with updated features slot.
#'
#' @importMethodsFrom patRoon as.data.table
#' @importFrom data.table setnames setorder
#' @importFrom dplyr select everything
#'
buildFeatureTable <- function(object) {

  cat("Building features table... \n")

  pat <- object@pat

  feat <- patRoon::as.data.table(pat, average = TRUE)

  feat <- setnames(feat, c("ret", "group"), c("rt", "id"))

  feat_b <- patRoon::as.data.table(pat, average = FALSE)

  rpl <- unique(replicates(object))

  rpl_sp <- lapply(rpl, function(x, st) {
    st$sample[st$replicate == x]
  }, st = samplesTable(object))

  feat_sd <- lapply(rpl_sp, function(x, feat_b) {
    temp <- feat_b[, x, with = FALSE]
    temp <- apply(temp, 1, function(x) sd(x) / mean(x) * 100)
    temp[is.nan(temp)] <- 0
    temp <- round(temp, digits = 0)
    return(temp)
  }, feat_b = feat_b)

  names(feat_sd) <- paste0(rpl, "_sd")
  feat <- cbind(feat, as.data.table(feat_sd))

  pk <- peaks(object)

  pk$index <- seq_len(nrow(pk))
  index <- lapply(feat$id, function(x) pk[feature == x, index])
  feat$mzmin <- unlist(lapply(index, function(x) min(pk[x, mzmin])))
  feat$mzmax <- unlist(lapply(index, function(x) max(pk[x, mzmax])))
  feat$rtmin <- unlist(lapply(index, function(x) min(pk[x, rtmin])))
  feat$rtmax <- unlist(lapply(index, function(x) max(pk[x, rtmax])))

  feat$d_ppm <- round(((feat$mzmax - feat$mzmin) / feat$mz) * 1E6, digits = 1)

  feat$d_sec <- round(feat$rtmax - feat$rtmin, digits = 0)

  feat$p_id <- I(lapply(feat$id, function(x) pk[feature == x, id]))

  pk$replicate <- factor(pk$replicate, levels = rpl)

  feat$npeaks <- lapply(index, function(x) {
    as.data.frame(table(pk$replicate[x]))$Freq
  })

  if ("is_filled" %in% colnames(pk)) {
    feat$hasFilled <- unlist(lapply(index, function(x) 1 %in% pk$is_filled[x]))
  } else {
    feat$hasFilled <- FALSE
  }

  feat <- select(feat, id, mz, rt, d_ppm, d_sec, everything())

  object@features <- feat

  return(object)
}
