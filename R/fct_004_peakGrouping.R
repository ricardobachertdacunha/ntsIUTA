

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
                         algorithm = NA_character_,
                         settings = NULL,
                         simplify = TRUE,
                         save = FALSE) {

  checkmate::assertClass(object, "ntsData")

  pat <- object@pat

  if (checkmate::testClass(pat, "featureGroups")) pat <- pat@features

  if (is.na(algorithm)) algorithm <- groupingParameters(object)@algorithm

  if (is.na(algorithm)) {
    warning("Peak grouping algorihtm not defined!")
    return(object)
  }

  if (is.null(settings)) settings <- groupingParameters(object)@settings

  if (algorithm == "xcms3") {
      settings$groupParam@sampleGroups <- pat@analysisInfo$group
    if (settings$rtalign) {
      settings$preGroupParam@sampleGroups <- pat@analysisInfo$group
    }
  }

  ag <- list(obj = pat, algorithm = algorithm)

  pat <- do.call(groupFeatures, c(ag, settings, verbose = TRUE))

  object <- groupingParameters(object, algorithm = algorithm, settings = settings)

  object@pat <- pat

  object <- buildPeaksTable(object)

  object <- buildFeatureTable(object)

  object <- addAdjustedRetentionTime(object)

  if (simplify) {
    newFeats <- new("featuresSIRIUS",
      analysisInfo = pat@features@analysisInfo,
      features = pat@features@features
    )

    newFeatGroups <- new("featureGroupsSIRIUS",
      groups = pat@groups,
      groupInfo = pat@groupInfo,
      analysisInfo = pat@analysisInfo,
      features = newFeats,
      ftindex = pat@ftindex
    )

    object@pat <- newFeatGroups
  }

  if (is.logical(save)) if (save) saveObject(object = object)

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

  cat("Building table with features... ")

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
  feat$mzmin <- unlist(lapply(index, function(x) min(pk[x, mz])))
  feat$mzmax <- unlist(lapply(index, function(x) max(pk[x, mz])))
  feat$rtmin <- unlist(lapply(index, function(x) min(pk[x, rt])))
  feat$rtmax <- unlist(lapply(index, function(x) max(pk[x, rt])))

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

  cat("Done! \n")
  return(object)
}




#' @title addAdjustedRetentionTime
#'
#' @description Function to add adjusted retention time information
#' to the scans slot of the \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @importFrom checkmate testClass assertClass
#' @importMethodsFrom xcms adjustedRtime processHistory peakGroupsMatrix hasAdjustedRtime
#' @importClassesFrom xcms XCMSnExp PeakGroupsParam
#' @importClassesFrom patRoon featureGroupsXCMS3
#'
addAdjustedRetentionTime <- function(object) {

  checkmate::assertClass(object, "ntsData")

  pat <- object@pat

  if (checkmate::testClass(pat, "featureGroupsXCMS3") & xcms::hasAdjustedRtime(pat@xdata)) {

    cat("Adding adjusted retention time values... ")

    rtAdj <- xcms::adjustedRtime(pat@xdata)

    pkAdj <- xcms::processHistory(
      pat@xdata,
      type = "Retention time correction"
    )[[1]]
    pkAdj <- pkAdj@param

    addAdjPoints <- FALSE
    if (checkmate::testClass(pkAdj, "PeakGroupsParam")) {
      addAdjPoints <- TRUE
      pkAdj <- xcms::peakGroupsMatrix(pkAdj)
    }

    scans <- object@scans
    scans <- lapply(seq_len(length(scans)), function(x, scans, rtAdj, addAdjPoints, pkAdj) {

      rts <- names(rtAdj)
      rts <- stringr::str_detect(rts, paste0("F", x))
      rts <- rtAdj[rts]

      temp <- scans[[x]]
      temp[, adjustedRetentionTime := rts]
      temp[, adjustment := adjustedRetentionTime - retentionTime]

      if (addAdjPoints) {
        pk_rts <- unique(pkAdj[, x])
        pk_rts <- pk_rts[pk_rts %in% temp$retentionTime]
        temp[retentionTime %in% pk_rts, adjPoints := pk_rts]
      }

      return(temp)
    },
    scans = scans,
    rtAdj = rtAdj,
    addAdjPoints = addAdjPoints,
    pkAdj = pkAdj)

    names(scans) <- samples(object)

    object@scans <- scans

    cat("Done! \n")
    return(object)
  }

  return(object)
}




#' @title updateFeatureTable
#'
#' @description Function to update feature table from
#' peaks remaining in the \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param fast Logical, set to \code{TRUE} for a lazy update of features.
#' When \code{TRUE}, only the mz, rt, and intensity are updated.
#' When \code{FALSE}, the feature table is completly updated.
#' Note that UFIs, quality data and annotation info are lost.
#'
#' @return An \linkS4class{ntsData} object with updated features slot.
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.table groupNames
#' @importFrom data.table setnames setorder copy
#' @importFrom dplyr select everything
#'
updateFeatureTable <- function(object, fast = TRUE) {

  cat("Updating features... ")

  pat <- object@pat

  feat <- patRoon::as.data.table(pat, average = TRUE)

  feat <- setnames(feat, c("ret", "group"), c("rt", "id"))

  rpl <- unique(replicates(object))

  hasRemoved <- FALSE

  feat_org <- features(object)

  feat_org <- feat_org[id %in% feat$id, ]

  feat2org <- feat[id %in% feat_org$id, ]

  feat_org[, rt := feat2org$rt]
  feat_org[, mz := feat2org$mz]

  if (!TRUE %in% all.equal(feat_org$id, feat2org$id)) {
    warning("Mistach of features order during update
      of feature table! Feature table not updated."
    )
    return(object)
  }

  feat_org[, (rpl) := feat2org[, rpl, with = FALSE]]

  colRem <- colnames(feat_org)
  colRem <- c(which(colRem == "d_sec") + 1, which(colRem == "mzmin") - 1)
  colRem <- colnames(feat_org)[colRem[1]:colRem[2]]
  colRem <- colRem[!grepl(rpl, colRem, fixed = TRUE)]

  if (length(colRem) > 0) feat_org[, (colRem) := NULL]

  if (nrow(object@removed) > 0) {
    hasRemoved <- TRUE
    feat_rem <- object@removed
    feat_rem <- feat_rem[id %in% feat$id, ]
    feat2rem <- feat[id %in% feat_rem$id, ]
    feat_rem[, rt := feat2rem$rt]
    feat_rem[, mz := feat2rem$mz]
    feat_rem[, (rpl) := feat2rem[, rpl, with = FALSE]]
    if (length(colRem) > 0) feat_rem[, (colRem) := NULL]
  }

  if (fast) {

    object@features <- copy(feat_org)
    if (hasRemoved) object@removed <- copy(feat_rem)

    cat("Done with the fast method! \n")

    return(object)
  }

  #not fast
  if (hasRemoved) feat_org <- rbind(feat_org, feat_rem)

  ID <- patRoon::groupNames(object@pat)
  if (length(ID) != nrow(feat_org)) {
    warning("There is a mismatch in the number of features between
      ntsData and the patRoon object! Features not updated."
    )
    return(object)
  }
  feat_org <- feat_org[data.table::data.table(id = ID), on = "id"]

  rpl_sp <- lapply(rpl, function(x, st) {
    st$sample[st$replicate == x]
  }, st = samplesTable(object))

  feat_b <- patRoon::as.data.table(pat, average = FALSE)

  #update intensities sd
  feat_sd <- lapply(rpl_sp, function(x, feat_b) {
    temp <- feat_b[, x, with = FALSE]
    temp <- apply(temp, 1, function(x) sd(x) / mean(x) * 100)
    temp[is.nan(temp)] <- 0
    temp <- round(temp, digits = 0)
    return(temp)
  }, feat_b = feat_b)

  names(feat_sd) <- paste0(rpl, "_sd")
  feat_sd <- as.data.table(feat_sd)
  feat_sd[, id := feat$id]

  sd_cols <- paste0(rpl, "_sd")
  feat_org[, (sd_cols) := feat_sd[which(id %in% feat_org$id), sd_cols, with = FALSE]]

  #update mz and rt ranges
  pk <- peaks(object)

  pk$index <- seq_len(nrow(pk))
  index <- lapply(feat_org$id, function(x) pk[feature == x, index])
  feat_org$mzmin <- unlist(lapply(index, function(x) min(pk[x, mz])))
  feat_org$mzmax <- unlist(lapply(index, function(x) max(pk[x, mz])))
  feat_org$rtmin <- unlist(lapply(index, function(x) min(pk[x, rt])))
  feat_org$rtmax <- unlist(lapply(index, function(x) max(pk[x, rt])))

  feat_org$d_ppm <- round(((feat_org$mzmax - feat_org$mzmin) / feat_org$mz) * 1E6, digits = 1)

  feat_org$d_sec <- round(feat_org$rtmax - feat_org$rtmin, digits = 0)

  feat_org$p_id <- I(lapply(feat_org$id, function(x) pk[feature == x, id]))

  pk$replicate <- factor(pk$replicate, levels = rpl)

  feat_org$npeaks <- lapply(index, function(x) {
    as.data.frame(table(pk$replicate[x]))$Freq
  })

  if ("is_filled" %in% colnames(pk)) {
    feat_org$hasFilled <- unlist(lapply(index, function(x) 1 %in% pk$is_filled[x]))
  } else {
    feat_org$hasFilled <- FALSE
  }

  object@features <- feat_org

  if (hasRemoved) object <- removeFilteredFeatures(object)

  cat("Done with complete updating method! \n")
  
  return(object)
}
