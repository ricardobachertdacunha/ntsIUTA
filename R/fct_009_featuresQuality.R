

#' @title calculateSNR
#'
#' @description Calculates the signal-to-noise (sn) ratio for features in a given \linkS4class{ntsData} object
#' or specified via the argument \code{targets}. The noise is extimated by the maximum intensity of the background
#' centroids (i.e., traces belonging to the same mass bin as the peaks in the feature but taken from a wider and predefined
#' time window (\code{rtExpand}) before and after the feature limits. The sn is then calculated by
#' dividind the maximum intensity of the feature in the \linkS4class{ntsData} object with the estimated noise.
#' Centroids that belong to other peaks within the same mass and time windows are excluded before noise estimation.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param targets Optionally, a character vector with the \emph{id} of the features to calculate the sn ratio.
#' @param rtExpand The time, in seconds,
#' expansion before and after the time window of each feature
#' to extract centroids for estimation of the noise level.
#'
#' @return An \linkS4class{ntsData} object with the noise, noise sd and sn ratio columns ammended to the
#' features slot.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom stats sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importFrom data.table copy rbindlist
#'
calculateSNR <-  function(object, targets = NULL, rtExpand = 200) {

  assertClass(object, "ntsData")

  feat_org <- features(object)

  peak_org <- peaks(object)

  if (!("sn_pN" %in% colnames(feat_org))) feat_org[, sn_pN := NA]
  if (!("sn_nN" %in% colnames(feat_org))) feat_org[, sn_nN := NA]
  if (!("sn_noise" %in% colnames(feat_org))) feat_org[, sn_noise := NA]
  if (!("sn_noise_sd" %in% colnames(feat_org))) feat_org[, sn_noise_sd := NA]
  if (!("sn_value" %in% colnames(feat_org))) feat_org[, sn_value := NA]

  if (!("sn_pN" %in% colnames(peak_org))) peak_org[, sn_pN := NA]
  if (!("sn_nN" %in% colnames(peak_org))) peak_org[, sn_nN := NA]
  if (!("sn_noise" %in% colnames(peak_org))) peak_org[, sn_noise := NA]
  if (!("sn_noise_sd" %in% colnames(peak_org))) peak_org[, sn_noise_sd := NA]
  if (!("sn_value" %in% colnames(peak_org))) peak_org[, sn_value := NA]

  peak_org$sn_value <- as.numeric(peak_org$sn_value)
  peak_org$sn_pN <- as.numeric(peak_org$sn_pN)
  peak_org$sn_nN <- as.numeric(peak_org$sn_nN)
  peak_org$sn_noise <- as.numeric(peak_org$sn_noise)
  peak_org$sn_noise_sd <- as.numeric(peak_org$sn_noise_sd)

  feat_org$sn_value <- as.numeric(feat_org$sn_value)
  feat_org$sn_pN <- as.numeric(feat_org$sn_pN)
  feat_org$sn_nN <- as.numeric(feat_org$sn_nN)
  feat_org$sn_noise <- as.numeric(feat_org$sn_noise)
  feat_org$sn_noise_sd <- as.numeric(feat_org$sn_noise_sd)

  feat_sn <- copy(feat_org)
  peak_sn <- copy(peak_org)

  if (!is.null(targets)) feat_sn <- feat_sn[id %in% targets, ]
  if (!is.null(targets)) peak_sn <- peak_sn[feature %in% targets, ]

  #extract centroids from each peak in each sample, expanding the rt
  sName <- unique(peak_sn$sample)
  pks <- peak_sn[, .(id, mzmin, mzmax, rtmin, rtmax, sample)]
  pks <- pks[, `:=`(rtmin = rtmin - rtExpand, rtmax = rtmax + rtExpand)]

  eic <- lapply(sName, function(x, pks, object) {
    temp <- extractEICs(
      object,
      samples =  x,
      mz = pks[sample %in% x, ]
    )
  }, pks = pks, object = object)
  eic <- rbindlist(eic)
  ncent <- split(eic, by = "id")
  cat("\n")

  cat(paste0("Calculating signal-to-noise ratio for ", nrow(feat_sn), " features...\n"))

  pb <- txtProgressBar(
    min = 0,
    max = nrow(feat_sn),
    style = 3,
    width = 50,
    char = "="
  )

  for (jj in  seq_len(nrow(feat_sn))) {

    idf <- feat_sn$id[jj]
    idp <- peak_sn[feature %in% idf, id]

    int <- lapply(idp, function(x) peak_sn[id == x, intensity])
    names(int) <- idp

    noise_data_raw <- lapply(idp, function(x, ncent, peak_sn) {

      temp <- peak_sn[id %in% x, ]
      temp2 <- ncent[[x]]
      temp2 <- temp2[!(temp2$rt >= temp$rtmin & temp2$rt <= temp$rtmax), ]

      if (nrow(temp2) > 1) {
        #remove other peaks within the same mass and time deviation
        others <- peak_sn[
          sample %in% temp$sample &
          rt >= min(temp2$rt) &
          rt <= max(temp2$rt) &
          mz >= temp$mzmin &
          mz <= temp$mzmax &
          !(id %in% x),
        ]

        if (nrow(others) > 0) {
          for (oo in seq_len(nrow(others))) {
            temp2 <- temp2[!(rt > others$rtmin[oo] & rt < others$rtmax[oo]), ]
          }
        }
      }

      temp2 <- temp2$intensity

      return(temp2)
    }, peak_sn = peak_sn, ncent = ncent)

    names(noise_data_raw) <- idp

    noise_data_sd <- lapply(noise_data_raw, function(x) sd(x, na.rm = TRUE))

    noise_data <- lapply(noise_data_raw, function(x) ifelse(length(x) > 0, max(x, na.rm = TRUE), 0))


    #calculate sn and add values to peak_sn
    for (pp in idp) {

      if (noise_data[[pp]] != 0) {
        peak_sn[id %in% pp, sn_value := round(int[[pp]] / noise_data[[pp]], digits = 0)]
      } else {
        peak_sn[id %in% pp, sn_value := NA]
      }

      peakcentN <- nrow(eic[id %in% pp & rt >= peak_sn[id %in% pp, rtmin] & rt <= peak_sn[id %in% pp, rtmax], ])
      noisecentN <- sapply(noise_data_raw, function(x) length(x))
      peak_sn[id %in% pp, sn_pN := peakcentN]
      peak_sn[id %in% pp, sn_nN := noisecentN[pp]]

      peak_sn[id %in% pp, sn_noise_sd := round(noise_data_sd[[pp]], digits = 0)]
      peak_sn[id %in% pp, sn_noise := round(noise_data[[pp]], digits = 0)]
    }

    toFeat <- peak_sn[id %in% idp, .(sn_pN, sn_nN, sn_noise, sn_noise_sd, sn_value)]
    toFeat <- toFeat[sn_value == max(sn_value, na.rm = TRUE), ]
    toFeat <- toFeat[1, ]

    feat_sn[id %in% idf, c("sn_pN", "sn_nN", "sn_noise", "sn_noise_sd", "sn_value")] <- toFeat

    setTxtProgressBar(pb, jj)
  }
  close(pb)
  cat("Done! \n")

  peak_org[id %in% peak_sn$id, ] <- peak_sn
  feat_org[id %in% feat_sn$id, ] <- feat_sn

  object@peaks <- peak_org
  object@features <- feat_org

  return(object)
}




#' @title calculateFeaturesMetadata
#'
#' @description Function to calculate feature quality indicators
#' using the \pkg{MetaClean} package via the \pkg{patRoon} interface.
#' See \code{?MetaClean} for further details.
#'
#' @param object An \linkS4class{ntsData} object containing features.
#' @param targets Optional character vector with features \emph{id}/s
#' for calculating metadata for the specified features.
#'
#' @return An \linkS4class{ntsData} object containing features and peaks
#' with amended quality indicators.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importMethodsFrom patRoon calculatePeakQualities featureTable
#' @importFrom data.table rbindlist
#'
calculateFeaturesMetadata <- function(object, targets = NULL) {

  checkmate::assertClass(object, "ntsData")

  object2 <- object

  #runs calculation for peaks of certain features
  if (!is.null(targets)) object2 <- object2[, targets]

  pat <- object2@pat
  patpeaks <- pat@features

  patpeaks <- patRoon::calculatePeakQualities(
    patpeaks,
    weights = NULL,
    flatnessFactor = 0.05, #Passed to MetaClean as the flatness.factor argument to calculateJaggedness and calculateModality.
    #avgFunc = mean, #not used for features object, --- mean additional parameter for handling featureGroups
    parallel = TRUE
  )

  pkq <- patRoon::featureTable(patpeaks)
  pkq <- rbindlist(pkq, idcol = "sample")
  newCols <- colnames(pkq)[c(which(colnames(pkq) == "ApexBoundaryRatio"):ncol(pkq))]

  pks <- peaks(object2)
  pks <- pks[, colnames(pks)[!colnames(pks) %in% newCols], with = FALSE]

  if (FALSE %in% (all.equal(pks$feature, pkq$group) & all.equal(pks$mz, pkq$mz))) {
    warning("Peaks table do not match between patRoon and ntsIUTA!")
    return(object)
  }

  pks <- cbind(pks, pkq[, newCols, with = FALSE])

  pks2 <- copy(pks)

  pks2 <- pks2[, .(
    ApexBoundaryRatio = min(ApexBoundaryRatio, na.rm = TRUE),
    FWHM2Base = max(FWHM2Base, na.rm = TRUE),
    Jaggedness = min(Jaggedness, na.rm = TRUE),
    Modality = min(Modality, na.rm = TRUE),
    Symmetry = max(Symmetry, na.rm = TRUE),
    GaussianSimilarity = max(GaussianSimilarity, na.rm = TRUE),
    Sharpness = max(Sharpness, na.rm = TRUE),
    TPASR = max(TPASR, na.rm = TRUE),
    ZigZag = min(ZigZag, na.rm = TRUE),
    ApexBoundaryRatioScore = max(ApexBoundaryRatioScore, na.rm = TRUE),
    FWHM2BaseScore = max(FWHM2BaseScore, na.rm = TRUE),
    JaggednessScore = max(JaggednessScore, na.rm = TRUE),
    ModalityScore = max(ModalityScore, na.rm = TRUE),
    SymmetryScore = max(SymmetryScore, na.rm = TRUE),
    GaussianSimilarityScore = max(GaussianSimilarityScore, na.rm = TRUE),
    SharpnessScore = max(SharpnessScore, na.rm = TRUE),
    TPASRScore = max(TPASRScore, na.rm = TRUE),
    ZigZagScore = max(ZigZagScore, na.rm = TRUE),
    totalScore = max(totalScore, na.rm = TRUE)
  ), by = feature]

  #update peaks newCols
  peaks_org <- copy(object@peaks)

  if (!TRUE %in% newCols %in% colnames(peaks_org)) {
    peaks_org[, (newCols) := as.numeric(NA)]
  }

  peaks_org[which(id %in% pks$id), (newCols) := pks[, newCols, with = FALSE]]

  object@peaks <- copy(peaks_org)

  #update features newCols
  feats_org <- copy(object@features)

  if (!TRUE %in% newCols %in% colnames(feats_org)) {
    feats_org[, (newCols) := as.numeric(NA)]
  }

  feats_org[which(feats_org$id %in% pks2$feature), (newCols) := pks2[, newCols, with = FALSE]]

  object@features <- copy(feats_org)

  return(object)
}
