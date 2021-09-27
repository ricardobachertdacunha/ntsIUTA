

#' calculateFeaturesMetadata
#'
#' @param obj An \linkS4class{ntsData} object containing a features.
#' @param ID Optional, vector with features IDs or indexes
#' for calculating metadata for the specified features.
#'
#' @return An \linkS4class{ntsData} object containing a features and peaks
#' quality metadata.
#'
#' @export
#'
calculateFeaturesMetadata <- function(obj, ID = NULL) {

  obj2 <- obj

  #runs calculation for peaks of certain features
  if (!is.null(ID)) obj2 <- obj2[, ID]

  pat <- obj2@patdata

  patpeaks <- pat@features

  #calculates qualities of each peak in each sample
  patpeaks <- patRoon::calculatePeakQualities(patpeaks,
                                weights = NULL,
                                flatnessFactor = 0.05, #Passed to MetaClean as the flatness.factor argument to calculateJaggedness and calculateModality.
                                #avgFunc = mean, #not used for features object, --- mean additional parameter for handling featureGroups
                                parallel = TRUE)

  pkq <- patRoon::featureTable(patpeaks)

  pkq <- data.table::rbindlist(pkq, idcol = "sample")

  pkq <- dplyr::rename(pkq, feature = group)

  pkq$ID <- obj2@peaks$ID

  pkq <- pkq[, colnames(pkq)[c(which(colnames(pkq) == "ID"):ncol(pkq))], with = FALSE]

  #consolidate to features, applying either max or min for each parameter
  # TODO does not handle well NAs, needs solution

  ft <- obj2@features

  #Parameters
  ft$ApexBoundaryRatio <- pkq$ApexBoundaryRatio[pkq[, .I[ifelse(length(na.omit(ApexBoundaryRatio)) == 0, unique(ApexBoundaryRatio), which.min(ApexBoundaryRatio))], by = feature]$V1]

  ft$FWHM2Base <- pkq$FWHM2Base[pkq[, .I[ifelse(length(na.omit(FWHM2Base)) == 0, unique(FWHM2Base), which.max(FWHM2Base))], by = feature]$V1]

  ft$Jaggedness <- pkq$Jaggedness[pkq[, .I[ifelse(length(na.omit(Jaggedness)) == 0, unique(Jaggedness), which.min(Jaggedness))], by = feature]$V1]

  ft$Modality <- pkq$Modality[pkq[, .I[ifelse(length(na.omit(Modality)) == 0, unique(Modality), which.min(Modality))], by = feature]$V1]

  ft$Symmetry <- pkq$Symmetry[pkq[, .I[ifelse(length(na.omit(Symmetry)) == 0, unique(Symmetry), which.max(Symmetry))], by = feature]$V1]

  ft$GaussianSimilarity <- pkq$GaussianSimilarity[pkq[, .I[ifelse(length(na.omit(GaussianSimilarity)) == 0, unique(GaussianSimilarity), which.max(GaussianSimilarity))], by = feature]$V1]

  ft$Sharpness <- pkq$Sharpness[pkq[, .I[ifelse(length(na.omit(Sharpness)) == 0, unique(Sharpness), which.max(Sharpness))], by = feature]$V1]

  #Triangle Peak Area Similarity Ratio (TPASR)
  ft$TPASR <- pkq$TPASR[pkq[, .I[ifelse(length(na.omit(TPASR)) == 0, unique(TPASR), which.max(TPASR))], by = feature]$V1]

  ft$ZigZag <- pkq$ZigZag[pkq[, .I[ifelse(length(na.omit(ZigZag)) == 0, unique(ZigZag), which.max(ZigZag))], by = feature]$V1]

  #Scores (max)
  ft$ApexBoundaryRatioScore <- pkq$ApexBoundaryRatioScore[pkq[, .I[ifelse(length(na.omit(ApexBoundaryRatioScore)) == 0, unique(ApexBoundaryRatioScore), which.max(ApexBoundaryRatioScore))], by = feature]$V1]

  ft$FWHM2BaseScore <- pkq$FWHM2BaseScore[pkq[, .I[ifelse(length(na.omit(FWHM2BaseScore)) == 0, unique(FWHM2BaseScore), which.max(FWHM2BaseScore))], by = feature]$V1]

  ft$JaggednessScore <- pkq$JaggednessScore[pkq[, .I[ifelse(length(na.omit(JaggednessScore)) == 0, unique(JaggednessScore), which.max(JaggednessScore))], by = feature]$V1]

  ft$ModalityScore <- pkq$ModalityScore[pkq[, .I[ifelse(length(na.omit(ModalityScore)) == 0, unique(ModalityScore), which.max(ModalityScore))], by = feature]$V1]

  ft$SymmetryScore <- pkq$SymmetryScore[pkq[, .I[ifelse(length(na.omit(SymmetryScore)) == 0, unique(SymmetryScore), which.max(SymmetryScore))], by = feature]$V1]

  ft$GaussianSimilarityScore <- pkq$GaussianSimilarityScore[pkq[, .I[ifelse(length(na.omit(GaussianSimilarityScore)) == 0, unique(GaussianSimilarityScore), which.max(GaussianSimilarityScore))], by = feature]$V1]

  ft$SharpnessScore <- pkq$SharpnessScore[pkq[, .I[ifelse(length(na.omit(SharpnessScore)) == 0, unique(SharpnessScore), which.max(SharpnessScore))], by = feature]$V1]

  ft$TPASRScore <- pkq$TPASRScore[pkq[, .I[ifelse(length(na.omit(TPASRScore)) == 0, unique(TPASRScore), which.max(TPASRScore))], by = feature]$V1]

  ft$ZigZagScore <- pkq$ZigZagScore[pkq[, .I[ifelse(length(na.omit(ZigZagScore)) == 0, unique(ZigZagScore), which.max(ZigZagScore))], by = feature]$V1]

  #total score
  ft$totalScore <- pkq$totalScore[pkq[, .I[ifelse(length(na.omit(totalScore)) == 0, unique(totalScore), which.max(totalScore))], by = feature]$V1]

  # TODO test to make a better total score value 
  ft$selfScore <- ft$JaggednessScore + ft$SharpnessScore + ft$GaussianSimilarityScore

  # TODO add signal-to-noise calculation from martin here


  ##updates peaks
  pk <- obj@peaks

  #add missing cols with dummy NA
  if (FALSE %in% (colnames(pkq) %in% colnames(pk))) {
    pk[colnames(pkq)[!colnames(pkq) %in% colnames(pk)]] <- NA
  }

  pk[pk$ID %in% pkq$ID, which(colnames(pk) %in% colnames(pkq))] <- pkq

  obj@peaks <- pk


  ##updates features
  ft2 <- obj@features

  #update ft2 with missing column and add NAs
  if (FALSE %in% (colnames(ft) %in% colnames(ft2))) {
    ft2[colnames(ft)[!colnames(ft) %in% colnames(ft2)]] <- NA
  }

  ft2[ft2$ID %in% ft$ID, ] <- ft

  obj@features <- ft2


  return(obj)

}
