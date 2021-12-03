

#' @title calculateSNR
#' 
#' @description Calculates the signal-to-noise (sn) ratio for features in a given \linkS4class{ntsData} object
#' or specified via the argument \code{ID}. The noise is extimated by the maximum intensity of the background
#' centroids, which belong to the same mass bin as the peaks in the feature but are taken from a predefined
#' (\code{rtWindow}) time window before and after the feature limits. The sn is then calculated by
#' dividind the maximum intensity of the feature in the \linkS4class{ntsData} object with the estimated noise.
#' Centroids that belong to other peaks within the same mass and time windows are excluded before noise estimation.
#' 
#' @param obj An \linkS4class{ntsData} object.
#' @param ID Optionally, a character vector with the ID of the features to calculate the sn ratio.
#' @param rtWindow The time window (before and after), in seconds,
#' to find centroids for estimation of the noise level.
#'
#' @return An \linkS4class{ntsData} object with the noise, noise sd and sn ratio columns ammended to the
#' features slot.
#' 
#' @export
#' 
#' @importFrom checkmate assertClass
#' @importFrom stats sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
calculateSNR <-  function(obj, ID = NULL, rtWindow = 120) {
  
  assertClass(obj, "ntsData")
  
  feat_org <- obj@features
  
  peak_org <- obj@peaks
  
  if (!("ncent" %in% colnames(feat_org))) feat_org$ncent <- NA
  if (!("noise" %in% colnames(feat_org))) feat_org$noise <- NA
  if (!("noise_sd" %in% colnames(feat_org))) feat_org$noise_sd <- NA
  if (!("sn" %in% colnames(feat_org))) feat_org$sn <- NA
  
  if (!("ncent" %in% colnames(peak_org))) peak_org$ncent <- NA
  if (!("noise" %in% colnames(peak_org))) peak_org$noise <- NA
  if (!("noise_sd" %in% colnames(peak_org))) peak_org$noise_sd <- NA
  if (!("sn2" %in% colnames(peak_org))) peak_org$sn2 <- NA
  
  feat_sn <- feat_org
  
  peak_sn <- peak_org
  
  if (!is.null(ID)) feat_sn <- feat_org[feat_org$ID %in% ID, ]
  
  if (!is.null(ID)) peak_sn <- peak_org[peak_org$feature %in% ID, ]
  
  pb <- txtProgressBar(min = 0, max = nrow(feat_sn), style = 3)
  
  for (jj in  seq_len(nrow(feat_sn))) {
    
    idf <- feat_sn$ID[jj]
    idp <- peak_sn$ID[peak_sn$feature %in% idf]
    loop <- seq_len(length(idp))
    
    featEIC <- extractEIC(obj = obj,
                         samples = NULL,
                         mz = c(feat_sn$mzmin[jj], feat_sn$mzmax[jj]), 
                         ppm = NULL,
                         rt = NULL,
                         rtWindow = c((feat_sn$rtmin[jj] - rtWindow), (feat_sn$rtmax[jj] + rtWindow)),
                         rtUnit = "sec")
    
    featEIC <- split(featEIC, featEIC$file)
    names(featEIC) <- samples(obj)[as.numeric(names(featEIC))]
    
    int <- lapply(loop, function(x) peak_sn$intensity[peak_sn$ID == idp[x]])
    names(int) <- idp
    
    ncent <- lapply(loop, function(x, idp, peak_sn, featEIC) {
      cents <- peak_sn[peak_sn$ID %in% idp[x], ]
      cents2 <- featEIC[[cents$sample]]
      cents2 <- cents2[cents2$rt >= cents$rtmin & cents2$rt <= cents$rtmax, ]
      cents2 <- nrow(cents2)
      return(cents2)
      
    }, idp = idp, peak_sn = peak_sn, featEIC = featEIC)
    
    noise <- lapply(loop, function(x, idp, peak_sn, featEIC, all_peaks) {
      temp <- peak_sn[peak_sn$ID %in% idp[x], ]
      temp2 <- featEIC[[temp$sample]]
      temp2 <- temp2[!(temp2$rt >= temp$rtmin & temp2$rt <= temp$rtmax), ]
      
      if (nrow(temp2) > 1) {
        #remove other peaks within the same mass and time deviation
        others <- all_peaks[(all_peaks$sample %in% temp$sample &
                             all_peaks$rt > min(temp2$rt) &
                             all_peaks$rt < max(temp2$rt) &
                             all_peaks$mz > min(temp2$mz) &
                             all_peaks$mz < max(temp2$mz)), ]
  
        others <- others[others$ID != idp[x], ]
        
        if (nrow(others) > 0) {
          for (o in seq_len(nrow(others))) {
            temp2 <- temp2[!(temp2$rt > others$rtmin[o] & temp2$rt < others$rtmax[o]), ]
          }
        }
      }
      
      temp2 <- temp2$i
      
      return(temp2)
      
    }, idp = idp, peak_sn = peak_sn, featEIC = featEIC, all_peaks = obj@peaks)
    names(noise) <- idp
    
    
    noise_sd <- lapply(noise, function(x) sd(x, na.rm = TRUE))
    
    noise <- lapply(noise, function(x) ifelse(length(x) > 0, max(x, na.rm = TRUE), 0))
    
    
    #calculate sn and add values to peak_sn
    sn <- lapply(loop, function(x, idp, int, ncent, noise, noise_sd) {
      
      if (noise[[x]] != 0) {
        temp <- int[[x]] / noise[[x]]
      } else {
        temp <- NA
      }
      
      temp <- round(temp, digits = 0)
      
      peak_sn$ncent[peak_sn$ID == idp[x]] <<- ncent[[x]]
      peak_sn$noise_sd[peak_sn$ID == idp[x]] <<- round(noise_sd[[x]], digits = 0)
      peak_sn$noise[peak_sn$ID == idp[x]] <<- round(noise[[x]], digits = 0)
      peak_sn$sn2[peak_sn$ID == idp[x]] <<- temp
      
      return(temp)
      
    },idp = idp, int = int, ncent = ncent, noise = noise, noise_sd = noise_sd)
    
    
    toFeat <- peak_sn[peak_sn$ID %in% idp, c("ncent", "noise", "noise_sd", "sn2")]
    toFeat <- toFeat[toFeat$sn2 == max(toFeat$sn2), ]
    toFeat <- toFeat[1, ]
    
    feat_sn[jj, c("ncent", "noise", "noise_sd", "sn")] <- toFeat
    
    setTxtProgressBar(pb, jj)
    
  }
  
  
  peak_org[peak_org$ID %in% peak_sn$ID, ] <- peak_sn
  feat_org[feat_org$ID %in% feat_sn$ID, ] <- feat_sn
  
  
  obj@peaks <- peak_org
  obj@features <- feat_org

  
  return(obj)
  
}




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
