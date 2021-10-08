

#' filterMinInt
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param intensityThreshold A numerical value set at the desired minimum intensity for features.
#'
#' @return
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.frame 
#' 
filterMinInt <- function(obj, intensityThreshold = 500) {
  
  feats_org <- obj@features
  
  if (!("isFiltered" %in% colnames(feats_org))) feats_org$isFiltered <- FALSE
  
  if (!("filteredBy" %in% colnames(feats_org))) feats_org$filteredBy <- NA_character_
  
  feats <- feats_org[feats_org$isFiltered == FALSE, ]
  
  check <- apply(feats[, unique(sampleGroups(obj))], MARGIN = 1, function(x){
    all(x < intensityThreshold, na.rm = TRUE)
  })
  
  feats$isFiltered <- check
  
  feats$filteredBy[check] <- "MinInt"
  
  feats_org[feats_org$ID %in% feats$ID, ] <- feats
  
  obj@features <- feats_org
  
  return(obj)
}




#' filterBlank
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param blankMultiplier A numerical value to multiply the assigned blank intensity.
#'
#' @return
#'
#' @export
#' 
filterBlank <- function(obj, blankMultiplier = 3) {
  
  feats_org <- obj@features
  
  if (!("isFiltered" %in% colnames(feats_org))) feats_org$isFiltered <- FALSE
  
  if (!("filteredBy" %in% colnames(feats_org))) feats_org$filteredBy <- NA_character_
  
  feats <- feats_org[feats_org$isFiltered == FALSE, ]
  
  blk <- obj@samples[, c("group", "blank")]
  blk <- unique(blk)
  
  rg <- blk$group[!blk$group %in% blk$blank]
  
  check <- feats[, rg, drop = FALSE]
  
  for (r in seq_len(length(rg))) {
    
    blkInt <- feats[, blk$blank[blk$group == rg[r]], drop = TRUE]
    
    check[, rg[r]] <- check[, rg[r], drop = TRUE] < (blkInt * blankMultiplier)
    
  }
  
  check <- apply(check, MARGIN = 1, function(x) all(x, na.rm = TRUE))
  
  feats$isFiltered <- check
  
  feats$filteredBy[check] <- "Blank"
  
  feats_org[feats_org$ID %in% feats$ID, ] <- feats
  
  obj@features <- feats_org
  
  return(obj)
  
}



#' filterSNR
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param snRatio A numerical value set at the desired signal-to-noise ratio for features.
#'
#' @return
#'
#' @export
#' 
filterSNR <- function(obj, snRatio = 5) {
  
  feats_org <- obj@features
  
  if (!("sn" %in% colnames(feats_org))) {
    warning("sn values not found!")
    return(obj)
  }
  
  if (!("isFiltered" %in% colnames(feats_org))) feats_org$isFiltered <- FALSE
  
  if (!("filteredBy" %in% colnames(feats_org))) feats_org$filteredBy <- NA_character_
  
  feats <- feats_org[feats_org$isFiltered == FALSE, ]
  
  feats$isFiltered <- feats$sn < snRatio
  
  feats$isFiltered[is.na(feats$isFiltered)] <- FALSE
  
  feats$filteredBy[feats$isFiltered] <- "SN"
  
  feats_org[feats_org$ID %in% feats$ID, ] <- feats
  
  obj@features <- feats_org
  
  return(obj)
  
}




#' filterMinReplicateAbundance
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param blankThreshold A numerical value set at the desired minimum intensity for features.
#'
#' @return
#'
#' @export
#' 
filterMinReplicateAbundance <- function(obj, ReplicateThreshold = 2) {
  
  feats_org <- obj@features
  
  if (!("isFiltered" %in% colnames(feats_org))) feats_org$isFiltered <- FALSE
  
  if (!("filteredBy" %in% colnames(feats_org))) feats_org$filteredBy <- NA_character_
  
  feats <- feats_org[feats_org$isFiltered == FALSE, ]
  
  check <- feats$npeaks
  
  feats$isFiltered <- unlist(sapply(check, function(x) all(x < ReplicateThreshold)))
  
  feats$filteredBy[feats$isFiltered] <- "abundance"
  
  feats_org[feats_org$ID %in% feats$ID, ] <- feats
  
  obj@features <- feats_org
  
  return(obj)
  
}




#' filterIntSD
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param sdThresholdPerc A numerical value set at the desired sd percentage maximum.
#'
#' @return
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.frame
#' 
filterIntSD <- function(obj, sdThresholdPerc = 22) {
  
  feats_org <- obj@features
  
  if (!("isFiltered" %in% colnames(feats_org))) feats_org$isFiltered <- FALSE
  
  if (!("filteredBy" %in% colnames(feats_org))) feats_org$filteredBy <- NA_character_
  
  feats <- feats_org[feats_org$isFiltered == FALSE, ]
  
  rgSD <- paste0(unique(sampleGroups(obj)), "_sd%")
  
  check <- apply(feats[, rgSD], MARGIN = 1, function(x){
    all(x > sdThresholdPerc, na.rm = TRUE)
  })
  
  feats$isFiltered <- check
  
  feats$filteredBy[check] <- "SD"
  
  feats_org[feats_org$ID %in% feats$ID, ] <- feats
  
  obj@features <- feats_org
  
  return(obj)
}
