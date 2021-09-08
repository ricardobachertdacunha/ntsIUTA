
#' filterMinInt
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param IntThreshold A numerical value set at the desired minimum intensity set for features.
#'
#' @return
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.frame 
#' 
filterMinInt <- function(obj, IntThreshold = 500) {

  feats <- patRoon::as.data.frame(obj@features)

  if (!("isFiltered" %in% colnames(feats))) {
    feats$isFiltered <- FALSE
  }
  
  if (!("filteredBy" %in% colnames(feats))) {
    feats$filteredBy <- ""
  }

  lastCol <- length(unique(sampleGroups(obj)))+5

  feats <- feats[feats$isFiltered == FALSE,]

  # for( j in 1:nrow(feats)) { 
  #   vari <- feats[j,i] <= IntThreshold
  #   if(!(is.na(vari))){
  #     if (vari) { feats$filtered[j] <- TRUE
  #                 feats$filteredBy[j] <- "MinInt" }
  #   }
  # }
  
    # for(j in 1:nrow(feats)) { 
    # 
    # vari <- feats[j,6:lastCol]
    # 
    #   if(!(is.na(vari))){
    #   
    #    if ((TRUE %in% (vari < IntThreshold))) { 
    #      feats$isFiltered[j] <- TRUE
    #      feats$filteredBy[j] <- "MinInt"}
    #   }
    # }
  
    for(i in 1:nrow(feats)) {
      if(all((feats[i,6:lastCol] < IntThreshold), na.rm = TRUE)) {
        feats$isFiltered[i] <- TRUE
        feats$filteredBy[i] <- "MinInt"
        }
    }


  obj@features <- feats

  return(obj)
}


#' filterBlank
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param blankThreshold A numerical value set at the desired minimum intensity set for features.
#'
#' @return
#'
#' @export
#' 
#' @examples
#'

filterBlank <- function(obj, blankThreshold = 3) {
  
  feats <-obj@features
  
  if (!("isFiltered" %in% colnames(feats))) {
    feats$isFiltered <- FALSE
  }
  
  if (!("filteredBy" %in% colnames(feats))) {
    feats$filteredBy <- ""
  }
  
  lastCol <- length(unique(sampleGroups(obj)))+5
  
  feats <- feats[feats$isFiltered == FALSE,]

  for(i in 1:nrow(feats)) {
    
 #   blnk <- feats[i,6] * blankThreshold
    

    if(all((feats[i,7:lastCol] < (feats[i,6] * blankThreshold)), na.rm = TRUE)) {
      
      feats$isFiltered[i] <- TRUE
      feats$filteredBy[i] <- "blankThreshold"
    }
  }
  
  
  obj@features <- feats
  
  return(obj)
}




# filterPeaks <- function(peaks = peaks,
#                         fileIndex = NULL,
#                         mz = NULL,
#                         ppm = NULL,
#                         rt = NULL,
#                         rtWindow = NULL,
#                         rtUnit = "sec",
#                         absMinIntensity = NULL,
#                         chromWidthRange = NULL,
#                         negate = FALSE,
#                         dataframe = TRUE) {
#
#   x <- peaks
#
#   x <- x[fileIndex]
#
#   mzr <- mzrBuilder(mz = mz, ppm = ppm)
#
#   rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
#
#   x <- patRoon::filter(obj = x,
#                        absMinIntensity = absMinIntensity,
#                        relMinIntensity = NULL,
#                        retentionRange = rtr,
#                        mzRange = mzr,
#                        mzDefectRange = NULL,
#                        chromWidthRange = chromWidthRange,
#                        negate = negate)
#
#   if (dataframe) x <- patRoon::as.data.frame(x)
#
#   return(x)
#
# }
