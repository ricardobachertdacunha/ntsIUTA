
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
#' @examples
#'

filterMinInt <- function(obj, IntThreshold = 500) {

  feats <- patRoon::as.data.frame(obj@features)

  if (!("filtered" %in% colnames(feats))) {
    feats$filtered <- FALSE
  }
  if (!("filteredBy" %in% colnames(feats))) {
    feats$filteredBy <- ""
  }
# TODO select only unfiltered data
for (i in 6:(length(unique(sampleGroups(obj)))-1)) {
  
  for( j in 1:nrow(feats)) { 
    
    vari <- feats[j,i] <= IntThreshold
    
    if(!(is.na(vari))){
      
      if (vari) { feats$filtered[j] <- TRUE
                  feats$filteredBy[j] <- "MinInt" }
    }
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
