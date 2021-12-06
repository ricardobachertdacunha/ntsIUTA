

#' filterFeatures
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param ... A sequence of filter arguments to applied to the features data frame.
#' The filter sequence order is user defined and is applied as given in the arguments.
#' See details for possible filter arguments.
#'
#' @return An \linkS4class{ntsData} object with filtered features annotated in the features slot.
#' 
#' @details The available filters are as folllows:
#' \itemize{
#'  \item \code{filterMinInt} Filter features below a minimum intensity threshold. For example, filterMinInt = 3000, removing features below 3000 counts.
#'  \item \code{filterBlank} Filter features that are not more intense than a defined multiplier times the assigned blank intensity. For example, filterBlank = 3, removing features that are not higher than 3 times the blank intensity.
#'  \item \code{filterSNR} Filter features below a minimum signal-to-noise ratio (SNR) threshold. For example, filterSNR = 3, filtering features below a SNR of 3.
#'  \item \code{filterMinReplicateAbundance} Filter features that are not present with at least a specified frequency in one sample replicate group. For example, filterMinReplicateAbundance = 2, meaning that featues should be at least represented in two samples of at least one sample replicate group.
#'  \item \code{filterIntSD} Filter features based on a maximum standard deviation (SD) among replicate samples. The SD should be at least below the defined threshold, as percentage, in at least one sample replicate group. For example, filterIntSD = 30, filtering features that do not have the SD below 30% in at least one sample replicate group.
#' }
#'
#' @export
#'
#' 
filterFeatures <- function(obj, ...) {
  
  filterList <- list(...)
  
  if (length(filterList) == 0) {
    warning("No filters selected or recognized.")
    return(obj)
  }
  
  filters <- names(filterList)
  
  listOfViableFilters <- c("filterMinInt",
                           "filterBlank",
                           "filterSNR",
                           "filterMinReplicateAbundance",
                           "filterIntSD")
  
  if (!all(filters %in% listOfViableFilters)) {
    warning("At least one filters is not recognized.")
    return(obj)
  }
  
  for(i in seq_len(length(filters))) {
    switch(names(filterList)[i],
      filterMinInt = (obj <- filterMinInt(obj, unlist(filterList[i]))),
      filterBlank = (obj <- filterBlank(obj, unlist(filterList[i]))),
      filterSNR = (obj <- filterSNR(obj, unlist(filterList[i]))),
      filterMinReplicateAbundance = (obj <- filterMinReplicateAbundance(obj, unlist(filterList[i]))),
      filterIntSD = (obj <- filterIntSD(obj, unlist(filterList[i])))
    )
  }

  obj@filters <- c(obj@filters, filterList)
  
  return(obj)

}




#' removeFilteredFeatures
#'
#' @param obj An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with filtered features moved to the removed slot.
#' 
#' @export
#' 
removeFilteredFeatures <- function(obj) {
  
  feats_org <- obj@features
  
  removed <- feats_org[feats_org$isFiltered, ]
  
  feats <- feats_org[!feats_org$isFiltered, ]
  
  if (nrow(obj@removed) > 0) { #when there are removed features already
    
    #add columns in to obj removed that were not there before
    obj@removed[, colnames(removed)[!(colnames(removed) %in% colnames(obj@removed))]] <- NA
    
    obj@removed <- rbind(obj@removed, removed)
    
  } else {
    
    obj@removed <- removed
    
  }
  
  obj@features <- feats
  
  return(obj)
  
}




#' restoreFilteredFeatures
#'
#' @param obj An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with filtered features restored to the features slot.
#' 
#' @importMethodsFrom patRoon groupNames
#' 
#' @export
#' 
restoreFilteredFeatures <- function(obj) {
  
  feats <- obj@features
  
  removed <- obj@removed
  
  if (nrow(removed) == 0) return(obj)
  
  feats_org <- rbind(feats, removed)
  
  ID <- patRoon::groupNames(obj@patdata)
  
  feats_org <- feats_org[match(ID, feats_org$ID), ]
  
  obj@removed <- data.frame()
  
  obj@features <- feats_org
  
  obj@filters <- list()
  
  return(obj)
  
}
