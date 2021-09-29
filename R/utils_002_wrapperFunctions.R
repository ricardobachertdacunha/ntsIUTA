

#' createFeatures
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param excludeBlanks Logical, set to \code{TRUE} to exclude
#' blank replicate groups from annotation.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#' 
#' @note The  \linkS4class{ntsData} object should contain the parameters
#' required for creation and annotation of features.
#'
#' @return An \linkS4class{ntsData} object with peaks
#' grouped and annotated across samples as features.
#' 
#' @export
#'
createFeatures <- function(obj, excludeBlanks = FALSE, save = TRUE) {
  
  obj <- peakPicking(obj, save = save)
  
  obj <- makeFeatures(obj, save = save)
  
  obj <- annotateFeatures(obj, excludeBlanks = excludeBlanks, save = save)
  
  return(obj)
  
}
