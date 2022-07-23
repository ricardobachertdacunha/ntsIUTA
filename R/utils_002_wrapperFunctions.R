

#' @title createFeatures
#'
#' @description Wrapper function for peak picking, grouping, filling, annotation
#' and creation of UFIs when settings are already defined in the \linkS4class{ntsParameters}.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @note The \linkS4class{ntsData} object should contain the parameter settings
#' required for creation and annotation of features.
#'
#' @return An \linkS4class{ntsData} object with peaks
#' grouped and annotated across samples as features.
#'
#' @export
#'
createFeatures <- function(object, save = FALSE) {

  object <- peakPicking(object, save = save)

  object <- peakGrouping(object, save = save)

  object <- peakFilling(object, save = save)

  object <- peakAnnotation(object, save = save)

  object <- makeUFI(object, save = save)

  if (is.logical(save)) if (save) saveObject(object = object)

  if (is.character(save)) saveObject(object = object, filename = save)

  return(object)
}
