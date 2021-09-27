

#' samples
#'
#' @param object Object to collect sample names.
#'
#' @return A vector with sample names.
#'
setGeneric("samples", function(object) standardGeneric("samples"))

#' sampleGroups
#'
#' @param object Object to collect sample replicate groups. 
#'
#' @return A vector with sample replicate group names.
#'
setGeneric("sampleGroups", function(object) standardGeneric("sampleGroups"))

#' sampleGroups<-
#'
#' @param object Object to assign sample replicate groups. 
#' @param value A character vector with one or more sample replicate group names.
#'
#' @return An object updated with sample replicate group names.
#'
setGeneric("sampleGroups<-", function(object, value) standardGeneric("sampleGroups<-"))

#' blanks
#'
#' @param object Object to collect blank sample replicate groups.
#'
#' @return A vector with the blank sample replicate group names.
#'
setGeneric("blanks", function(object) standardGeneric("blanks"))

#' blanks<-
#'
#' @param object Object to assign blank sample replicate groups.
#' @param value A character vector with names of the sample replicate groups
#' to be used as blanks.
#'
#' @return An object updated with blank sample replicate groups.
#'
setGeneric("blanks<-", function(object, value) standardGeneric("blanks<-"))

#' QC
#'
#' @param object Object to collect the samples used for QC.
#'
#' @return A vector with sample names used for QC.
#'
setGeneric("QC", function(object) standardGeneric("QC"))

#' QC<-
#'
#' @param object Object to collect the samples to be used as QC.
#' @param value A character vector the names of the sample replicate group/s
#' to be used as QC.
#'
#' @return Ac object with experimental samples assigned to QC.
#'
setGeneric("QC<-", function(object, value) standardGeneric("QC<-"))

#' peakPickingParameters
#'
#' @param object An object to add or get peak picking parameters.
#' @param ... Other method dependent parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("peakPickingParameters", function(object, algorithm, param) standardGeneric("peakPickingParameters"))

#' peakGroupingParameters
#'
#' @param object An object to add or get peak grouping parameters.
#' @param ... Other method dependent parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("peakGroupingParameters", function(object, algorithm, param) standardGeneric("peakGroupingParameters"))

#' fillMissingParameters
#'
#' @param object An object to add or get fill missing parameters.
#' @param ... Other method dependent parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("fillMissingParameters", function(object, algorithm, param) standardGeneric("fillMissingParameters"))

#' annotationParameters
#'
#' @param object An object to add or get annotationParameters parameters.
#' @param ... Other method dependent parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("annotationParameters", function(object, algorithm, param) standardGeneric("annotationParameters"))

#' peaks
#'
#' @param object An object to extract chromatographic peaks.
#' @param ... Other method dependent parameters.
#'
#' @return A data.frame with peak details.
#'
setGeneric("peaks", function(object, ...) standardGeneric("peaks"))

#' features
#'
#' @param object An object to extract features.
#' @param ... Other method dependent parameters.
#'
#' @return A data.frame with feature details.
#'
setGeneric("features", function(object, ...) standardGeneric("features"))

#' components
#'
#' @param object An object to extract feature components.
#' @param ... Other method dependent parameters.
#'
#' @return A data.frame with feature component details.
#'
setGeneric("components", function(object, ...) standardGeneric("components"))
