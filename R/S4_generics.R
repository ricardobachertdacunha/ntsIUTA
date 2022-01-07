
#' projectInfo
#'
#' @param object Object to change or get the information.
#'
#' @return Get or change title, description and/or date of the a project.
#'
setGeneric("projectInfo", function(object, ...) standardGeneric("projectInfo"))

#' samples
#'
#' @param object Object to collect sample names.
#'
#' @return A vector with sample names.
#'
setGeneric("samples", function(object) standardGeneric("samples"))

#' replicates
#'
#' @param object Object to collect sample replicate names.
#'
#' @return A vector with sample replicate names.
#'
setGeneric("replicates", function(object) standardGeneric("replicates"))

#' replicates<-
#'
#' @param object Object to assign sample replicate names.
#' @param value A character vector with one or more sample replicate names.
#'
#' @return An object updated with sample replicate names.
#'
setGeneric("replicates<-", function(object, value) standardGeneric("replicates<-"))

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

#' acquisitionMethods
#'
#' @param object Object to collect method names used for data acquisition.
#'
#' @return A vector with method names.
#'
setGeneric("acquisitionMethods", function(object) standardGeneric("acquisitionMethods"))

#' acquisitionMethods<-
#'
#' @param object Object to assign acquisition method names to an object.
#' @param value A character vector with acquisition method names.
#'
#' @return An object updated with method names.
#'
setGeneric("acquisitionMethods<-", function(object, value) standardGeneric("acquisitionMethods<-"))

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
#' @param ... Other argumments.
#'
#' @return Ac object with experimental samples assigned to QC.
#'
setGeneric("QC<-", function(object, value, ...) standardGeneric("QC<-"))

#' peakPickingParameters
#'
#' @param object An object to add or get peak picking parameters.
#' @param algorithm Algorithm for the set of parameters.
#' @param param List of parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("peakPickingParameters", function(object, algorithm, param) standardGeneric("peakPickingParameters"))

#' peakGroupingParameters
#'
#' @param object An object to add or get peak grouping parameters.
#' @param algorithm Algorithm for the set of parameters.
#' @param param List of parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("peakGroupingParameters", function(object, algorithm, param) standardGeneric("peakGroupingParameters"))

#' fillMissingParameters
#'
#' @param object An object to add or get fill missing parameters.
#' @param algorithm Algorithm for the set of parameters.
#' @param param List of parameters.
#'
#' @return A list of parameters or the object with added parameters.
#'
setGeneric("fillMissingParameters", function(object, algorithm, param) standardGeneric("fillMissingParameters"))

#' annotationParameters
#'
#' @param object An object to add or get annotationParameters parameters.
#' @param algorithm Algorithm for the set of parameters.
#' @param param List of parameters.
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

#' exportPlots
#'
#' @param object A workflow object.
#' @param path The path to export the plots.
#' @param ... Other parameters used in the production of plots.
#'
#' @return A set of plots saved into the given path.
#'
setGeneric("exportPlots", function(object, path, ...) standardGeneric("exportPlots"))

#' exportCSVs
#'
#' @param object An workflow object.
#' @param path The path to export the CSV files.
#' @param ... Other parameters used in the production of CSV files.
#'
#' @return A set of CSVs saved into the given path.
#'
setGeneric("exportCSVs", function(object, path, ...) standardGeneric("exportCSVs"))
