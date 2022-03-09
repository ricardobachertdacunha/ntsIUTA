
#' projectInfo
#'
#' @param object Object to change or get the information.
#'
#' @return Get or change title, description and/or date of the a project.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("projectInfo", function(object, ...) standardGeneric("projectInfo"))

#' path
#'
#' @param object Object to change or get the information.
#' @param ... Other arguments.
#'
#' @return Get the path of the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("path", function(object, ...) standardGeneric("path"))

#' samplesTable
#'
#' @param object Object to collect samples table.
#'
#' @return A data table with sample information.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("samplesTable", function(object) standardGeneric("samplesTable"))

#' filePaths
#'
#' @param object Object to collect file paths.
#'
#' @return A vector with file paths.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("filePaths", function(object) standardGeneric("filePaths"))

#' samples
#'
#' @param object Object to collect sample names.
#'
#' @return A vector with sample names.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("samples", function(object) standardGeneric("samples"))

#' replicates
#'
#' @param object Object to collect sample replicate names.
#'
#' @return A vector with sample replicate names.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("replicates", function(object) standardGeneric("replicates"))

#' replicates<-
#'
#' @param object Object to assign sample replicate names.
#' @param value A character vector with one or more sample replicate names.
#'
#' @return An object updated with sample replicate names.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("replicates<-", function(object, value) standardGeneric("replicates<-"))

#' blanks
#'
#' @param object Object to collect blank sample replicate groups.
#'
#' @return A vector with the blank sample replicate group names.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
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
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("blanks<-", function(object, value) standardGeneric("blanks<-"))

#' acquisitionMethods
#'
#' @param object Object to collect method names used for data acquisition.
#'
#' @return A vector with method names.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("acquisitionMethods", function(object) standardGeneric("acquisitionMethods"))

#' acquisitionMethods<-
#'
#' @param object Object to assign acquisition method names to an object.
#' @param value A character vector with acquisition method names.
#'
#' @return An object updated with method names.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("acquisitionMethods<-", function(object, value) standardGeneric("acquisitionMethods<-"))

#' polarity
#'
#' @param object Object to collect polarities.
#' @param ... Other arguments.
#'
#' @return A vector with polarity modes for entries in the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("polarity", function(object, ...) standardGeneric("polarity"))

#' polarity<-
#'
#' @param object Object to assign polarity modes.
#' @param value A character vector with polarities mode to assign to the object.
#'
#' @return An object updated with polarity modes.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("polarity<-", function(object, value) standardGeneric("polarity<-"))

#' metadata
#'
#' @param x Object to collect metadata.
#' @param ... Other argumments.
#'
#' @return A \link[data.table]{data.table} with metadata in x.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("metadata", function(x, ...) standardGeneric("metadata"))

#' QC
#'
#' @param object Object to collect the samples used for QC.
#'
#' @return A vector with sample names used for QC.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
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
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("QC<-", function(object, value, ...) standardGeneric("QC<-"))

#' EICs
#'
#' @param object Object to get extracted ion chromatograms (EICs).
#'
#' @return A \link[data.table]{data.table} with EICs in the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("EICs", function(object, ...) standardGeneric("EICs"))

#' plotEICs
#'
#' @param object Object to extract extract ion chromatograms (EICs).
#'
#' @return A plot with EICs for defined samples in the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("plotEICs", function(object, ...) standardGeneric("plotEICs"))

#' TICs
#'
#' @param object Object to extract total ion chromatograms (TICs).
#'
#' @return A \link[data.table]{data.table} with TICs from the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("TICs", function(object, ...) standardGeneric("TICs"))

#' plotTICs
#'
#' @param object Object to extract total ion chromatograms (TICs).
#'
#' @return A plot with TICs from the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("plotTICs", function(object, ...) standardGeneric("plotTICs"))

#' XICs
#'
#' @param object Object to get three dimentional (\emph{m/z}, time and intensity)
#' extracted ion chromatograms (XICs).
#'
#' @return A \link[data.table]{data.table} with XICs from the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("XICs", function(object, ...) standardGeneric("XICs"))

#' plotXICs
#'
#' @param object Object to plot three dimentional (\emph{m/z}, time and intensity)
#' extracted ion chromatograms (XICs).
#'
#' @return A plot with XICs from the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("plotXICs", function(object, ...) standardGeneric("plotXICs"))

#' MS2s
#'
#' @param object Object to get level 2 fragmentation data (MS/MS or MS2).
#'
#' @return A \link[data.table]{data.table} with MS2 from the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("MS2s", function(object, ...) standardGeneric("MS2s"))

#' plotMS2s
#'
#' @param object Object to get level 2 fragmentation data (MS/MS or MS2).
#'
#' @return A plot with MS2 from the object.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("plotMS2s", function(object, ...) standardGeneric("plotMS2s"))

#' pickingParameters
#'
#' @param object An object to add or get peak picking parameter settings.
#' @param algorithm Algorithm for the parameters.
#' @param settings List of respective settings.
#'
#' @return A list of parameter settings or the object with added parameter settings.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("pickingParameters", function(object, algorithm, settings) standardGeneric("pickingParameters"))

#' groupingParameters
#'
#' @param object An object to add or get peak grouping parameter settings.
#' @param algorithm Algorithm for the parameters.
#' @param settings List of respective settings.
#'
#' @return A list of parameter settings or the object with added the parameter settings.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("groupingParameters", function(object, algorithm, settings) standardGeneric("groupingParameters"))

#' fragmentsParameters
#'
#' @param object An object to add or get fragments extraction parameter settings.
#' @param algorithm Algorithm for the parameters.
#' @param settings List of parameter settings.
#'
#' @return A list of parameters or the object with added parameters.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("fragmentsParameters", function(object, algorithm, settings) standardGeneric("fragmentsParameters"))

#' fillingParameters
#'
#' @param object An object to add or get filling parameter settings.
#' @param algorithm Algorithm for the parameters.
#' @param settings List of parameter settings.
#'
#' @return A list of parameters or the object with added parameters.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("fillingParameters", function(object, algorithm, settings) standardGeneric("fillingParameters"))

#' annotationParameters
#'
#' @param object An object to add or get annotation parameter settings.
#' @param algorithm Algorithm for the parameters.
#' @param settings List of parameter settings.
#'
#' @return A list of parameters or the object with added parameters.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("annotationParameters", function(object, algorithm, settings) standardGeneric("annotationParameters"))

#' peaks
#'
#' @param object An object to extract chromatographic peaks.
#' @param ... Other method dependent parameters.
#'
#' @return A data.table with peak details.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("peaks", function(object, ...) standardGeneric("peaks"))

#' plotPeaks
#'
#' @param object An object to with EICs for plotting chromatographic peaks.
#' @param ... Other method dependent parameters.
#'
#' @return A plot with chromatographic peaks, including integrated areas.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("plotPeaks", function(object, ...) standardGeneric("plotPeaks"))

#' mapPeaks
#'
#' @param object An object to with EICs for mapping chromatographic peaks.
#' @param ... Other method dependent parameters.
#'
#' @return A plot with mapped peaks, including time and mass peak space.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("mapPeaks", function(object, ...) standardGeneric("mapPeaks"))

#' features
#'
#' @param object An object to extract features.
#' @param ... Other method dependent parameters.
#'
#' @return A data.frame with feature details.
#'
setGeneric("features", function(object, ...) standardGeneric("features"))

#' plotFeatures
#'
#' @param object An object with features for plotting the respective peaks.
#' @param ... Other method dependent parameters.
#'
#' @return A plot with chromatographic peaks, including integrated areas,
#' belonging to given features.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("plotFeatures", function(object, ...) standardGeneric("plotFeatures"))

#' mapFeatures
#'
#' @param object An object to with EICs for mapping features (i.e., grouped chromatographic peaks).
#' @param ... Other method dependent parameters.
#'
#' @return A plot with mapped peaks, including time and mass peak space.
#'
#' @seealso See method details in the \linkS4class{ntsData} help page.
#'
setGeneric("mapFeatures", function(object, ...) standardGeneric("mapFeatures"))

#' hasAdjustedRetentionTime
#'
#' @param object A object to check for existence of adjusted retention time.
#'
#' @return A logical value.
#'
setGeneric("hasAdjustedRetentionTime", function(object) standardGeneric("hasAdjustedRetentionTime"))

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
