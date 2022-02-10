

### projectInfo ---------------------------------------------------------------------------------------------

#' @describeIn ntsData setter for project basic information.
#' When the \code{title}, \code{description} and \code{date} arguments
#' are missing it returns a list with the project basic information.
#' If arguments are given, the project basic information is updated
#' based on the specified arguments, returning the \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param title A character string to be used as title.
#' @param description A character string with a description for the project.
#' @param date \link{Date} object.
#'
#'
#' @export
#'
setMethod("projectInfo", "ntsData", function(object, title = NULL, description = NULL, date = NULL) {

  if (missing(title) & missing(description) & missing(date)) {
    info <- list(
      title = object@title,
      description = object@description,
      date = object@date
    )
    return(info)
  }

  if (!missing(title) & !is.null(title)) object@title <- title
  if (!missing(description) & !is.null(description)) object@description <- description
  if (!missing(date) & !is.null(date)) object@date <- as.Date(date)
  return(object)
})




### path ----------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for project path.
#'
#' @export
#'
setMethod("path", "ntsData", function(object) object@path)




### samplesTable --------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for samples data table.
#'
#' @export
#'
setMethod("samplesTable", "ntsData", function(object) object@samples)




### filePaths -----------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for file paths.
#'
#' @export
#'
setMethod("filePaths", "ntsData", function(object) object@samples$file)




### samples -------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for sample names.
#'
#' @export
#'
setMethod("samples", "ntsData", function(object) object@samples$sample)




### replicates ----------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for sample replicate names.
#'
#' @export
#'
setMethod("replicates", "ntsData", function(object) object@samples$replicate)

#' @describeIn ntsData Setter for sample replicate names.
#' The \code{value} is a character vector with the same length as the number of samples
#' in the \code{object}, containing sample replicate names for each sample.
#'
#' @param value A character vector applicable to the respective method.
#'
#' @export
#'
#' @importFrom data.table data.table
#'
setMethod("replicates<-", signature("ntsData", "ANY"), function(object, value) {

  if (length(value) != length(object@samples$sample)) {
    warning("Length of value does not match the number of samples.")
    return(object)
  }

  object@samples$replicate <- value
  object@metadata <- data.table(replicate = unique(value))
  return(object)
})




### blanks --------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for blank replicate names.
#'
#' @export
#'
setMethod("blanks", "ntsData", function(object) object@samples$blank)

#' @describeIn ntsData Setter for blank replicate groups.
#' The \code{value} is a character vector with the same length as the number of samples
#' in the \code{object}, containing blank replicate name to associate to each sample.
#'
#' @export
#'
setMethod("blanks<-", signature("ntsData", "ANY"), function(object, value) {

  if (FALSE %in% unique(value %in% object@samples$replicate)) {
    warning("Blank replicate sample groups must be one or more sample replicate groups.")
    return(object)
  }

  object@samples$blank <- value
  return(object)
})




### acquisitionMethods --------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for acquisition method names.
#'
#' @export
#'
setMethod("acquisitionMethods", "ntsData", function(object) {

  m <- object@samples$method

  if (length(unique(m)) == 1) {
    return(unique(m))
  } else {
    names(m) <- object@samples$sample
    return(m)
  }
})

#' @describeIn ntsData Setter for acquisition method names.
#' The \code{value} is a character vector with the same length as the number of samples
#' in the \code{object}, containing the name of the aquisition method used for each sample.
#' When the length is one, the name is applied to all samples.
#'
#' @export
#'
setMethod("acquisitionMethods<-", signature("ntsData", "ANY"), function(object, value) {

  if (length(value) != 1 & length(value) != length(samples(object))) {
    warning("The number of method names is not equal to the number of samples in the object!")
    return(object)
  }

  object@samples$method <- value
  return(object)
})




### polarity ------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for the polarity of each sample/replicate.
#' The \code{groupBy} argument is either \code{samples} (the default) or \code{replicates}
#' as character string to return the polarites either for each sample or each replicate.
#'
#' @param groupBy A length one character string.
#'
#' @export
#'
setMethod("polarity", "ntsData", function(object, colorBy) {

  if (missing(colorBy)) colorBy <- "samples"

  if (colorBy == "samples") {
    m <- object@samples$polarity
    names(m) <- object@samples$sample
    return(m)
  } else {
    m <- object@samples$polarity
    names(m) <- object@samples$replicate
    return(m[!duplicated(names(m))])
  }
})

#' @describeIn ntsData Setter for the polarity mode of the samples (i.e., files).
#' The \code{value} is a character vector with either \emph{positive} or \emph{negative} strings
#' with the same length as the number of samples in the set. If only one polarity mode is used for all the samples,
#' the \code{value} can be of length one and either \emph{positive} or \emph{negative}.
#' The polarity mode is used for all the samples.
#'
#' @export
#'
setMethod("polarity<-", "ntsData", function(object, value) {

  if (length(value) != 1 & length(value) != length(samples(object))) {
    warning("The length of the value argument is not equal to the number of samples in the object!")
    return(object)
  }

  if (FALSE %in% (value %in% c("positive", "negative"))) {
    warning("Only positive and negative polarity modes are possible!")
    return(object)
  }

  object@samples$polarity <- value
  return(object)
})




### metadata ------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for metadata as \link[data.table]{data.table}.
#' The \code{varname} argument is a character string to specify
#' which metadata variable/s to extract.
#'
#' @param x An \linkS4class{ntsData} object.
#' @param varname A character vector applicable to the respective method.
#'
#' @export
#'
setMethod("metadata", "ntsData", function(x, varname) {

  if (!missing(varname)) {
    if (all(varname %in% colnames(x@metadata))) {
      if (!"replicate" %in% varname) {varname <- c("replicate", varname)}
      m <- x@metadata[, varname, with = FALSE]
    } else {
      warning("varname not found in the metadata.")
      m <- x@metadata
    }
  } else {
    m <- x@metadata
  }

  return(m)
})




### QC ------------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for the QC samples \link[data.table]{data.table}.
#'
#' @export
#'
setMethod("QC", "ntsData", function(object) object@QC@samples)

#' @describeIn ntsData Setter for QC samples or sample replicates.
#' The \code{value} is a character vector with the names of the samples or sample replicate name/s
#' to be used for QC as predefined by the \code{nameType} argument when using
#' "samples" or "replicates", respectively. If the \code{remove} argument is \code{TRUE}
#' the specified sample or replicate names are moved from
#' the QC to the samples of the \linkS4class{ntsData} object.
#'
#' @param nameType A character string of length one applicable to the respective method.
#' @param remove A logical value applicable to the respective method.
#'
#' @export
#'
setMethod("QC<-", "ntsData", function(object, value, remove = FALSE, nameType = "replicates") {

  if (missing(nameType)) nameType <- "replicates"

  if (!missing(remove) & remove) {
    if (nameType == "replicates") {
      if (FALSE %in% unique(value %in% object@QC@samples$replicate)) {
        cat("Given replicate name/s in value not found in the QC slot of the object.")
      } else {
        object@samples <- rbind(object@samples, object@QC@samples[object@QC@samples$replicate %in% value, ])
        object@samples <- object@samples[order(sample)]
        object@QC@samples <- object@QC@samples[!object@QC@samples$replicate %in% value, ]
      }
    } else {
      if (FALSE %in% unique(value %in% object@QC@samples$sample)) {
        cat("Given sample name/s in value not found in the QC slot of the object.")
      } else {
        object@samples <- rbind(object@samples, object@QC@samples[object@QC@samples$sample %in% value, ])
        object@samples <- object@samples[order(sample)]
        object@QC@samples <- object@QC@samples[!object@QC@samples$sample %in% value, ]
      }
    }
    return(object)
  }

  if (nameType == "replicates") {
    if (FALSE %in% unique(value %in% object@samples$replicate)) {
      cat("Given replicate name/s in value not found in the object.")
    } else {
      object@QC@samples <- rbind(object@QC@samples, object@samples[object@samples$replicate %in% value, ])
      object@samples <- object@samples[!(object@samples$replicate %in% value), ]
    }
  } else {
    if (FALSE %in% unique(value %in% object@samples$sample)) {
      cat("Given sample name/s in value not found in the object.")
    } else {
      object@QC@samples <- rbind(object@QC@samples, object@samples[object@samples$sample %in% value, ])
      object@samples <- object@samples[!(object@samples$sample %in% value), ]
    }
  }

  return(object)
})




### [ sub-setting samples -----------------------------------------------------------------------------------

#' @describeIn ntsData Subset on samples, using sample index or name.
#'
#' @param i The indice/s or name/s of the samples to keep in the \code{x} object.
#'
#' @export
#'
#' @importMethodsFrom MSnbase filterFile fileNames
#' @importMethodsFrom xcms filterFile
#' @importMethodsFrom patRoon analyses
#'
setMethod("[", c("ntsData", "ANY", "missing", "missing"), function(x, i, ...) {

  if (!missing(i)) {

    if (!is.character(i)) {
      sn <- x@samples$sample[i]
      sidx <- i
    } else {
      if (FALSE %in% (i %in% x@samples$sample)) {
        warning("Given sample name/s not found in the ntsData object.")
        return(x)
      }
      sn <- i
      sidx <- which(x@samples$sample %in% sn)
    }

    x@samples <- x@samples[sample %in% sn, ]

    x@metadata <- x@metadata[replicate %in% replicates(x), ]

    # if (length(analyses(x@patdata)) > 0) {

    #   x@patdata <- x@patdata[sidx]

    #   x@peaks <- x@peaks[x@peaks$sample %in% sn, ]

    #   if (nrow(x@features) > 0) x <- updateFeatureList(x)

    # }

    # #annotation, remove replicate groups without samples
    # if (nrow(x@annotation$comp) > 0) {
    #   rg <- unique(x@samples$group)
    #   x@annotation$comp <- x@annotation$comp[x@annotation$comp$group %in% rg, ]
    #   x@annotation$raw <- x@annotation$raw[names(x@annotation$raw) %in% rg]
    #   if (nrow(x@features) > 0) x@annotation$comp <- x@annotation$comp[x@annotation$comp$ID %in% x@features$ID, ]
    # }
  }

  return(x)

})




### EICs -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData get extracted ion chromatograms (EICs)
#' for specified \emph{m/z} and retention time (seconds) targets
#' in given samples. The arguments \code{mz}, \code{ppm}, \code{rt}
#' and \code{sec} are used to construct the targets.
#' See ?\link{makeTargets} for more information.
#'
#' @param samples A numeric or character vector with the indice/s or name/s
#' of samples from the \code{object}.
#' @param mz A numeric vector or data.table/data.frame to make targets.
#' See ?\link{makeTargets} for more information.
#' @param ppm A numeric vector of length one with the mass deviation, in ppm, to calculate ranges.
#' See ?\link{makeTargets} for more information.
#' @param rt A numeric vector or data.table/data.frame to make targets.
#' See ?\link{makeTargets} for more information.
#' @param sec A numeric vector of length one with the time deviation, in seconds, to calculate ranges.
#' See ?\link{makeTargets} for more information.
#'
#' @export
#'
setMethod("EICs", "ntsData", function(object,
                                      samples = NULL,
                                      mz = NULL, ppm = 20,
                                      rt = NULL, sec = 60) {

  eic <- extractEICs(
    object,
    samples,
    mz,
    ppm,
    rt,
    sec
  )

  return(eic)
})




### plotEICs ------------------------------------------------------------------------------------------------

#' @describeIn ntsData A method for plotting extracted ion chromatograms (EICs)
#' of data in an \linkS4class{ntsData} object.
#' The \code{colorBy} argument can be be \code{"samples"}, \code{replicates} or \code{targets}
#' (the default), for colouring by samples, replicates or EICs targets, respectively.
#' The \code{legendNames} is a character vector with the same length as targets for plotting and
#' can be used to lengend the plot. Note that, by setting \code{legendNames} the \code{colorBy}
#' is set to "targets".
#'
#' @param colorBy A length one character vector applicable to the respective method.
#' @param legendNames A character vector with the same length as the number of targets in the respective function.
#' @param interactive Logical value, set to \code{TRUE} to use
#' the \pkg{plotly} instead of \pkg{base}. The default is \code{FALSE}.
#'
#' @export
#'
setMethod("plotEICs", "ntsData", function(object,
                                          samples = NULL,
                                          mz = NULL, ppm = 20,
                                          rt = NULL, sec = 30,
                                          colorBy = "targets",
                                          legendNames = NULL,
                                          title = NULL,
                                          interactive = FALSE) {

  eic <- extractEICs(
    object,
    samples = samples,
    mz = mz,
    rt = rt,
    ppm = ppm,
    sec = sec
  )

  if (nrow(eic) < 1) return(cat("Data was not found for any of the targets!"))

  if (colorBy == "samples") {
    leg <- unique(eic$sample)
    varkey <- eic$sample
  } else if (colorBy == "replicates") {
    leg <- unique(eic[, .(sample, replicate)])
    leg <- leg$replicate
    varkey <- eic$replicate
  } else if (!is.null(legendNames) & length(legendNames) == length(sp)) {
    leg <- legendNames
    names(leg) <- unique(eic$id)
    varkey <- sapply(eic$id, function(x) leg[[x]])
  } else {
    leg <- unique(eic$id)
    varkey <- eic$id
  }

  eic[, var := varkey][]

  if (!interactive) {

    win.metafile()
    dev.control("enable")
    plotStaticEICs(
      eic,
      title
    )
    plot <- recordPlot()
    dev.off()

  } else {

    plot <- plotInteractiveEICs(eic, title, colorBy)

  }

  return(plot)
})




### plotEICs - data.table -----------------------------------------------------------------------------------

#' @describeIn ntsData A method for plotting extracted ion chromatograms (EICs)
#' of data in a \link[data.table]{data.table} object obtained with the \link{EICs} method.
#' The \code{colorBy} argument can be be \code{"samples"}, \code{replicates} or \code{targets}
#' (the default), for colouring by samples, replicates or EICs targets, respectively.
#' The \code{legendNames} is a character vector with the same length as targets for plotting and
#' can be used to lengend the plot. Note that, by setting \code{legendNames} the \code{colorBy}
#' is set to "targets". A subset of the targets in the object
#' can be plotted using the \code{targets} argument.
#'
#' @param targets A character vector with target names.
#'
#' @export
#'
setMethod("plotEICs", "data.table", function(object,
                                             samples = NULL,
                                             colorBy = "targets",
                                             legendNames = NULL,
                                             targets = NULL,
                                             title = NULL,
                                             interactive = FALSE) {

  eic <- object

  if (!is.null(samples)) {
    if (class(samples) == "numeric") samples <- unique(eic$sample)[samples]
    eic[sample %in% samples, ]
  }
  if (!is.null(targets)) eic[id %in% targets, ]

  if (nrow(eic) < 1) return(cat("Data was not found for any of the targets!"))

  if (colorBy == "samples") {
    leg <- unique(eic$sample)
    varkey <- eic$sample
  } else if (colorBy == "replicates") {
    leg <- unique(eic[, .(sample, replicate)])
    leg <- leg$replicate
    varkey <- eic$replicate
  } else if (!is.null(legendNames) & length(legendNames) == length(unique(eic$id))) {
    leg <- legendNames
    names(leg) <- unique(eic$id)
    varkey <- sapply(eic$id, function(x) leg[[x]])
  } else {
    leg <- unique(eic$id)
    varkey <- eic$id
  }

  eic[, var := varkey][]

  if (!interactive) {

    win.metafile()
    dev.control("enable")
    plotStaticEICs(
      eic,
      title
    )
    plot <- recordPlot()
    dev.off()

  } else {

    plot <- plotInteractiveEICs(eic, title, colorBy)

  }

  return(plot)
})




### TICs -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData Extract total ion chromatograms (TICs)
#' for samples in an \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("TICs", "ntsData", function(object, samples = NULL) {

  tic <- extractEICs(
    object,
    samples = samples,
    mz = NULL,
    rt = NULL
  )

  return(tic)
})




### plotTICs ------------------------------------------------------------------------------------------------

#' @describeIn ntsData Plots a total ion chromatogram (TIC)
#' from each sample in an \linkS4class{ntsData} object.
#' \code{colorBy} can be \code{"samples"} (the default)
#' or \code{replicates} for colouring by samples or replicates, respectively.
#'
#'
#' @export
#'
setMethod("plotTICs", "ntsData", function(object,
                                          samples = NULL,
                                          colorBy = "samples",
                                          title = NULL,
                                          interactive = FALSE) {

  ticplot <- plotEICs(
    object,
    samples = samples,
    mz = NULL,
    rt = NULL,
    colorBy = colorBy,
    title = title,
    interactive = interactive
  )

  return(ticplot)
})




### plotTICs - data.table -----------------------------------------------------------------------------------

#' @describeIn ntsData Plots a total ion chromatogram (TIC)
#' from each sample in a \link[data.table]{data.table} object as produced
#' by the \link{TICs} method.
#'
#' @export
#'
setMethod("plotTICs", "data.table", function(object,
                                             samples = NULL,
                                             colorBy = "samples",
                                             title = NULL,
                                             interactive = FALSE) {

  ticplot <- plotEICs(
    object,
    samples = samples,
    colorBy = colorBy,
    title = title,
    interactive = interactive
  )

  return(ticplot)
})




### XICs -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData get three dimentional (\emph{m/z}, time and intensity)
#' extracted ion chromatograms (XICs) for specified \emph{m/z} and retention time pair targets
#' in samples of an \linkS4class{ntsData} object. The arguments \code{mz}, \code{ppm}, \code{rt}
#' and \code{sec} are used to construct the targets.
#' See ?\link{makeTargets} for more information.
#'
#' @export
#'
setMethod("XICs", "ntsData", function(object,
                                      samples = NULL,
                                      mz = NULL, ppm = 20,
                                      rt = NULL, sec = 60) {

  xic <- extractXICs(
    object,
    samples,
    mz,
    ppm,
    rt,
    sec
  )

  return(xic)
})




### plotXICs ------------------------------------------------------------------------------------------------

#' @describeIn ntsData plots three dimentional (\emph{m/z}, time and intensity)
#' extracted ion chromatograms (XICs) for specified \emph{m/z} and retention time pair targets
#' in samples of an \linkS4class{ntsData} object. The arguments \code{mz}, \code{ppm}, \code{rt}
#' and \code{sec} are used to construct the targets. See ?\link{makeTargets} for more information.
#' When \code{plotTargetMark} is \code{TRUE} a target is plotted representing the deviations as defined
#' by the arguments \code{ppmMark} and \code{secMark} in ppm and seconds, respectively.
#' When ranges were given to build the XIC, exact \emph{m/z} and time targets can be specified with
#' the argument \code{targetsMark}. \code{targetsMark} should be a two column table named mz and rt with
#' exact \emph{m/z} and time targets. Note that the number of rows should be the same as the number of target
#' in the XIC. The number of rows to plot multiple targets can be defined by the \code{numberRows} argument.
#'
#' @param plotTargetMark Logical, set to \code{TRUE} to plot a target mark.
#' @param targetsMark A two columns \link[data.table]{data.table} or \link{data.frame} with
#' \emph{m/z} and time targets. The column must be named with "mz" and "rt" for
#' \emph{m/z} and time values, respectively.
#' @param ppmMark A numeric vector of length one to define the mass deviation, in ppm,
#' of the target mark.
#' @param secMark A numeric vector of length one to define the time deviation, in seconds,
#' of the target mark.
#' @param numberRows A numeric vector of length one to define
#' the number of rows to grid the plots.
#'
#' @export
#'
#' @importFrom data.table is.data.table
#'
setMethod("plotXICs", "ntsData", function(object,
                                          samples = NULL,
                                          mz = NULL, ppm = 20,
                                          rt = NULL, sec = 60,
                                          legendNames = NULL,
                                          plotTargetMark = TRUE,
                                          targetsMark = NULL,
                                          ppmMark = 5,
                                          secMark = 10,
                                          numberRows = 1) {

  xic <- extractXICs(
    object,
    samples,
    mz,
    ppm,
    rt,
    sec
  )

  if (nrow(xic) < 1) return(cat("Data was not found for any of the targets!"))

  #change id by legendNames
  ids <- unique(xic$id)
  if (!is.null(legendNames) & length(legendNames) == length(ids)) {
    names(legendNames) <- unique(xic$id)
    xic$id <- sapply(xic$id, function(x) legendNames[[x]])
  }

  if (plotTargetMark) {
    otherTargets <- FALSE
    if (!is.null(targetsMark)) {
      if ((is.data.table(targetsMark) | is.data.frame(targetsMark))) {
        if (nrow(targetsMark) == length(ids) & "mz" %in% colnames(targetsMark) & "rt" %in% colnames(targetsMark)) {
          tgmMZ <- as.numeric(targetsMark$mz)
          names(tgmMZ) <- unique(xic$id)
          tgmRT <- as.numeric(targetsMark$rt)
          names(tgmRT) <- unique(xic$id)
          xic[, mz_id := tgmMZ[xic$id]]
          xic[, rt_id := tgmRT[xic$id]]
          otherTargets <- TRUE
        }
      }
    }

    if (!otherTargets & class(xic$mz_id) == "character") {
      tgmMZ <- sapply(xic$mz_id, function(x) mean(as.numeric(stringr::str_split(x, "-", simplify = TRUE)[1, ])))
      tgmRT <- sapply(xic$rt_id, function(x) mean(as.numeric(stringr::str_split(x, "-", simplify = TRUE)[1, ])))
      xic[, mz_id := tgmMZ]
      xic[, rt_id := tgmRT]
    }
  }

  plot <- plotXIC_base(
    xic,
    plotTargetMark = plotTargetMark,
    ppmMark = ppmMark,
    secMark = secMark,
    numberRows = numberRows
  )

  return(plot)
})




### plotXICs - data.table -----------------------------------------------------------------------------------

#' @describeIn ntsData plots three dimentional (\emph{m/z}, time and intensity)
#' extracted ion chromatograms (XICs) for specified \emph{m/z} and retention time pair targets
#' in samples of a \link[data.table]{data.table} object as produced
#' by the \link{XICs} method. \code{samples} and \code{targets} can be used to filter the XIC table.
#' When \code{plotTargetMark} is \code{TRUE} a target is plotted representing the deviations as defined
#' by the arguments \code{ppmMark} and \code{secMark} in ppm and seconds, respectively.
#' When ranges were given to build the XIC, exact \emph{m/z} and time targets can be specified with
#' the argument \code{targetsMark}. \code{targetsMark} should be a two column table named mz and rt with
#' exact \emph{m/z} and time targets. Note that the number of rows should be the same as the number of target
#' in the XIC. The number of rows to plot multiple targets can be defined by the \code{numberRows} argument.
#'
#' @export
#'
#' @importFrom data.table is.data.table
#'
setMethod("plotXICs", "data.table", function(object,
                                          samples = NULL,
                                          targets = NULL,
                                          legendNames = NULL,
                                          plotTargetMark = TRUE,
                                          targetsMark = NULL,
                                          ppmMark = 5,
                                          secMark = 10,
                                          numberRows = 1) {

  xic <- object

  if (!is.null(samples)) {
    if (is.numeric(samples)) {
      samples <- unique(xic$sample)[samples]
    }
    xic <- xic[sample %in% samples, ]
  }

  if (!is.null(targets)) xic <- xic[id %in% targets, ]

  if (nrow(xic) < 1) return(cat("Data was not found for any of the targets!"))

  ids <- unique(xic$id)
  if (!is.null(legendNames) & length(legendNames) == length(ids)) {
    names(legendNames) <- unique(xic$id)
    xic$id <- sapply(xic$id, function(x) legendNames[[x]])
  }

  if (plotTargetMark) {
    otherTargets <- FALSE
    if (!is.null(targetsMark)) {
      if ((!is.data.table(targetsMark) | is.data.frame(targetsMark))) {
        if (nrow(targetsMark) == length(ids) & "mz" %in% colnames(targetsMark) & "rt" %in% colnames(targetsMark)) {
          tgmMZ <- targetsMark$mz
          names(tgmMZ) <- unique(xic$id)
          tgmRT <- targetsMark$rt
          names(tgmRT) <- unique(xic$id)
          xic[, mz_id := tgmMZ[xic$id]]
          xic[, rt_id := tgmRT[xic$id]]
          otherTargets <- TRUE
        }
      }
    }

    if (!otherTargets & class(xic$mz_id) == "character") {
      tgmMZ <- sapply(xic$mz_id, function(x) mean(as.numeric(stringr::str_split(x, "-", simplify = TRUE)[1, ])))
      tgmRT <- sapply(xic$rt_id, function(x) mean(as.numeric(stringr::str_split(x, "-", simplify = TRUE)[1, ])))
      xic[, mz_id := tgmMZ]
      xic[, rt_id := tgmRT]
    }
  }

  plot <- plotXIC_base(
    xic,
    plotTargetMark = plotTargetMark,
    ppmMark = ppmMark,
    secMark = secMark,
    numberRows = numberRows
  )

  return(plot)
})




### MS2s -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData get MS2 data for specified \emph{m/z} and retention time (seconds) targets
#' in samples of an \linkS4class{ntsData} object. The \code{clusteringUnit} defines the method used for clustering.
#' Possible values are \emph{euclidean} (the default) or \emph{distance}.
#' The \code{clusteringUnit} and \code{clusteringWindow} define
#' the mass deviation unit and deviation to cluster mass traces from different spectra, respectively.
#' For the \code{clusteringUnit}, possible values are \emph{mz} (the default) or \emph{ppm}.
#' The \code{minIntensityPre} and \code{minIntensityPost}
#' define the minimum intensity for mass traces before and after clustering, respectively.
#' Set \code{mergeCEs} to \code{TRUE} for merging spectra acquired with different collision energies.
#' The \code{mergeBy} argument is used to merge spectra by "samples" or "replicates".
#' When \code{NULL}, MS2 is given per target and per sample.
#'
#' @param clusteringMethod A character vector specifying the clustering unit.
#' @param clusteringUnit A character vector specifying the clustering unit.
#' @param clusteringWindow A length one numeric vector with the mass deviation for clustering.
#' @param minIntensityPre A length one numeric vector with the minimum intensity.
#' @param minIntensityPost A length one numeric vector with the minimum intensity.
#' @param mergeCEs Logical, set to TRUE to cluster different collision energies.
#' @param mergeBy A character string applicable to the respective method.
#'
#' @export
#'
setMethod("MS2s", "ntsData", function(object = NULL,
                                      samples = NULL,
                                      mz = NULL, ppm = 20,
                                      rt = NULL, sec = 60,
                                      clusteringMethod = "euclidean",
                                      clusteringUnit = "mz",
                                      clusteringWindow = 0.008,
                                      minIntensityPre = 250,
                                      minIntensityPost = 100,
                                      mergeCEs = FALSE,
                                      mergeBy = "samples") {

  level <- 2

  ms2 <- extractMSn(
    object,
    samples,
    level,
    mz, ppm,
    rt, sec,
    clusteringMethod,
    clusteringUnit,
    clusteringWindow,
    minIntensityPre,
    minIntensityPost,
    mergeCEs,
    mergeBy
  )

  return(ms2)
})




### plotMS2s -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData plots MS2 data for specified \emph{m/z} and retention time (seconds) targets
#' in samples of an \linkS4class{ntsData} object. The \code{clusteringUnit} defines the method used for clustering.
#' Possible values are \emph{euclidean} (the default) or \emph{distance}.
#' The \code{clusteringUnit} and \code{clusteringWindow} define
#' the mass deviation unit and deviation to cluster mass traces from different spectra, respectively.
#' For the \code{clusteringUnit}, possible values are \emph{mz} (the default) or \emph{ppm}.
#' The \code{minIntensityPre} and \code{minIntensityPost}
#' define the minimum intensity for mass traces before and after clustering, respectively.
#' Set \code{mergeCEs} to \code{TRUE} for merging spectra acquired with different collision energies.
#' The \code{mergeBy} argument is used to merge spectra by "samples" or "replicates".
#' When \code{NULL}, MS2 is given per target and per sample. The possible values for the
#' \code{colorBy} argument are "targets", "samples", "replicates" and "CEs" to color by
#' each target, sample, replicate or collision energy, respectively.
#'
#' @export
#'
setMethod("plotMS2s", "ntsData", function(object = NULL,
                                          samples = NULL,
                                          mz = NULL, ppm = 20,
                                          rt = NULL, sec = 60,
                                          clusteringMethod = "euclidean",
                                          clusteringUnit = "mz",
                                          clusteringWindow = 0.008,
                                          minIntensityPre = 250,
                                          minIntensityPost = 100,
                                          mergeCEs = FALSE,
                                          mergeBy = "samples",
                                          legendNames = NULL,
                                          title = NULL,
                                          colorBy = NULL,
                                          interactive = FALSE) {

  level <- 2

  ms2 <- extractMSn(
    object,
    samples,
    level,
    mz, ppm,
    rt, sec,
    clusteringMethod,
    clusteringUnit,
    clusteringWindow,
    minIntensityPre,
    minIntensityPost,
    mergeCEs,
    mergeBy
  )

  if (nrow(ms2) < 1) return(cat("Data was not found for any of the targets!"))

  if (colorBy == "samples" & "sample" %in% colnames(ms2)) {
    leg <- unique(ms2$sample)
    varkey <- ms2$sample
  } else if (colorBy == "replicates" & "replicate" %in% colnames(ms2)) {
    leg <- unique(ms2$replicate)
    varkey <- ms2$replicate
  } else if (colorBy == "CEs" & "CE" %in% colnames(ms2)) {
    leg <- unique(ms2$CE)
    varkey <- ms2$CE
  } else if (!is.null(legendNames) & length(legendNames) == length(unique(ms2$id))) {
    leg <- legendNames
    names(leg) <- unique(ms2$id)
    varkey <- sapply(ms2$id, function(x) leg[[x]])
  } else {
    leg <- unique(ms2$id)
    varkey <- ms2$id
  }

  ms2[, var := varkey]
  ms2$var <- factor(ms2$var, levels = unique(ms2$var), labels = unique(ms2$var))

  if (!interactive) {

    win.metafile()
    dev.control("enable")
    plotStaticMSn(
      ms2,
      title
    )
    plot <- recordPlot()
    dev.off()

  } else {

    plot <- plotInteractiveMSn(ms2, title)

  }

  return(plot)
})




### plotMS2s - data.table -----------------------------------------------------------------------------------

#' @describeIn ntsData plots MS2 data for specified \emph{m/z} and retention time (seconds) targets
#' in a \link[data.table]{data.table} as obtained by the \link{MS2s}. The targets in the
#' object can be filtered using the \code{targets} argument. Also, "samples" and "replicates"
#' can be filtered using the \code{samples} and \code{replicates} arguments, respectively.
#' Note that the column sample/replicate should be present.
#' The possible values for the \code{colorBy} argument are
#' "targets", "samples", "replicates" and "CEs" to color by
#' each target, sample, replicate or collision energy, respectively.
#'
#' @param replicates A numeric or character vector with the indice/s or name/s
#' of replicates from the \code{object}.
#'
#' @export
#'
setMethod("plotMS2s", "data.table", function(object = NULL,
                                             samples = NULL,
                                             replicates = NULL,
                                             targets = NULL,
                                             legendNames = NULL,
                                             title = NULL,
                                             colorBy = "targets",
                                             interactive = FALSE) {

  ms2 <- object

  if (!is.null(samples) & "sample" %in% colnames(ms2)) {
    if (class(samples) == "numeric") samples <- unique(ms2$sample)[samples]
    ms2[sample %in% samples, ]
  }

  if (!is.null(replicates) & "replicate" %in% colnames(ms2)) {
    if (class(replicates) == "numeric") replicates <- unique(ms2$replicate)[replicates]
    ms2[replicate %in% replicates, ]
  }

  if (!is.null(targets)) ms2[id %in% targets, ]

  if (nrow(ms2) < 1) return(cat("Data was not found for any of the targets!"))

  if (colorBy == "samples" & "sample" %in% colnames(ms2)) {
    leg <- unique(ms2$sample)
    varkey <- ms2$sample
  } else if (colorBy == "replicates" & "replicate" %in% colnames(ms2)) {
    leg <- unique(ms2$replicate)
    varkey <- ms2$replicate
  } else if (colorBy == "CEs" & "CE" %in% colnames(ms2)) {
    leg <- unique(ms2$CE)
    varkey <- ms2$CE
  } else if (!is.null(legendNames) & length(legendNames) == length(unique(ms2$id))) {
    leg <- legendNames
    names(leg) <- unique(ms2$id)
    varkey <- sapply(ms2$id, function(x) leg[[x]])
  } else {
    leg <- unique(ms2$id)
    varkey <- ms2$id
  }

  ms2[, var := varkey]
  ms2$var <- factor(ms2$var, levels = unique(ms2$var), labels = unique(ms2$var))

  if (!interactive) {

    win.metafile()
    dev.control("enable")
    plotStaticMSn(
      ms2,
      title
    )
    plot <- recordPlot()
    dev.off()

  } else {

    plot <- plotInteractiveMSn(ms2, title)

  }

  return(plot)
})





































### [ sub-setting features ----------------------------------------------------------------------------------

#' @describeIn ntsData Subset on samples, using sample index or name.
#'
#' @param x An \linkS4class{ntsData} object.
#' @param i The indice/s or name/s of the samples to keep in the \code{x} object.
#' @param j Ignored.
#' @param drop Ignored.
#' @param \dots Ignored.
#'
#' @export
#'
#' @importMethodsFrom MSnbase filterFile fileNames
#' @importMethodsFrom xcms filterFile
#' @importMethodsFrom patRoon analyses
#'
setMethod("[", c("ntsData", "ANY", "ANY", "missing"), function(x, i, j, ...) {

  if (!missing(i)) x <- x[i]

  if (!missing(j)) {

    if (nrow(x@features) == 0) {
      warning("There are no features in the ntsData object!")
      return(x)
    }

    if (!is.character(j)) {
      id <- x@features$ID[j]
      idn <- j
    } else {
      id <- j
      idn <- which(x@features$ID %in% j)
    }

    x@patdata <- x@patdata[, id]

    x@peaks <- x@peaks[x@peaks$feature %in% id, ]

    x@features <- x@features[idn, ]

    if (nrow(x@annotation$comp) != 0) {
      x@annotation$comp <- x@annotation$comp[x@annotation$comp$ID %in% id, ]
    }

  }

  return(x)

})




### peaks ---------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for chromatographic peaks.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples The indices or names of samples to keep.
#' @param ID The ID of peaks of interest.
#' @param mz Alternatively to \code{ID}, the \emph{m/z} of interest to find peaks.
#' Can be of length two, defining a mass range to find peaks.
#' @param ppm The mass deviation, in ppm, of a given \code{mz}.
#' @param rt The retention time to find peaks.
#' @param rtWindow The time deviation. Can be of length two, defining a time range.
#' @param rtUnit The unit of the time arguments. Possible values are "sec" and "min".
#'
#' @export
#'
#' @importFrom checkmate assertSubset
#' @importFrom dplyr filter between
#'
setMethod("peaks", "ntsData", function(object,
                                       samples = NULL, ID = NULL,
                                       mz = NULL, ppm = 20,
                                       rt = NULL, rtWindow = NULL,
                                       rtUnit = "sec") {

  if (nrow(object@peaks) == 0) return(object@peaks)

  if (!missing(samples)) object <- filterFileFaster(object, samples)

  if (missing(ID)) ID <- NULL

  if (missing(mz)) mz <- NULL

  rtr <- NULL

  if (!is.null(ID)) {
    pk <- object@peaks[object@peaks$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      if (missing(rt)) rt <- NULL
      if (missing(rtWindow)) rtWindow <- NULL
      if (missing(rtUnit)) rtUnit <- "sec"
      if (missing(ppm)) ppm <- 20
      assertSubset(rtUnit, c("sec", "min"))
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      pk <- dplyr::filter(object@peaks,
                          between(mz, mzr[1], mzr[2]),
                          between(rt, rtr[1], rtr[2]))
    } else {
      return(cat("One of ID or mz should be given."))
    }
  }

  return(pk)

})




### features -----

#' @describeIn ntsData Getter for features (i.e., grouped peaks).
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples The indices or names of samples to keep in the \code{object}.
#' @param ID The ID of features of interest.
#' @param mz Alternatively to \code{ID}, the \emph{m/z} of interest.
#' can be of length two, defining a mass range.
#' @param ppm The mass deviation, in ppm, of a given \code{mz}.
#' @param rt The retention time to find features.
#' @param rtWindow The time deviation. Can be of length two, defining a time range.
#' @param rtUnit The unit of the time arguments. Possible values are "sec" and "min".
#'
#' @export
#'
#' @importFrom checkmate assertSubset
#' @importFrom dplyr filter between
#'
setMethod("features", "ntsData", function(object,
                                          samples = NULL, ID = NULL,
                                          mz = NULL, ppm = 20,
                                          rt = NULL, rtWindow = NULL,
                                          rtUnit = "sec") {

  if (nrow(object@features) == 0) return(object@features)

  if (missing(ID)) ID <- NULL

  if (missing(mz)) mz <- NULL

  if (!missing(samples)) object <- filterFileFaster(object, samples)

  rtr <- NULL

  if (!is.null(ID)) {
    ft <- object@features[object@features$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      if (missing(rt)) rt <- NULL
      if (missing(rtWindow)) rtWindow <- NULL
      if (missing(rtUnit)) rtUnit <- "sec"
      if (missing(ppm)) ppm <- 20
      assertSubset(rtUnit, c("sec", "min"))
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      ft <- dplyr::filter(object@features,
                          between(mz, mzr[1], mzr[2]),
                          between(rt, rtr[1], rtr[2]))
    } else {
      return(cat("One of ID or mz should be given."))
    }
  }

  return(ft)

})




### components -----

#' @describeIn ntsData Getter for components (i.e., annotated features).
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples The indice/s or name/s of samples to keep in the \code{object}.
#' @param ID The ID of features of interest.
#' @param mz Alternatively to \code{ID}, the \emph{m/z} of interest.
#' can be of length two, defining a mass range.
#' @param ppm The mass deviation, in ppm, of a given \code{mz}.
#' @param rt The retention time to find features.
#' @param rtWindow The time deviation. Can be of length two, defining a time range.
#' @param rtUnit The unit of the time arguments. Possible values are "sec" and "min".
#' @param compNumber Alternatively, the component number to find features.
#' @param entireComponents Logical, set to \code{TRUE} to extract all features
#' from the represented components.
#' @param onlyAnnotated Logical, set to \code{TRUE} to extract only annotated features.
#' @param onlyRelated Logical, set to \code{TRUE} to extract only features that are related
#' to the features of interest.
#'
#' @export
#'
#' @importFrom checkmate assertSubset
#' @importFrom dplyr filter between
#' @importFrom stats na.omit
#'
setMethod("components", "ntsData", function(object,
                                            samples = NULL,
                                            ID = NULL,
                                            mz = NULL, ppm = 5,
                                            rt = NULL, rtWindow = 1, rtUnit = "sec",
                                            compNumber = NULL,
                                            entireComponents = TRUE,
                                            onlyAnnotated = FALSE,
                                            onlyRelated = TRUE) {

  if (missing(samples)) samples <- NULL

  if (missing(ID)) ID <- NULL

  if (missing(mz)) mz <- NULL

  if (missing(compNumber)) compNumber <- NULL

  comp <- object@annotation$comp

  if (nrow(comp) == 0) return(comp)

  rg <- unique(object@samples$group)

  if (!is.null(samples)) {
    if (is.character(samples)) {
      rg <- unique(object@samples$group[object@samples$sample %in% samples])
    } else {
      rg <- unique(object@samples$group[samples])
    }
  }

  comp <- comp[comp$group %in% rg, ]

  if (nrow(comp) == 0) return(comp)

  if (!is.null(ID)) {
    comp <- comp[comp$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      if (missing(rt)) rt <- NULL
      if (missing(rtWindow)) rtWindow <- NULL
      if (missing(rtUnit)) rtUnit <- "sec"
      if (missing(ppm)) ppm <- 20
      assertSubset(rtUnit, c("sec", "min"))
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      comp <- dplyr::filter(comp,
                            between(mz, mzr[1], mzr[2]),
                            between(rt, rtr[1], rtr[2]))
    } else {
      if (!is.null(compNumber)) {
        comp <- comp[comp$comp %in% compNumber, ]
      }
    }
  }

  if (nrow(comp) == 0) return(comp)

  comp2 <- comp

  if (!missing(entireComponents)) {
    if (entireComponents) {
      comp2 <- object@annotation$comp[object@annotation$comp$comp %in% comp2$comp |
                                        object@annotation$comp$isogroup %in% comp2$isogroup, ]
      comp2 <- comp2[comp2$group %in% rg, ]
    }
  }

  if (!missing(onlyAnnotated)) {
    if (onlyAnnotated) comp2 <- comp2[!is.na(comp2$isoclass) | (comp2$nAdducts > 0), ]
  }

  if (!missing(onlyRelated)) {
    if (onlyRelated) {
      isos <- comp2$ID[comp2$Mion %in% stats::na.omit(comp$Mion)]
      comp2 <- comp2[comp2$ID %in% unique(isos), ]
    }
  }

  return(comp2)

})




### show -----

#' @describeIn ntsData Informative printing of the \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#' @importFrom dplyr count
setMethod("show", "ntsData", function(object) {

  st <- object@samples[, c("sample", "replicate", "blank", "polarity")]
  if (length(object@peaks) > 0) {
    st$peaks <- count(object@peaks, sample)$n
  }

  cat(
    is(object)[[1]], "\n",
    "  Project: ", object@title, "\n",
    "  Date:  ", as.character(object@date), "\n",
    "  Path:  ", object@path, "\n",
    #"  Samples:  ", "\n", "\n", sep = ""
    "\n", sep = ""
  )

  if (nrow(st) > 0) rownames(st) <- seq_len(nrow(st))
  print(st)

  cat("\n")

  mtd <- colnames(object@metadata)
  mtd <- mtd[mtd != "replicate"]
  if (length(mtd) >= 1) { cat("  Metadata Variables:  ", mtd, "\n", sep = " ") }
  if (length(mtd) < 1) { cat("  Metadata Variables:  ", "empty", "\n", sep = " ") }

  cat("\n")

  cat("  QC samples:  ",
    ifelse(nrow(object@QC@samples) < 1, "empty", paste(object@QC@samples$sample, collapse = ", ")),
    "\n", " QC results:  ",
    ifelse(length(object@QC@results) < 1, "empty", paste(nrow(object@QC@results), " compounds evaluated")),
    "\n", set = ""
  )

  cat("\n")

  cat(
    "  Parameters:  ", "\n",
    "    ", "Peak Picking: ", ifelse(!is.na(object@parameters@peakPicking@algorithm),
                                      object@parameters@peakPicking@algorithm,
                                      "n.a."), "\n",

    "    ", "Grouping: ", ifelse(!is.na(object@parameters@peakGrouping@algorithm),
                                  object@parameters@peakGrouping@algorithm,
                                  "n.a."), "\n",

    "    ", "Fill Missing: ", ifelse(!is.na(object@parameters@fillMissing@algorithm),
                                      object@parameters@fillMissing@algorithm,
                                      "n.a."), "\n",

    "    ", "Annotation: ", ifelse(!is.na(object@parameters@annotation@algorithm),
                                    object@parameters@annotation@algorithm,
                                    "n.a."), "\n",

    "    ", "MS2: ", ifelse(length(object@parameters@MS2) > 0,
                            class(object@parameters@MS2),
                            "empty"), "\n",
    sep = ""
  )

  cat("\n")

  cat(
    # "  MSnExp spectra:  ", length(object@MSnExp), "\n",
    # "  MS level(s):  ", ifelse(length(object@MSnExp) > 0,
    #                             paste(sort(unique(msLevel(object@MSnExp))), collapse = ", "),
    #                             "-"), "\n",
    "  patRoon-class:  ", ifelse(length(object@pat) > 0,
                                  class(object@pat), "empty"), "\n",
    sep = ""
  )

  # if (nrow(object@annotation$comp) > 0) {
  #   cat("  Annotation:  ", "Yes", "\n", sep = "")
  # }

  if (length(object@features) > 0) {
    cat("\n")
    cat("  Number of features:  ", nrow(object@features), "\n", sep = "")
  }

  if (length(object@features) > 0) {
    cat("  Filtered features:  ", nrow(object@removed), "\n", sep = "")
  }

  if (length(object@workflows) > 0) {
    cat("\n")
    cat("  Workflow Objects:\n", paste(names(object@workflows), "\n", sep = ""))
  }

})
