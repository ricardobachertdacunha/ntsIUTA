

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
  object@metadata$replicate <- object@samples$replicate
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

#' @describeIn ntsData Getter for the polarity of the \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("polarity", "ntsData", function(object) {

  return(object@polarity)
})

#' @describeIn ntsData Setter for the polarity mode of the samples (i.e., files).
#' The \code{value} is a character vector with length one with either \emph{positive} or \emph{negative} strings.
#'
#' @export
#' 
#' @importFrom checkmate testChoice
#'
setMethod("polarity<-", "ntsData", function(object, value) {

  if (testChoice(value, c("positive", "negative"))) object@polarity <- value

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
      if (!"replicate" %in% varname) {
        varname <- c("replicate", varname)
      }
      if (!"sample" %in% varname) {
        varname <- c("sample", varname)
      }
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
        object@metadata <- rbind(object@metadata, object@QC@samples[object@QC@samples$replicate %in% value, .(sample, replicate)], fill = TRUE)
        object@samples <- object@samples[order(sample)]
        object@metadata <- object@metadata[order(sample)]

        logVec <- object@QC@samples$replicate %in% value
        object@scans <- c(object@scans, object@QC@scans[logVec])
        object@scans <- object@scans[order(names(object@scans))]

        object@QC@samples <- object@QC@samples[!object@QC@samples$replicate %in% value, ]
        object@QC@scans <- object@QC@scans[!logVec]
      }
    } else {
      if (FALSE %in% unique(value %in% object@QC@samples$sample)) {
        cat("Given sample name/s in value not found in the QC slot of the object.")
      } else {
        object@samples <- rbind(object@samples, object@QC@samples[object@QC@samples$sample %in% value, ])
        object@metadata <- rbind(object@metadata, object@QC@samples[object@QC@samples$sample %in% value, .(sample, replicate)], fill = TRUE)
        object@samples <- object@samples[order(sample)]
        object@metadata <- object@metadata[order(sample)]

        logVec <- object@QC@samples$sample %in% value
        object@scans <- c(object@scans, object@QC@scans[logVec])
        object@scans <- object@scans[order(names(object@scans))]

        object@QC@samples <- object@QC@samples[!object@QC@samples$sample %in% value, ]
        object@QC@scans <- object@QC@scans[!logVec]
      }
    }
    return(object)
  }

  if (nameType == "replicates") {
    if (FALSE %in% unique(value %in% object@samples$replicate)) {
      cat("Given replicate name/s in value not found in the object.")
    } else {
      object@QC@samples <- rbind(object@QC@samples, object@samples[object@samples$replicate %in% value, ])

      logVec <- object@samples$replicate %in% value
      object@QC@scans <- c(object@QC@scans, object@scans[logVec])

      object@samples <- object@samples[!(object@samples$replicate %in% value), ]
      object@metadata <- object@metadata[!(object@metadata$replicate %in% value), ]
      object@scans <- object@scans[!logVec]
    }
  } else {
    if (FALSE %in% unique(value %in% object@samples$sample)) {
      cat("Given sample name/s in value not found in the object.")
    } else {
      object@QC@samples <- rbind(object@QC@samples, object@samples[object@samples$sample %in% value, ])

      logVec <- object@samples$sample %in% value
      object@QC@scans <- c(object@QC@scans, object@scans[logVec])

      object@samples <- object@samples[!(object@samples$sample %in% value), ]
      object@metadata <- object@metadata[!(object@metadata$sample %in% value), ]
      object@scans <- object@scans[!logVec]
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
#' @importMethodsFrom patRoon analyses
#'
setMethod("[", c("ntsData", "ANY", "missing", "missing"), function(x, i, ...) {

  if (!missing(i)) {

    if (!is.character(i)) {
      sname <- x@samples$sample[i]
      sidx <- i
    } else {
      if (FALSE %in% (i %in% x@samples$sample)) {
        warning("Given sample name/s not found in the ntsData object.")
        return(x)
      }
      sname <- i
      sidx <- which(x@samples$sample %in% sname)
    }

    x@samples <- x@samples[sample %in% sname, ]

    x@metadata <- x@metadata[sample %in% sname, ]

    x@scans <- x@scans[sname]

    if (length(analyses(x@pat)) > 0) {

      x@pat <- x@pat[sidx]

      if (nrow(x@peaks) > 0) x@peaks <- x@peaks[sample %in% sname, ]

      if (nrow(x@features) > 0) x <- buildFeatureTable(x)
      # TODO update features ranges based on remaining peaks, not creating a new feature table to avoid losing info

    }
  }

  return(x)
})



### hasAdjustedRetentionTime -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData checks if the \linkS4class{ntsData} has adjusted
#' retention time derived from alighment of peaks across samples.
#'
#' @export
#'
setMethod("hasAdjustedRetentionTime", "ntsData", function(object) {

  hasAdj <- object@scans
  hasAdj <- lapply(hasAdj, function(x) "adjustedRetentionTime" %in% colnames(x))
  hasAdj <- unlist(hasAdj)

  return(all(hasAdj))
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

    # win.metafile()
    # dev.control("enable")
    # plotStaticEICs(
    #   eic,
    #   title
    # )
    # plot <- recordPlot()
    # dev.off()

    return(
      plotStaticEICs(
        eic,
        title
      )
    )

  } else {

    plot <- plotInteractiveEICs(eic, title, colorBy)

    return(plot)
  }
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

    # win.metafile()
    # dev.control("enable")
    # plotStaticEICs(
    #   eic,
    #   title
    # )
    # plot <- recordPlot()
    # dev.off()

    return(
      plotStaticEICs(
        eic,
        title
      )
    )

  } else {

    plot <- plotInteractiveEICs(eic, title, colorBy)

    return(plot)
  }
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

  return(
    plotEICs(
      object,
      samples = samples,
      mz = NULL,
      rt = NULL,
      colorBy = colorBy,
      title = title,
      interactive = interactive
    )
  )
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

  return(
    plotEICs(
      object,
      samples = samples,
      colorBy = colorBy,
      title = title,
      interactive = interactive
    )
  )
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
#' The mass (in Da) and time (in seconds) isolation windows to screen for the respective precursors
#' are defined with the arguments \code{isolationMassWindow} and \code{isolationTimeWindow}, respectively.
#' The \code{clusteringUnit} and \code{clusteringWindow} define
#' the mass deviation unit and deviation to cluster mass traces from different spectra, respectively.
#' For the \code{clusteringUnit}, possible values are \emph{mz} (the default) or \emph{ppm}.
#' The \code{minIntensityPre} and \code{minIntensityPost}
#' define the minimum intensity for mass traces before and after clustering, respectively.
#' Set \code{mergeCEs} to \code{TRUE} for merging spectra acquired with different collision energies.
#' The \code{mergeBy} argument is used to merge spectra by "samples" or "replicates".
#' When \code{NULL}, MS2 is given per target and per sample.
#' 
#' @param isolationTimeWindow A character vector of length one with the time isolation window, in seconds.
#' @param isolationMassWindow A character vector of length one with the mass isolation window, in Da.
#' @param clusteringMethod A character vector specifying the clustering unit.
#' @param clusteringUnit A character vector specifying the clustering unit.
#' @param clusteringWindow A length one numeric vector with the mass deviation for clustering.
#' @param minIntensityPre A length one numeric vector with the minimum intensity.
#' @param minIntensityPost A length one numeric vector with the minimum intensity.
#' @param asPatRoon Logical, set to \code{TRUE} for return a pkg{patRoon} class object.
#' @param mergeCEs Logical, set to TRUE to cluster different collision energies.
#' @param mergeBy A character string applicable to the respective method.
#'
#' @export
#'
setMethod("MS2s", "ntsData", function(object = NULL,
                                      samples = NULL,
                                      mz = NULL, ppm = 20,
                                      rt = NULL, sec = 60,
                                      isolationTimeWindow = 10,
                                      isolationMassWindow = 1.3,
                                      clusteringMethod = "distance",
                                      clusteringUnit = "mz",
                                      clusteringWindow = 0.005,
                                      minIntensityPre = 200,
                                      minIntensityPost = 200,
                                      asPatRoon = FALSE,
                                      mergeCEs = FALSE,
                                      mergeBy = "samples") {

  level <- 2

  ms2 <- extractMSn(
    object,
    samples,
    level,
    mz, ppm,
    rt, sec,
    isolationTimeWindow,
    isolationMassWindow,
    clusteringMethod,
    clusteringUnit,
    clusteringWindow,
    minIntensityPre,
    minIntensityPost,
    asPatRoon,
    mergeCEs,
    mergeBy
  )

  return(ms2)
})




### plotMS2s -----------------------------------------------------------------------------------------------------

#' @describeIn ntsData plots MS2 data for specified \emph{m/z} and retention time (seconds) targets
#' in samples of an \linkS4class{ntsData} object. The \code{clusteringUnit} defines the method used for clustering.
#' Possible values are \emph{euclidean} (the default) or \emph{distance}.
#' The mass (in Da) and time (in seconds) isolation windows to screen for the respective precursors
#' are defined with the arguments \code{isolationMassWindow} and \code{isolationTimeWindow}, respectively.
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
                                          isolationTimeWindow = 10,
                                          isolationMassWindow = 1.3,
                                          clusteringMethod = "distance",
                                          clusteringUnit = "mz",
                                          clusteringWindow = 0.005,
                                          minIntensityPre = 200,
                                          minIntensityPost = 200,
                                          mergeCEs = FALSE,
                                          mergeBy = "samples",
                                          legendNames = NULL,
                                          title = NULL,
                                          colorBy = NULL,
                                          interactive = FALSE) {

  level <- 2

  asPatRoon <- FALSE

  ms2 <- extractMSn(
    object,
    samples,
    level,
    mz, ppm,
    rt, sec,
    isolationTimeWindow,
    isolationMassWindow,
    clusteringMethod,
    clusteringUnit,
    clusteringWindow,
    minIntensityPre,
    minIntensityPost,
    asPatRoon,
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




### peaks ---------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for chromatographic peaks.
#' The arguments \code{targets} and \code{mz}/\code{rt} can be used
#' to select specific peaks. The \emph{id} of peaks and/or features can be given in the \code{targets}
#' argument to select the respective peaks. Also, samples can be selected using the
#' \code{samples} argument.
#'
#' @export
#'
#' @importFrom dplyr between
#'
setMethod("peaks", "ntsData", function(object,
                                       samples = NULL,
                                       targets = NULL,
                                       mz = NULL, ppm = 20,
                                       rt = NULL, sec = 60) {

  pks <- object@peaks

  if (!is.null(samples)) {
    if (!class(samples) == "character") {
      samples <- samples(object)[samples]
    }
    pks <- pks[sample %in% samples, ]
  }

  if (!is.null(targets) & "feature" %in% colnames(pks)) {
    pks <- pks[id %in% targets | feature %in% targets, ]
    return(pks)
  } else if (!is.null(targets)) {
    pks <- pks[id %in% targets, ]
    return(pks)
  }

  if (!is.null(mz)) {
    targets <- makeTargets(mz, rt, ppm, sec)

    sel <- rep(FALSE, nrow(pks))
    for (i in seq_len(nrow(targets))) {
      sel[between(pks$mz, targets$mzmin[i], targets$mzmax[i]) &
          between(pks$rt, targets$rtmin[i], targets$rtmax[i])] <- TRUE
    }

    return(pks[sel])
  }

  return(pks)
})




### plotPeaks ------------------------------------------------------------------------------------------------

#' @describeIn ntsData A method for plotting chromatographic peaks
#' in an \linkS4class{ntsData} object.
#' The \code{colorBy} argument can be be \code{"samples"}, \code{replicates} or \code{targets}
#' (the default), for colouring by samples, replicates or peak targets, respectively.
#' The \code{legendNames} is a character vector with the same length as targets for plotting and
#' can be used to lengend the plot. Note that, by setting \code{legendNames} the \code{colorBy}
#' is set to "targets" automatically.
#'
#'
#' @export
#'
#' @importFrom data.table rbindlist
#'
setMethod("plotPeaks", "ntsData", function(object,
                                           samples = NULL,
                                           targets = NULL,
                                           mz = NULL, ppm = 20,
                                           rt = NULL, sec = 30,
                                           colorBy = "targets",
                                           legendNames = NULL,
                                           title = NULL,
                                           interactive = FALSE) {

  pks <- peaks(
    object,
    samples,
    targets,
    mz, ppm,
    rt, sec
  )

  eic <- lapply(split(pks, pks$sample), function(x, object) {
    return(
      extractEICs(
        object,
        samples = unique(x$sample),
        mz = x
      )
    )
  }, object = object)

  eic <- rbindlist(eic)

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
    varkey <- sapply(eic$id, function(x) leg[x])
  } else {
    leg <- paste0(pks$id, " - ", round(pks$mz, digits = 4), "/", round(pks$rt, digits = 0))
    names(leg) <- unique(eic$id)
    varkey <- sapply(eic$id, function(x) leg[[x]])
  }

  eic[, var := varkey][]

  if (!interactive) {

    return(
      plotPeaksStatic(
        eic,
        pks,
        title
      )
    )

  } else {

    plot <- plotPeaksInteractive(eic, pks, title, colorBy)

    return(plot)
  }
})




### mapPeaks ------------------------------------------------------------------------------------------------

#' @describeIn ntsData A method for mapping peaks mass and time space
#' in an \linkS4class{ntsData} object.
#' The \code{colorBy} argument can be be \code{"samples"}, \code{replicates} or \code{targets}
#' (the default), for colouring by samples, replicates or peak targets, respectively.
#' The \code{legendNames} is a character vector with the same length as targets for plotting and
#' can be used to lengend the plot. Note that, by setting \code{legendNames} the \code{colorBy}
#' is set to "targets" automatically.
#'
#' @param xlim A length one or two numeric vector for setting the \emph{x} limits of a plot.
#' @param ylim A length one or two numeric vector for setting the \emph{y} limits of a plot.
#'
#' @export
#'
setMethod("mapPeaks", "ntsData", function(object,
                                          samples = NULL,
                                          targets = NULL,
                                          mz = NULL, ppm = 20,
                                          rt = NULL, sec = 30,
                                          colorBy = "targets",
                                          legendNames = NULL,
                                          xlim = 30,
                                          ylim = 0.05,
                                          title = NULL) {

  pks <- peaks(
    object,
    samples,
    targets,
    mz, ppm,
    rt, sec
  )

  if (nrow(pks) < 1) return(cat("Requested peaks were not found!"))

  if (colorBy == "samples") {
    leg <- unique(pks$sample)
    varkey <- pks$sample
  } else if (colorBy == "replicates") {
    leg <- unique(pks[, .(sample, replicate)])
    leg <- leg$replicate
    varkey <- pks$replicate
  } else if (!is.null(legendNames) & length(legendNames) == length(unique(pks$id))) {
    leg <- legendNames
    names(leg) <- unique(pks$id)
    varkey <- sapply(pks$id, function(x) leg[x])
  } else {
    leg <- paste0(pks$id, " - ", round(pks$mz, digits = 4), "/", round(pks$rt, digits = 0))
    names(leg) <- pks$id
    varkey <- sapply(pks$id, function(x) leg[names(leg) == x])
  }

  pks[, var := varkey][]

  plot <- mapPeaksInteractive(pks, xlim, ylim, title)

  return(plot)
})




### [ sub-setting features ----------------------------------------------------------------------------------

#' @describeIn ntsData Subset on samples, using sample index or name.
#'
#' @param j The indice/s or \emph{id}/s for subsettting on features.
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

    if (!is.character(j)) j <- x@features$id[j]

    x@pat <- x@pat[, j]

    x@peaks <- x@peaks[feature %in% j, ]

    x@features <- x@features[id %in% j, ]
  }

  return(x)
})




### features ------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for features (i.e., grouped peaks).
#'
#' @export
#'
#' @importFrom checkmate assertSubset
#' @importFrom dplyr filter between
#'
setMethod("features", "ntsData", function(object,
                                          samples = NULL,
                                          targets = NULL,
                                          mz = NULL, ppm = 20,
                                          rt = NULL, sec = 60) {

  feats <- object@features

  if (!is.null(targets)) return(feats[id %in% targets, ])

  if (!is.null(mz)) {
    targets <- makeTargets(mz, rt, ppm, sec)

    sel <- rep(FALSE, nrow(feats))
    for (i in seq_len(nrow(targets))) {
      sel[between(feats$mz, targets$mzmin[i], targets$mzmax[i]) &
          between(feats$rt, targets$rtmin[i], targets$rtmax[i])] <- TRUE
    }

    return(feats[sel])
  }

  return(feats)
})




### plotFeatures --------------------------------------------------------------------------------------------

#' @describeIn ntsData A method for plotting peaks from given features
#' in an \linkS4class{ntsData} object.
#' The \code{colorBy} argument can be be \code{"samples"}, \code{replicates} or \code{targets}
#' (the default), for colouring by samples, replicates or peak targets, respectively.
#' The \code{legendNames} is a character vector with the same length as targets for plotting and
#' can be used to lengend the plot. Note that, by setting \code{legendNames} the \code{colorBy}
#' is set to "targets" automatically.
#'
#'
#' @export
#'
#' @importFrom data.table rbindlist
#'
setMethod("plotFeatures", "ntsData", function(object,
                                              samples = NULL,
                                              targets = NULL,
                                              mz = NULL, ppm = 20,
                                              rt = NULL, sec = 30,
                                              colorBy = "targets",
                                              legendNames = NULL,
                                              title = NULL,
                                              interactive = FALSE) {

  feats <- features(
    object,
    samples,
    targets,
    mz,
    ppm,
    rt,
    sec
  )

  pks <- peaks(
    object,
    samples,
    targets = feats$id
  )

  if (!is.null(legendNames) & length(legendNames) == length(unique(pks$feature))) {
    names(legendNames) <- unique(pks$feature)
    pks$feature <- sapply(pks$feature, function(x) legendNames[x])
  } else {
    legendNames <- pks$feature
  }

  eic <- lapply(split(pks, pks$sample), function(x, object) {
    return(
      extractEICs(
        object,
        samples = unique(x$sample),
        mz = x
      )
    )
  }, object = object)

  eic <- rbindlist(eic)

  if (hasAdjustedRetentionTime(object)) {
    spls <- unique(eic$sample)
    for (i in spls) {
      eic[sample == i, rt := sapply(rt, function(x, object, i) {
        object@scans[[i]][retentionTime == x, adjustedRetentionTime]
      }, object = object, i = i)]
    }
  }

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
    varkey <- sapply(eic$id, function(x) leg[x])
  } else {
    leg <- paste0(pks$id, " - ", round(pks$mz, digits = 4), "/", round(pks$rt, digits = 0))
    names(leg) <- unique(eic$id)
    varkey <- sapply(eic$id, function(x) leg[[x]])
  }

  eic[, var := varkey][]

  for (f_tar in unique(pks$feature)) {
    pks[feature %in% f_tar, rt := feats[id %in% f_tar, rt]]
  }

  if (!interactive) {

    return(
      plotPeaksStatic(
        eic,
        pks,
        title
      )
    )

  } else {

    return(
      plotPeaksInteractive(
        eic,
        pks,
        title,
        colorBy
      )
    )
  }

})




### show ----------------------------------------------------------------------------------------------------

#' @describeIn ntsData Informative printing of the \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#' @importFrom dplyr count
setMethod("show", "ntsData", function(object) {

  st <- object@samples[, .(sample, replicate, blank)]
  if (length(object@peaks) > 0) {
    st$peaks <- count(object@peaks, sample)$n
  }

  cat(
    is(object)[[1]], "\n",
    "  Project: ", object@title, "\n",
    "  Date:  ", as.character(object@date), "\n",
    "  Polarity:  ", object@polarity, "\n",
    "  Path:  ", object@path, "\n",
    #"  Samples:  ", "\n", "\n", sep = ""
    "\n", sep = ""
  )

  if (nrow(st) > 0) rownames(st) <- seq_len(nrow(st))
  print(st)

  cat("\n")

  mtd <- colnames(object@metadata)
  mtd <- mtd[!mtd %in% c("sample", "replicate")]
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
    "    ", "Picking:    ", ifelse(!is.na(object@parameters@picking@algorithm),
                                      object@parameters@picking@algorithm,
                                      "n.a."), "\n",

    "    ", "Grouping:   ", ifelse(!is.na(object@parameters@grouping@algorithm),
                                  object@parameters@grouping@algorithm,
                                  "n.a."), "\n",

    "    ", "Filling:    ", ifelse(!is.na(object@parameters@filling@algorithm),
                                      object@parameters@filling@algorithm,
                                      "n.a."), "\n",

    "    ", "Annotation: ", ifelse(!is.na(object@parameters@annotation@algorithm),
                                    object@parameters@annotation@algorithm,
                                    "n.a."), "\n",

    "    ", "Fragments:  ", ifelse(!is.na(object@parameters@fragments@algorithm),
                                   object@parameters@fragments@algorithm,
                                   "n.a."), "\n",
    sep = ""
  )

  cat("\n")

  cat(
    "  patRoon-class:  ", ifelse(length(object@pat) > 0,
                                 class(object@pat), "empty"), "\n",
    sep = ""
  )

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
