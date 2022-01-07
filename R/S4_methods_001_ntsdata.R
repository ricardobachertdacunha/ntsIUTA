

### projectInfo ---------------------------------------------------------------------------------------------

#' @describeIn ntsData setter and getter for project information.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param title The project title.
#' @param description The project description.
#' @param date The project date.
#'
#' @return An \linkS4class{ntsData} object with updated project information.
#'
#' @note When no arguments are given in the call for \code{projectInfo} a list
#' with the title, description and date is returned instead of the \code{object}.
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
#' @param object An \linkS4class{ntsData} object.
#'
#' @return A character vector with project path.
#'
#' @export
#' 
setMethod("path", "ntsData", function(object) object@path)




### samples -------------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for sample names.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @return A character vector with sample names.
#'
#' @export
#' 
setMethod("samples", "ntsData", function(object) object@samples$sample)





### replicates ----------------------------------------------------------------------------------------------

#' @describeIn ntsData Getter for sample replicate names.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @return A character vector with sample replicate names.
#'
#' @export
#'
setMethod("replicates", "ntsData", function(object) object@samples$replicate)

#' @describeIn ntsData Setter for sample replicate names.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector to assign sample replicate groups.
#'
#' @return Adds the value vector with sample replicate names
#' to the replicate column of the samples data table in the object.
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
#' @param object An \linkS4class{ntsData} object.
#'
#' @return A character vector with blank sample replicate names.
#'
#' @export
#'
setMethod("blanks", "ntsData", function(object) object@samples$blank)

#' @describeIn ntsData Setter for blank replicate groups.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector with the name/s of the blank sample replicate group/s.
#' If more than one, the length of \code{value} should be equal to
#' the number of samples.
#'
#' @return Adds the value vector with blank sample replicate names
#' to the blank column of the samples data table in the object.
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
#' @param object An \linkS4class{ntsData} object.
#'
#' @return A vector with acquisition method names used for each sample
#' or when only one method is present, the single acquisition method name.
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
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector with the name/s of the acquisition methods used for each file.
#' If only one name is given, it is used for all samples. If more than one, the \code{value} should
#' have the same length as the number of samples in the \code{object}.
#'
#' @return Adds the \code{value} vector with acquisition method names
#' to the method column of the samples data table in the object.
#'
#' @export
#'
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
#'
#' @param object An \linkS4class{ntsData} object.
#' @param groupby Either \code{samples} (the default) or \code{replicates} as character string
#' to return the polarites either for each sample or each replicate.
#'
#' @return A vector with the polarity mode of each sample (i.e., file).
#'
#' @export
#'
setMethod("polarity", "ntsData", function(object, groupby) {

  if (missing(groupby)) {groupby <- "samples"}

  if (groupby == "samples") {
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
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector with either \emph{positive} or \emph{negative} strings
#' with the same length as the number of samples in the set. If only one polarity mode is used for all the samples,
#' the \code{value} can be of length 1 and either \emph{positive} or \emph{negative}.
#' The polarity mode is used for all the samples.
#'
#' @return Adds the \code{value} vector with polarity modes
#' to the polarity column of the samples data table in the object.
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

#' @describeIn ntsData Getter for metadata.
#'
#' @param x An \linkS4class{ntsData} object.
#' @param varname The name of the variables to extract.
#' 
#' @return A data table with the requested metadata.
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

#' @describeIn ntsData Getter for QC samples data table.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("QC", "ntsData", function(object) object@QC@samples)

#' @describeIn ntsData Setter for QC samples or sample replicates.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector with the names of the samples or sample replicate name/s
#' to be used for QC.
#' @param nametype A character string specifying if sample or replicate names should be used.
#' Posiible values are \emph{samples} or \emph{replicates} (the default).
#' @param remove A logical value. When \code{TRUE}, specified sample or replicate names
#' are moved from the QC slot to the samples slot of the \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("QC<-", "ntsData", function(object, value, remove = FALSE, nametype = "replicates") {

  if (missing(nametype)) nametype <- "replicates"

  if (!missing(remove) & remove) {
    if (nametype == "replicates") {
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

  if (nametype == "replicates") {
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
setMethod("[", c("ntsData", "ANY", "missing", "missing"), function(x, i, ...) {

  if (!missing(i)) {

    if (!is.character(i)) {
      sn <- x@samples$sample[i]
      sidx <- which(x@samples$sample %in% sn)
    } else {
      if (FALSE %in% (i %in% x@samples$sample)) {
        warning("Given sample name/s not found in the ntsData object.")
        return(x)
      }
        sn <- i
        sidx <- which(x@samples$sample %in% sn)
    }

    x@samples <- x@samples[x@samples$sample %in% sn,, drop = FALSE]

    x@metadata <- x@metadata[x@metadata$replicate %in% x@samples$replicate,, drop = FALSE]

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
    "  Samples:  ", "\n", "\n", sep = ""
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
