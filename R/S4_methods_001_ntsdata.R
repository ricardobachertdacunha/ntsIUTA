

### show -----

#' @describeIn ntsData Informative printing.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#' @importFrom dplyr count
setMethod("show", "ntsData", function(object) {

  st <- object@samples[, c("sample", "group", "blank")]
  if (nrow(object@peaks) > 0) {
    st$peaks <- count(object@peaks, sample)$n
  }

  cat(
    is(object)[[1]], "\n",
    "  Project: ", object@title, "\n",
    "  Date:  ", as.character(object@date), "\n",
    "  Path:  ", object@path, "\n",
    "  Polarity:  ", object@polarity, "\n",
    "  Samples:  ", "\n", "\n", sep = ""
  )

  if (nrow(st) > 0) rownames(st) <- seq_len(nrow(st))
  print(st)

  cat("\n")

  mtd <- colnames(object@metadata)
  mtd <- mtd[mtd != "sample"]
  cat("  Metadata Variables:  ", mtd, "\n", sep = " ")

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
    "  MSnExp spectra:  ", length(object@MSnExp), "\n",
    "  MS level(s):  ", ifelse(length(object@MSnExp) > 0,
                                paste(sort(unique(msLevel(object@MSnExp))), collapse = ", "),
                                "-"), "\n",
    "  patRoon-class:  ", ifelse(length(object@patdata) > 0,
                                  class(object@patdata), "empty"), "\n",
    sep = ""
  )

  if (nrow(object@annotation$comp) > 0) {
    cat("  Annotation:  ", "Yes", "\n", sep = "")
  }

  if (nrow(object@features) > 0) {
    cat("\n")
    cat("  Number of features:  ", nrow(object@features), "\n", sep = "")
  }

  if (nrow(object@features) > 0) {
    cat("  Filtered features:  ", nrow(object@removed), "\n", sep = "")
  }

  if (length(object@workflows) > 0) {
    cat("\n")
    cat("  Workflow Objects:\n", paste(names(object@workflows), "\n", sep = ""))
  }

})




### samples -----

#' @describeIn ntsData Getter for samples.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
setMethod("samples", "ntsData", function(object) object@samples$sample)




### sampleGroups -----

#' @describeIn ntsData Getter for sample replicate groups.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("sampleGroups", "ntsData", function(object) object@samples$group)

#' @describeIn ntsData Setter for sample replicate groups.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector to assign sample replicate groups.
#'
#' @export
#'
setMethod("sampleGroups<-", signature("ntsData", "ANY"), function(object, value) {

  if (length(value) != length(object@samples$sample)) {
    stop("Length of value does not match the number of samples.")
  }

  object@samples$group <- value
  
  object@metadata <- data.frame(group = unique(value))

  return(object)
})




### blanks -----

#' @describeIn ntsData Getter for blank replicate groups.
#'
#' @param object An \linkS4class{ntsData} object.
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
#' @export
#'
setMethod("blanks<-", signature("ntsData", "ANY"), function(object, value) {

  if (FALSE %in% unique(value %in% object@samples$group)) {
    error("Blank replicate sample groups must be one or more sample replicate groups.")
  }

  object@samples$blank <- value

  return(object)
})




### metadata -----

#' @describeIn ntsData Getter for metadata.
#'
#' @param x An \linkS4class{ntsData} object.
#' @param vars The name of the variables to extract.
#'
#' @export
#'
setMethod("metadata", "ntsData", function(x, vars) {
  
  if (!missing(vars)) {
    if (all(vars %in% colnames(x@metadata))) {
      m <- x@metadata[, c("group", vars)]
    } else {
      warning("vars not find in metadata.")
      m <- x@metadata
    }
  } else {
    m <- x@metadata
  }
  
  return(m)
})




### QC -----

#' @describeIn ntsData Getter for QC replicate sample groups.
#'
#' @param object An \linkS4class{ntsData} object.
#'
#' @export
#'
setMethod("QC", "ntsData", function(object) object@QC@samples$sample)

#' @describeIn ntsData Setter for QC replicate sample groups.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param value A character vector with the names of the sample replicate group/s
#' to be used for QC.
#'
#' @export
#'
setMethod("QC<-", "ntsData", function(object, value) {

  if (FALSE %in% unique(value %in% object@samples$group)) {
    cat("Sample replicate group/s not found in the sample replicate groups.")
  } else {
    object@QC@samples <- object@samples[object@samples$group %in% value, ]
    object@samples <- object@samples[!(object@samples$group %in% value), ]
  }

  return(object)
})




### [ sub-setting samples -----

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
        warning("Given sample names not found in the ntsData object!")
        return(x)
      }
        sn <- i
        sidx <- which(x@samples$sample %in% sn)
    }

    x@samples <- x@samples[x@samples$sample %in% sn,, drop = FALSE]

    x@metadata <- x@metadata[x@metadata$group %in% x@samples$group,, drop = FALSE]

    if (length(x@MSnExp) > 0) x@MSnExp <- MSnbase::filterFile(x@MSnExp, file = sidx)

    if (length(analyses(x@patdata)) > 0) {

      x@patdata <- x@patdata[sidx]

      x@peaks <- x@peaks[x@peaks$sample %in% sn, ]

      if (nrow(x@features) > 0) x <- updateFeatureList(x)

    }

    #annotation, remove replicate groups without samples
    if (nrow(x@annotation$comp) > 0) {
      rg <- unique(x@samples$group)
      x@annotation$comp <- x@annotation$comp[x@annotation$comp$group %in% rg, ]
      x@annotation$raw <- x@annotation$raw[names(x@annotation$raw) %in% rg]
      if (nrow(x@features) > 0) x@annotation$comp <- x@annotation$comp[x@annotation$comp$ID %in% x@features$ID, ]
    }
  }

  return(x)

})




### [ sub-setting features -----

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




### peaks -----

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
