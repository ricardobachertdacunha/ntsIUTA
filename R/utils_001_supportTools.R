

#' @title saveObject
#' @description Saves project objects as RDS in the rdata folder (the default)
#' or as other formats in the results folder, when specified.
#' If not present, the folders are created.
#'
#' @param ... An object to save.
#' @param format The saving format. The default is \code{rds},
#' which is saved in rdata folder.
#' Other formats are saved in the results folder.
#' @param filename A character string with the desired RDS file name.
#' @param path The project path to save the \code{obj}.
#' The default is the working directory as obtained by \code{getwd()}.
#' Not needed if \code{obj} is of class \linkS4class{ntsData}
#' since the path is taken from the slot \code{path}.
#'
#' @export
#'
saveObject <- function(format = "rds",
                       filename = NULL,
                       path = getwd(), ...) {

  dots <- list(...)

  if (class(dots[[1]]) == "ntsData") path <- dots[[1]]@path

  if (format == "rds") {
    rdata <- paste0(path, "\\rdata")
    if (!dir.exists(rdata)) dir.create(rdata)

    if (is.null(filename)) {
      if (class(dots[[1]]) == "ntsData") {
        filename <- "ntsData"
      } else {
        filename <- as.character(names(dots[[1]]))
      }
    }

    saveRDS(dots[[1]], file = paste0(rdata, "\\", filename, ".rds"))

  }

  # TODO Implement saving for other file types, ex. plots and tables

}



#' @title getColors
#'
#' @param x An \linkS4class{ntsData} object with one or more files or the number of colours to be produced.
#' @param which Possible entries are \code{samples}, \code{sampleGroups} and \code{groups} for getting
#' individual sample, sample replicate group and individual group colours, respectively.
#'
#' @return A vector of colours, named according to a \linkS4class{ntsData} object if given.
#'
#' @export
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @importFrom dplyr count
#'
getColors <- function(x, which = c("samples", "sampleGroups", "groups")) {

  colors <- c(RColorBrewer::brewer.pal(8, "Greys")[6],
              RColorBrewer::brewer.pal(8, "Greens")[6],
              RColorBrewer::brewer.pal(8, "Blues")[6],
              RColorBrewer::brewer.pal(8, "Oranges")[6],
              RColorBrewer::brewer.pal(8, "Purples")[6],
              RColorBrewer::brewer.pal(8, "PuRd")[6],
              RColorBrewer::brewer.pal(8, "YlOrRd")[6],
              RColorBrewer::brewer.pal(8, "PuBuGn")[6],
              RColorBrewer::brewer.pal(8, "GnBu")[6],
              RColorBrewer::brewer.pal(8, "BuPu")[6],
              RColorBrewer::brewer.pal(8, "Dark2"))

  if (!class(x) == "numeric" & !class(x) == "integer") {

    if (which != "samples") nameOfcolors <- unique(sampleGroups(x))
    if (which == "samples") nameOfcolors <- samples(x)

    numberOfGroups <- length(nameOfcolors)

    if (numberOfGroups > 18) {
      colors <- grDevices::colorRampPalette(colors)(numberOfGroups)
    }

    vec_colors <- colors[1:numberOfGroups]

    if (which == "groups") {
      names(vec_colors) <- nameOfcolors
      return(vec_colors)
    }

    if (which == "sampleGroups") {
      count <- count(x@samples[, c("sample", "group")], group)
      vec_colors <- rep(vec_colors, times = count[, "n"])
      names(vec_colors) <- samples(x)
      return(vec_colors)
    }

    if (which == "samples") {
      names(vec_colors) <- samples(x)
      return(vec_colors)
    }

  } else {

    numberOfGroups <- x

    if (numberOfGroups > 18) {
      colors <- grDevices::colorRampPalette(colors)(numberOfGroups)
    }

    vec_colors <- colors[1:numberOfGroups]

    return(vec_colors)
  }

}



#' @title mzrBuilder
#'
#' @param mz The target \emph{m/z} to calculate range.
#' @param ppm The mass deviation in ppm to calculate the \emph{m/z} range.
#'
mzrBuilder <- function(mz = NULL, ppm = NULL) {

  mzr <- NULL

  if (length(mz) == 1 && !is.null(mz)) {
    if (is.null(ppm)) ppm <- 20
    mzr <- c(mz - ((ppm / 1E6) * mz), mz + ((ppm / 1E6) * mz))
  }

  if (length(mz) == 2) mzr <- c(mz[1], mz[2])

  return(mzr)

}



#' @title rtrBuilder
#'
#' @param rt The target retention time to calculate range.
#' @param rtWindow The retention time window to calculate the range.
#' Can be length 1 (default to 1 minute), for a fixed deviation of the given \code{rt} or  length 2.
#' When length 2, it overwrites a given \code{rt} and uses the two values to calculate the range.
#' @param rtUnit The time unit for the time arguments. Possible values are "sec" and "min".
#'
rtrBuilder <- function(rt = NULL, rtWindow = NULL, rtUnit = "sec") {

  rtr <- NULL

  if (rtUnit == "min") if (!is.null(rt)) rt <- rt * 60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow * 60

  if (!is.null(rt)) {
    rtr <- c((rt) - ifelse(!is.null(rtWindow), rtWindow, 60),
            (rt) + ifelse(!is.null(rtWindow), rtWindow, 60))
  }

  if (unique(!is.null(rtWindow))) if (length(rtWindow) == 2) {
    rtr <- c(rtWindow[1], rtWindow[2])
  }

  return(rtr)

}



#' @title filterFileFaster
#' @description Filter files (i.e., samples) but it does not redo feature list.
#' Useful for plotting data or extracting EICs from certain samples.
#'
#' @param x An \linkS4class{ntsData} object.
#' @param i The indices or names of the samples to keep.
#'
#' @return The sub-setted \linkS4class{ntsData} object.
#'
#' @importMethodsFrom MSnbase filterFile
#'
filterFileFaster <- function(x, i) {

  if (is.null(i)) return(x)

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

  rgr <- unique(x@samples$group[!(x@samples$group %in% unique(x@samples$group[sidx]))])

  x@samples <- x@samples[x@samples$sample %in% sn,, drop = FALSE]

  x@metadata <- x@metadata[x@metadata$sample %in% sn,, drop = FALSE]

  x@MSnExp <- filterFile(x@MSnExp, file = sidx)
  
  if (length(analyses(x@patdata)) > 0) {
    x@patdata <- x@patdata[sidx]
  }
  
  if (nrow(x@peaks) > 0) x@peaks <- x@peaks[x@peaks$sample %in% sn, ]

  if (nrow(x@features) > 0) {
    x@features <- x@features[x@features$ID %in% names(x@patdata), ]
  }

  return(x)

}



#' filterFeatureGroups
#'
#' @param x A \linkS4class{featureGroups} object.
#' @param i The indices or the sample names to subset \code{x}.
#'
#' @return A subset of the \linkS4class{featureGroups} object.
#'
#' @importClassesFrom patRoon featureGroups featureGroupsScreening
#' @importClassesFrom xcms XCMSnExp
#' @importMethodsFrom patRoon analyses groupTable
#' @importMethodsFrom MSnbase fileNames
#' @importMethodsFrom xcms filterFile
#'
filterFeatureGroups <- function(x, i) {

  if (!is.character(i)) {
    sn <- patRoon::analyses(x)[i]
    sidx <- which(patRoon::analyses(x) %in% sn)
  } else {
    sn <- i
    sidx <- which(patRoon::analyses(x) %in% sn)
  }

  if (class(x) == "featureGroupsXCMS3" | class(x@features) == "featuresXCMS3") {

    x@groups <- x@groups[sidx]

    x@ftindex <- x@ftindex[sidx]

    #remove featuresGroups with 0 intensity
    zeros <- patRoon::groupTable(x)
    zeros <- apply(zeros, MARGIN = 2, FUN = function(z) sum(z))
    zeros <- !(zeros == 0)

    x@groups <- x@groups[, ..zeros]
    x@ftindex <- x@ftindex[, ..zeros]
    x@groupInfo <- x@groupInfo[zeros, ]

    x@analysisInfo <- x@analysisInfo[sidx, ]
    x@features@features <- x@features@features[sidx]
    x@features@analysisInfo <- x@features@analysisInfo[x@features@analysisInfo$analysis %in% sn, ]

    if (x@algorithm == "xcms3") {
      files <- basename(fileNames(x@xdata)[sidx])
      x@xdata <- filterFile(x@xdata, files, keepFeatures = TRUE)
      x@features@xdata <- filterFile(x@features@xdata, files, keepFeatures = TRUE)
    }

    if (class(x) == "featureGroupsScreening") {
      #x@screenInfo <- x@screenInfo[zeros, ] #x@screenInfo[x@screenInfo$groups %in% x@groups, ]
      x@screenInfo[x@screenInfo$groups %in% x@groups, ]
    }

  } else {

    x <- x[sidx, ]

  }

  return(x)

}



#' extractMS2
#'
#' @param obj A \linkS4class{featureGroups} or an \linkS4class{ntsData} object to
#' extract MS2 data according to percursor ions.
#' @param param A \linkS4class{MS2param} object with parameters for MS2 extraction.
#'
#' @return A \linkS4class{MSPeakLists} when \code{obj} is a \linkS4class{featureGroups}
#' or an \linkS4class{ntsData} updated with MS2 data
#' when \code{obj} is an \linkS4class{ntsData} object.
#'
#' @export
#'
#' @importFrom patRoon getDefAvgPListParams
#' @importMethodsFrom patRoon generateMSPeakLists
#'
#'
extractMS2 <- function(obj = NULL,
                       param = MS2param()) {

  pat <- obj

  if (class(obj) == "ntsData") pat <- pat@patdata

  control_avgPListParams <- getDefAvgPListParams(
    clusterMzWindow = param@clusterMzWindow,
    topMost = param@topMost,
    minIntensityPre = param@minIntensityPre,
    minIntensityPost = param@minIntensityPost
  )

  MS2 <- suppressWarnings(generateMSPeakLists(
    pat, "mzr",
    maxMSRtWindow = param@maxMSRtWindow,
    precursorMzWindow = param@precursorMzWindow,
    avgFeatParams = control_avgPListParams,
    avgFGroupParams = control_avgPListParams
  ))

  # TODO add option to add MS2 to the ntsData object

  return(MS2)

}
