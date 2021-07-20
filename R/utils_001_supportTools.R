

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



#' @title getScreeningListTemplate
#'
#' @param projPath The project folder location. Default is \code{setup$projPath}.
#'
#' @return Pastes a template .csv file of the screeningList into the projPath.
#'
#' @export
#'
getScreeningListTemplate <- function(projPath = setup$projPath) {
  base::file.copy(from = base::paste0(base::system.file(package = "ntsIUTA", dir = "extdata"), "/screeningList_template.csv"),
                  to = setup$projPath,
                  overwrite = FALSE)
}



#' @title getColors
#'
#' @param x An \linkS4class{ntsData} object with one or more files or the number of colours to be produced.
#' @param which Possible entries are \code{samples}, \code{samplegroups} and \code{groups} for getting
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
getColors <- function(x, which = c("samples", "samplegroups", "groups")) {

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

    if (which != "samples") nameOfcolors <- unique(pullSamplegroups(x))
    if (which == "samples") nameOfcolors <- pullSamples(x)

    numberOfGroups <- length(nameOfcolors)
    
    if (numberOfGroups > 18) {
      colors <- grDevices::colorRampPalette(colors)(numberOfGroups)
    }

    vec_colors <- colors[1:numberOfGroups]

    if (which == "groups") {
      names(vec_colors) <- nameOfcolors
      return(vec_colors)
    }

    if (which == "samplegroups") {
      count <- count(x@samples[,c("sample","group")], group)
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
#' @param x An \linkS4class{OnDiskMSnExp} object.
#' @param i The indices or names of the samples to keep.
#'
#' @return The sub-setted \linkS4class{OnDiskMSnExp} object.
#'
#' @importMethodsFrom MSnbase filterFile
#'
filterFileFaster <- function(x, i) {
  
  if (!is.character(i)) {
    sn <- x@samples$sample[i]
    sidx <- which(x@samples$sample %in% sn)
  } else {
    sn <- i
    sidx <- which(x@samples$sample %in% sn)
  }
  
  rgr <- unique(x@samples$group[!(x@samples$group %in% unique(x@samples$group[sidx]))])
  
  x@samples <- x@samples[x@samples$sample %in% sn, , drop = FALSE]
  
  x@metadata <- x@metadata[x@metadata$sample %in% sn, , drop = FALSE]
  
  x@MSnExp <- filterFile(x@MSnExp, file = sidx)
  
  if (nrow(x@peaks) > 0) x@peaks <- x@peaks[x@peaks$sample %in% sn, , drop = FALSE]
  
  if (nrow(x@features) > 0) {
    rg <- unique(x@samples$group)
    x@features <- x@features[, !(colnames(x@features) %in% rgr)]
    x@features <- x@features[!sapply(x@features[, rg] == 0, function(x) sum(x) == length(rg)), ]
  }
  
  return(x)

}
