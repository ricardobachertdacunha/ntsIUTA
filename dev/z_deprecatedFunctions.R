

#' @title getInstParam
#' @description Tranforms the list of parameters for peak picking (PP), alignment (\emph{i.e.} retention time adjustment across samples)
#' and grouping into the respetive class objects.
#' As input a list as follows should be given where for each data processing step the input should be given as described in
#' the documentation of the aimed function of the \pkg{xcms} package.
#' \itemize{
#'   \item \code{paramList <- list(PP = list("..."), preGrouping = list("..."), alignment = list("..."), grouping = list("..."))}
#' }
#' \code{preGrouping} is used for \code{\link[xcms]{adjustRtime}} function
#' when using the method \code{PeakGroups} as it requires peaks groups (\emph{i.e.} features) for alignment.
#' Note that currently for PP only \code{MassifquantParam} is possible, for alignment only \code{PeakGroupsParam} is possible
#' and for grouping only \code{PeakDensityParam}. Other methods will be integrated in the future.
#' 
#' @param paramList A list of parameters for PP, alignment and grouping as described above.
#' See \code{\link{paramListExample}} object for a structure example.
#'
#' @return A \code{list} with parameters for PP, aligment and grouping to be passed to \code{\link{peakPicking}} and \code{\link{makeFeatures}}. 
#' 
#'  @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#'
#' @examples
#' paramList <- ntsIUTA::paramListExample
#' getInstParam(paramList)
#' 
getInstParam <- function(paramList = paramList) {
  
  instParam <- list()
  
  instParam$PP <- NULL
  
  # Make instParam for PP
  if(!is.null(paramList$PP)) {
    
    if(paramList$PP$funcName == class(xcms::MassifquantParam())) {
      
      instParam$PP <- xcms::MassifquantParam(ppm = paramList$PP$ppm,
                                             peakwidth = paramList$PP$peakwidth,
                                             snthresh = paramList$PP$snthresh,
                                             prefilter = paramList$PP$prefilter,
                                             mzCenterFun = paramList$PP$mzCenterFun,
                                             integrate = paramList$PP$integrate,
                                             mzdiff = paramList$PP$mzdiff,
                                             fitgauss = paramList$PP$fitgauss,
                                             noise = paramList$PP$noise,
                                             criticalValue = paramList$PP$criticalValue,
                                             consecMissedLimit = paramList$PP$consecMissedLimit,
                                             unions = paramList$PP$unions,
                                             checkBack = paramList$PP$checkBack,
                                             withWave = paramList$PP$verboseColumns,
                                             verboseColumns = paramList$PP$verboseColumns)
      
    }
    
    #TODO Add parameters for other methods
    
  }
  
  instParam$preGrouping <- NULL
  instParam$alignment <- NULL
  
  # Make instParam for alignment, including preGrouping for PeakGroupsParam method
  if(!is.null(paramList$alignment)) {
    
    if(paramList$alignment$funcName == class(xcms::PeakGroupsParam())) {
      
      instParam$preGrouping <-  xcms::PeakDensityParam(sampleGroups = "placeHolder",
                                                       bw = paramList$preGrouping$bw,
                                                       minFraction = paramList$preGrouping$minFraction,
                                                       minSamples = paramList$preGrouping$minSamples,
                                                       binSize = paramList$preGrouping$binSize,
                                                       maxFeatures = paramList$preGrouping$maxFeatures)
      
      instParam$alignment <- xcms::PeakGroupsParam(minFraction = paramList$alignment$minFraction,
                                                   extraPeaks = paramList$alignment$extraPeaks,
                                                   smooth = paramList$alignment$smooth,
                                                   span = paramList$alignment$span,
                                                   family = paramList$alignment$family,
                                                   peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
                                                   subset = integer(),
                                                   subsetAdjust = "average")
      
    }
    
    #Add parameters for other methods
    
  }
  
  instParam$grouping <- NULL
  
  #Make instParam for grouping.
  if(!is.null(paramList$grouping)) {
    
    if(paramList$grouping$funcName == class(xcms::PeakDensityParam(sampleGroups = "placeHolder"))) {
      
      instParam$grouping <- PeakDensityParam(sampleGroups = "placeHolder",
                                             bw = paramList$grouping$bw,
                                             minFraction = paramList$grouping$minFraction,
                                             minSamples = paramList$grouping$minSamples,
                                             binSize = paramList$grouping$binSize,
                                             maxFeatures = paramList$grouping$maxFeatures)
      
    }
    
    #Add parameters for other methods
    
  }
  
  return(instParam)
  
}



### paramList -----

# paramListExample <- list(
#   
#   instName = "Generic",
#   
#   PP = list(
#     funcName = "MassifquantParam",
#     ppm = 20,
#     peakwidth = c(8, 60),
#     snthresh = 3,
#     prefilter = c(6, 500),
#     mzCenterFun = "wMean",
#     integrate = 1,
#     mzdiff = -0.001,
#     fitgauss = TRUE,
#     noise = 200,
#     criticalValue = 1.5,
#     consecMissedLimit = 2,
#     unions = 1,
#     checkBack = 1,
#     withWave = TRUE,
#     verboseColumns = TRUE
#   ),
#   
#   preGrouping = list(
#     funcName = "PeakDensityParam",
#     bw = 5,
#     minFraction = 0.5,
#     minSamples = 1,
#     binSize = 0.003,
#     maxFeatures = 100
#   ),
#   alignment = list(
#     funcName = "PeakGroupsParam",
#     minFraction = 1,
#     extraPeaks = 0,
#     smooth = "loess",
#     span = 0.2,
#     family = "gaussian"
#   ),
#   grouping = list(
#     funcName = "PeakDensityParam",
#     bw = 3,
#     minFraction = 0.5,
#     minSamples = 1,
#     binSize = 0.003,
#     maxFeatures = 100
#   ),
#   componentization = list(
#     
#   )
# )


#additional function for ploting adepted XIC from MSnbase 
.vertical_sub_layout_ex <- function(x, sub_plot = 2) {
  sqrt_x <- base::sqrt(base::length(x))
  ncol <- base::ceiling(sqrt_x)
  nrow <- base::round(sqrt_x)
  rws <- base::split(1:(ncol * nrow * sub_plot), f = base::rep(1:nrow, each = sub_plot * ncol))
  base::do.call(rbind, base::lapply(rws, matrix, ncol = ncol))
}



#additional function for ploting adepted XIC from MSnbase
.plotXIC_ex <- function(x, main = "", col = "white", colramp = colorRampPalette(c("#383E47", "#5E8CAA", "#16B9E5", "#16E5C9", "#16E54C")), #topo.colors colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdGy")))
                        grid.color = "lightgrey", pch = 21,
                        layout = base::matrix(1:2, ncol = 1), plotTargetMark = plotTargetMark,...) {
  # start edit
  plot_strip <- function(..., v, h) base::plot(...)
  dots <- base::list(...)
  #print(list(...))
  # end edit
  
  if (!plotTargetMark)
  {
    dots$h <- NULL
    dots$v <- NULL
  }
  
  targetColor <- "#0644E9"
  
  if (base::is.matrix(layout))
    graphics::layout(layout)
  
  ## Chromatogram.
  bpi <- base::unlist(base::lapply(base::split(x$i, x$rt), max, na.rm = TRUE))
  brks <- lattice::do.breaks(base::range(x$i), nint = 256)
  graphics::par(mar = c(0, 4, 2, 1))
  
  # start edit
  plot_strip(base::as.numeric(base::names(bpi)), bpi, xaxt = "n", col = col, main = main,
             bg = lattice::level.colors(bpi, at = brks, col.regions = colramp), xlab = "",
             pch = pch, cex = 1.5, ylab = "", las = 2, ...)
  # end edit
  
  if(!is.null(dots$v)) graphics::abline(v=dots$v,col=targetColor,lty=3)
  
  graphics::mtext(side = 4, line = 0, "Intensity", cex = graphics::par("cex.lab"))
  graphics::grid(col = grid.color)
  graphics::par(mar = c(3.5, 4, 0, 1))
  
  # start edit
  plot_strip(x$rt, x$mz, main = "", pch = pch, col = col, xlab = "", ylab = "", cex = 1.5,
             yaxt = "n", bg = lattice::level.colors(x$i, at = brks, col.regions = colramp),
             ...)
  
  if(!is.null(dots$h)) graphics::abline(h=dots$h,col=targetColor,lty=3)
  if(!is.null(dots$v)) graphics::abline(v=dots$v,col=targetColor,lty=3)
  
  if(!is.null(dots$v) & !is.null(dots$h))
  {
    graphics::rect(dots$v-10, dots$h-((5/1E6)*dots$h), dots$v+10, dots$h+((5/1E6)*dots$h),
                   col = NA, lty = 2, border = targetColor)
  }
  # end edit
  
  graphics::axis(side = 2, las = 2)
  graphics::grid(col = grid.color)
  graphics::mtext(side = 1, line = 2.5, "Retention time (sec.)", cex = par("cex.lab"))
  graphics::mtext(side = 4, line = 0, "m/z", cex = par("cex.lab"))
  
}


#' @title plotRawChrom_old
#' @description Plot BPC or TIC chromatograms of an \linkS4class{OnDiskMSnExp} object.
#'
#' @param x An \linkS4class{OnDiskMSnExp} object with one or more files.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param mz Optional target \emph{m/z} to obtain an extracted ion chromatogram (EIC).
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{sec}.
#' @param msLevel The MS level to extract the data. POssible values are 1 or 2.
#' @param type The type of chromatogram. Possible entries are "bpc" of base peak chromatogram or "tic" for total ion chromatogram.
#' The default is "tic".
#'
#' @return A plot.
#'
#' @examples
#' 
#' plotRawChrom_old(ntsIUTA::rawDataExample)
#' 
#' 
plotRawChrom_old <- function(x = rawData, fileIndex = NULL,
                             mz = NULL, ppm = 20,
                             rt = NULL, rtWindow = NULL,
                             rtUnit = "min",
                             msLevel = 1, type = "tic") {
  
  require(magrittr)
  
  if (rtUnit == "min") if (!is.null(rt)) rt <- rt*60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow*60
  
  if (!base::is.null(fileIndex)) {x <- MSnbase::filterFile(x, fileIndex)}
  
  rtr <- c(base::min(MSnbase::rtime(x)), base::max(MSnbase::rtime(x)))
  mzr <- NULL
  if (!is.null(mz)) {
    
    if (length(mz) == 1) { mzr <- c(mz - ((ppm/1E6)*mz), mz + ((ppm/1E6)*mz)) }
    if (length(mz) == 2) { mzr <- c(mz[1], mz[2]) }
    
  }
  
  if (!is.null(rt)) { rtr <- c((rt) - ifelse(!is.null(rtWindow), rtWindow, 1), (rt) + ifelse(!is.null(rtWindow), rtWindow, 1)) }
  if (is.null(rt)) if (unique(!is.null(rtWindow))) if (length(rtWindow) == 2) { rtr <- c(rtWindow[1], rtWindow[2]) }
  
  colors <- getColors(x, "samples")
  
  if (is.null(mzr)) {
    chrom <- MSnbase::chromatogram(x, aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))), rt = rtr, msLevel = msLevel)
  } else {
    chrom <- MSnbase::chromatogram(x, aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))), rt = rtr, mz = mzr, msLevel = msLevel)
  }
  
  if (!base::is.null(mzr))
  {
    main <- base::paste0("EIC of ", round(mz, digits = 4)," +/- ", round(ppm, digits = 0)," ppm")
  } else {
    main <- base::toupper(type)
  }
  
  # graphics::layout(base::matrix(1:1))
  # graphics::par(mar=c(5,4,4,2)+0.1)
  
  plot <- MSnbase::plot(chrom, col = colors, lwd = 0.5, main = main, xlab = "Retention time (sec.)", ylab = "Intensity",)
  
  plot <- plot %>% graphics::legend("topleft" , legend = x$sample_name, col = colors, lty = 1, lwd = 2, cex = 0.6, bty = "n", xjust = 0)
  
  
  return(plot)
  
}

