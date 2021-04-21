

#' @title plotRawChrom
#' @description Plot BPC or TIC chromatograms of an \linkS4class{OnDiskMSnExp} object using \pkg{plotly} to procude an iteractive plot.
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
#' @return An iterative 
#' 
#' @export
#'
#' @import magrittr
#' @importFrom MSnbase filterFile rtime chromatogram
#' @importFrom BiocGenerics as.data.frame
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
#' @examples
#' 
#' 
#' 
plotRawChrom <- function(x = rawData, fileIndex = NULL,
                         mz = NULL, ppm = 20,
                         rt = NULL, rtWindow = NULL,
                         rtUnit = "min",
                         msLevel = 1, type = "tic") {
  
  
  # x = rawDataExample
  # fileIndex = 1:2
  # mz <- 233.0243
  # rt <- 15.6809
  # rtUnit = "min"
  # ppm <- 20
  # rtWindow = 3
  # msLevel = 1
  
  #require(magrittr)
  
  if (rtUnit == "min") if (!is.null(rt)) rt <- rt*60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow*60
  
  if (!base::is.null(fileIndex)) {x <- MSnbase::filterFile(x, fileIndex)}
  
  rtr <- c(base::min(MSnbase::rtime(x)), base::max(MSnbase::rtime(x)))
  mzr <- NULL
  if (!base::is.null(mz)) {
    
    if (base::length(mz) == 1) { mzr <- c(mz - ((ppm/1E6)*mz), mz + ((ppm/1E6)*mz)) }
    if (base::length(mz) == 2) { mzr <- c(mz[1], mz[2]) }
    
  }
  
  if (!base::is.null(rt)) { rtr <- c((rt) - base::ifelse(!base::is.null(rtWindow), rtWindow, 1), (rt) + base::ifelse(!base::is.null(rtWindow), rtWindow, 1)) }
  if (base::is.null(rt)) if (base::unique(!base::is.null(rtWindow))) if (base::length(rtWindow) == 2) { rtr <- c(rtWindow[1], rtWindow[2]) }
  
  if (!base::is.null(mzr))
  {
    main <- base::paste0("EIC of ", base::round(mz, digits = 4)," +/- ", base::round(ppm, digits = 0)," ppm")
  } else {
    main <- base::toupper(type)
  }
  
  colors <- ntsIUTA::getColors(x, "samples")
  
  if (is.null(mzr)) {
    chrom <- MSnbase::chromatogram(x,
                                   aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))),
                                   rt = rtr, msLevel = msLevel)
  } else {
    chrom <- MSnbase::chromatogram(x,
                                   aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))),
                                   rt = rtr, mz = mzr, msLevel = msLevel)
  }
  
  df <- chrom[[1]]
  df <- BiocGenerics::as.data.frame(df)
  df$file <- x$sample_name[1]
  if (base::length(chrom)>1)
  {
    for (i in 2:base::length(chrom)) {
      temp <- chrom[[i]]
      temp <- BiocGenerics::as.data.frame(temp)
      temp$file <- x$sample_name[i]
      df <- base::rbind(df,temp)
    }
  }
  
  title <- base::list(text = main, x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))
  
  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                      titlefont = base::list(size = 12, color = "black"))
  
  yaxis = base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                     titlefont = base::list(size = 12, color = "black"))
  
  
  plot <- plotly::plot_ly(df,
                          x = df[df$file == x$sample_name[1],"rtime"],
                          y = df[df$file == x$sample_name[1],"intensity"],
                          type = "scatter", mode = "lines+markers",
                          line = base::list(width = 0.5, color = base::unname(ntsIUTA::getColors(x,"samples")[1])),
                          marker = base::list(size = 2, color = base::unname(ntsIUTA::getColors(x,"samples")[1])),
                          name = x$sample_name[1])
  
  if (base::length(chrom) > 1)
  {
    for (i in 2:base::length(chrom)) {
      plot  <- plot %>% plotly::add_trace(df,
                                          x = df[df$file == x$sample_name[i],"rtime"],
                                          y = df[df$file == x$sample_name[i],"intensity"],
                                          type = "scatter", mode = "lines+markers",
                                          line = base::list(width = 0.5, color = base::unname(ntsIUTA::getColors(x,"samples")[i])),
                                          marker = base::list(size = 2, color = base::unname(ntsIUTA::getColors(x,"samples")[i])),
                                          name = x$sample_name[i])
    }
  }
  
  plot <- plot %>% plotly::layout(legend = base::list(title = base::list(text='<b> Sample: </b>')),
                                  xaxis = xaxis,yaxis = yaxis, title = title)
  
  return(plot)
  
}




