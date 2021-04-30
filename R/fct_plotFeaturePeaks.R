


#' @title plotFeaturePeaks
#' @description Plots feature information for a given \linkS4class{XCMSnExp} object containing grouped peaks across samples.
#'
#' @param x An \linkS4class{XCMSnExp} or \linkS4class{featureGroups} object with one or more files and grouped peaks.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param features The identifier of the features of interest.
#' @param mz Optional target \emph{m/z} to find features using the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the features when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below. Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect features. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{min}.
#' @param plotBy Grouping for legend. Possible entries are: \code{samples}, \code{features} or \code{groups}. 
#' @param names A character string with names for each feature given in \code{features}.
#' Note that length should match between \code{names} and \code{features}.
#'
#' @return A double plot with peaks chromatograms on the top part and feature peak groups below.
#' 
#' @export
#' 
#' @import magrittr
#' @importFrom xcms filterFile hasFeatures featureDefinitions featureChromatograms chromPeaks
#' @importFrom MSnbase rtime
#' @importFrom BiocGenerics as.data.frame
#' @importFrom patRoon as.data.table
#' @importFrom plotly plot_ly add_trace layout hide_colorbar subplot toRGB
#' @importFrom stats setNames
#'
#' @examples
#' 
#' 
#' 
plotFeaturePeaks <- function(x = featData, fileIndex = NULL, 
                             features = NULL,
                             mz = NULL, ppm = 5,
                             rt = NULL, rtWindow = 1,
                             rtUnit = "min", plotBy = "samples",
                             names = NULL) {
  
  # x = qcPat #featData
  # fileIndex = NULL
  # features = qcdf$group #c("FT0084", "FT0029", "FT1113")
  # mz = NULL
  # ppm = NULL
  # rt = NULL
  # rtWindow = NULL
  # rtUnit = "min"
  # plotBy = "features"
  # names = qcdf$name
  
  #library(magrittr)
  
  #check if is x is XCMSnExp (xcms) or featureGroups (patRoon)
  if (base::class(x) == "featureGroupsXCMS3") {
    if (!base::is.null(fileIndex)) {x <- x[fileIndex,]}
    y <- x@xdata
  } else {
    if (base::class(x) == "XCMSnExp") {
      if (!base::is.null(fileIndex)) {x <- xcms::filterFile(x, fileIndex, keepFeatures = TRUE,  keepAdjustedRtime = TRUE)}
      y <- x
    } else {
      stop("x must be an XCMSnExp (xcms) or featureGroups (patRoon) object.")
    }
  }
  
  if (!xcms::hasFeatures(y)) stop ("Features not found in the given XCMSnExp object.")
  
  rtr <- c(base::min(MSnbase::rtime(y)), base::max(MSnbase::rtime(y)))
  if (!base::is.null(rt)) if (rtUnit == "min") rt <- rt*60
  if (!base::is.null(rtWindow)) if (rtUnit == "min") rtWindow <- rtWindow*60
  if (!base::is.null(rt) & !base::is.null(rtWindow)) { rtr <- c((rt) - rtWindow, (rt) + rtWindow) }
  if (base::is.null(rt)) if (!base::is.null(rtWindow)) if (base::length(rtWindow) == 2) { rtr <- rtWindow }
  
  #When priotity is for feature and than mz
  if (!base::is.null(features)) {
    if (base::unique(base::grepl("M*_R", features, fixed = FALSE))) {
      gKey <- base::cbind(patRoon::as.data.table(x, average = TRUE)[,.SD, .SDcols = "group"],
                          base::data.frame(FT = base::row.names(xcms::featureDefinitions(y))))
      gKey <- base::as.data.frame(gKey)
      FT <- gKey[gKey$group %in% features, "FT", drop = T]
    } else {
      FT <- features
    }
    #When features are no given but a specific mz +/- ppm  
  } else {
    if (base::is.null(mz)) stop("If features are not given, at least mz should be given.")
    FT <- base::row.names(xcms::featureDefinitions(y, mz = mz, ppm = ppm, rt = rtr, type = "within", msLevel = 1))
  }
  
  if (base::length(FT) == 0 | base::is.null(FT)) stop ("Features were not found with the given paramters.") 
  
  defFT <- xcms::featureDefinitions(y)
  defFT <- defFT[base::row.names(defFT) %in% FT,]
  
  #Make features legend
  features <- base::character()
  for (f in 1:base::length(FT)) {
    features <- c(features, base::paste("M", base::round(defFT$mzmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                        "_R", base::round(defFT$rtmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                        "_", FT[f],sep=""))
  }
  
  featuresID <- features
  
  if (!base::is.null(names) & base::length(names) == base::length(features)) {
    features <- names
  }
  
  chrom <- xcms::featureChromatograms(y, features = FT, aggregationFun = "sum",
                                      expandRt = rtr[2]-rtr[1], include = "any",
                                      filled = TRUE, missing = NA)
  
  df <- base::data.frame(rtime = base::as.numeric(),intensity = base::as.numeric(),sample = base::as.character(),FT = base::as.character())
  for (f in 1:base::nrow(defFT)) {
    for (s in 1:base::ncol(chrom)) {
      temp <- chrom[f,s]
      temp <- BiocGenerics::as.data.frame(temp)
      temp$sample <- y$sample_name[s]
      temp$FT <- FT[f]
      df <- base::rbind(df,temp)
    }
  }
  
  if (plotBy == "samples") {
    colors <- ntsIUTA::getColors(y, "samples")
  } else {
    if (plotBy == "features") {
      colors <- ntsIUTA::getColors(base::length(features))
    } else {
      colors = ntsIUTA::getColors(y, "groups")
    }
  }
  
  p1 <- plotly::plot_ly(df)
  showlegend = base::rep(0,base::length(features))
  for (s in 1:base::length(y$sample_name)) { 
    for (f in 1:base::length(FT)) {
      rtFT <- base::as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
      if(base::unique(TRUE %in% (rtFT$sample == s))) {
        showlegend[f] = 1 + showlegend[f]
        p1 <- p1 %>% plotly::add_trace(df,
                                       x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "rtime"],
                                       y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "intensity"],
                                       type = "scatter", mode = "lines",
                                       line = list(width = 0.5,
                                                   color = ifelse(plotBy == "samples",base::unname(colors[s]),
                                                                  ifelse(plotBy == "features", colors[f],base::unname(colors[s])))),
                                       connectgaps = TRUE,
                                       name = ifelse(plotBy == "samples",y$sample_name[s],
                                                     ifelse(plotBy == "features", features[f], y$sample_group[s])),
                                       legendgroup = ifelse(plotBy == "samples",y$sample_name[s],
                                                            ifelse(plotBy == "features", features[f], y$sample_group[s])),
                                       showlegend = ifelse(plotBy == "samples",base::ifelse(f == 1, T, F),
                                                           ifelse(plotBy == "features", base::ifelse(showlegend[f] == 1, T, F),
                                                                  base::ifelse(f == 1, T, F))))
        
        
        
        rtFT <- rtFT[rtFT$sample == s, c("rtmin","rtmax")]
        p1 <- p1 %>%  plotly::add_trace(df,
                                        x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2], "rtime"],
                                        y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2],"intensity"],
                                        type = "scatter", mode = "lines", fill = 'tozeroy', connectgaps = TRUE,
                                        fillcolor = ifelse(plotBy == "samples",paste(color = base::unname(colors[s]),50, sep = ""),
                                                           ifelse(plotBy == "features", paste(color = base::unname(colors[f]),50, sep = ""),
                                                                  paste(color = base::unname(colors[s]),50, sep = ""))),
                                        line = list(width = 0.1,
                                                    color = ifelse(plotBy == "samples",base::unname(colors[s]),
                                                                   ifelse(plotBy == "features", colors[f],base::unname(colors[s])))),
                                        name = ifelse(plotBy == "samples",y$sample_name[s],
                                                      ifelse(plotBy == "features", features[f],y$sample_group[s])),
                                        legendgroup = ifelse(plotBy == "samples",y$sample_name[s],
                                                             ifelse(plotBy == "features", features[f],y$sample_group[s])),
                                        showlegend = F) #base::ifelse(showlegend == 1, T, F)
      }
    }
  }
  
  p2 <- plotly::plot_ly()
  rect <- list()
  colorsRect <- ntsIUTA::getColors(x = base::length(features))
  dotsColor <- c("#000000","#FDFEFE00")
  dotsColor <- stats::setNames(dotsColor, c("0","1"))
  for (f in 1:base::length(FT)) {
    rtFT <- base::as.data.frame(xcms::chromPeaks(y, isFilledColumn = TRUE)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
    
    rect[[f]] <- list(type = "rect", fillcolor = base::paste(colorsRect[f],"15", sep = ""),
                      line = list(color = colorsRect[f], width = 1, dash = 'dash'), #'rgba(62, 186, 32,0.15)'
                      x0 = base::min(rtFT$rtmin, na.rm = T),
                      x1 = base::max(rtFT$rtmax, na.rm = T),
                      xref = "x", y0 = 1, y1 = base::max(base::length(y$sample_name)), yref = "y")
    
    p2 <- p2 %>% plotly::add_trace(rtFT,
                                   x = rtFT$rt,
                                   y = rtFT$sample, type = "scatter", mode = "markers",
                                   color = as.character(rtFT$is_filled), colors = dotsColor, size = 6,
                                   marker = list(line = list(color = colorsRect[f], width = 3)),
                                   showlegend = FALSE,
                                   hoverinfo = 'text', text = paste('</br> feature: ', featuresID[f],
                                                                    '</br> sample: ', y$sample_name[rtFT$sample],
                                                                    '</br> height: ', base::round(rtFT$maxo, digits = 0),
                                                                    '</br> width: ', base::round(rtFT$rtmax-rtFT$rtmin, digits = 0),
                                                                    '</br> dppm: ', rtFT$dppm,
                                                                    '</br> sn: ', rtFT$sn,
                                                                    '</br> egauss: ', base::round(rtFT$egauss, digits = 3),
                                                                    '</br> filled: ', base::ifelse(rtFT$is_filled == 1, "TRUE", "FALSE")))
  }
  p2 <- p2 %>% plotly::layout(shapes = rect)
  p2 <- plotly::hide_colorbar(p2)
  
  plotList <- list()
  plotList[[paste0("p1",s)]]<- p1
  plotList[[paste0("p2",s)]]<- p2
  
  title <- base::list(text = "Coisas", x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))
  
  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                      titlefont = base::list(size = 12, color = "black"),
                      range = c(rtr), autotick = T, ticks = "outside")
  
  yaxis1 = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                titlefont = list(size = 12, color = "black"))
  
  yaxis2 = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Sample",
                titlefont = list(size = 12, color = "black"), tick0 = 0, dtick = 1)
  
  plot <- plotly::subplot(plotList, nrows = 2, margin = 0.04, shareX = TRUE, which_layout = "merge")
  plot <- plot %>% plotly::layout(xaxis = xaxis, yaxis = yaxis1, yaxis2 = yaxis2)
  
  return(plot)
  
}




