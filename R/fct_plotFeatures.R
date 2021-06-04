

#' @title plotFeatures
#' @description Iterative plot for features from an \linkS4class{XCMSnExp} or \linkS4class{featureGroups} containing grouped peaks.
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
#' @param title The title for the plot, optional.
#' @param plotBy The grouping for the plot. Possible are "samples", "replicates" and "features"
#' for ploting by individual samples, relpicate samples or features, respectively. The default is "samples".
#'
#' @return An iterative plot.
#' 
#' @export
#' 
#' @import magrittr
#' @importFrom xcms filterFile featureDefinitions featureChromatograms chromPeaks
#' @importMethodsFrom MSnbase rtime
#' @importFrom BiocGenerics as.data.frame
#' @importFrom patRoon as.data.table
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
#' @examples
#' 
#' 
#' 
plotFeatures <- function(x = patData, fileIndex = NULL, 
                         features = NULL,
                         mz = NULL, ppm = 5,
                         rt = NULL, rtWindow = 1,
                         rtUnit = "min",
                         title = NULL, plotBy = c("samples")) {
  
  #Examples
  # plotFeatures(ntsIUTA::patDataExample, features = "M233_R941_36")
  
  # y = featData
  # x = patData
  # fileIndex = NULL
  # features = "M233_R941_107"
  # mz <- NULL
  # rt <- NULL
  # rtUnit = "min"
  # ppm <- NULL
  # rtWindow = 2
  # plotBy = "features"
  # title = NULL
  
  #library(magrittr)
  
  #check if is x is XCMSnExp (xcms) or featureGroups (patRoon)
  if (base::class(x) == "featureGroupsXCMS3") {
    if (!base::is.null(fileIndex)) {x <- x[fileIndex,]}
    y <- x@xdata
  } else {
    if (base::class(x) == "XCMSnExp") {
      if (!base::is.null(fileIndex)) {x <- xcms::filterFile(x, fileIndex, keepFeatures = TRUE, keepAdjustedRtime = TRUE)}
      y <- x
    } else {
      stop("x must be an XCMSnExp (xcms) or featureGroups (patRoon) object.")
    }
  }
  
  rtr <- c(base::min(MSnbase::rtime(y)), base::max(MSnbase::rtime(y)))
  if (!base::is.null(rt)) if (rtUnit == "min") rt <- rt*60
  if (!base::is.null(rtWindow)) if (rtUnit == "min") rtWindow <- rtWindow*60
  if (!base::is.null(rt) & !base::is.null(rtWindow)) { rtr <- c((rt) - rtWindow, (rt) + rtWindow) }
  if (base::is.null(rt)) if (!base::is.null(rtWindow)) if (base::length(rtWindow) == 2) { rtr <- rtWindow }
  
  #When priotity is for feature and than mz
  if (!base::is.null(features))
  {
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
    FT <- base::row.names(xcms::featureDefinitions(y, mz = mz, ppm = ppm, rt = rtr, type = "within", msLevel = 1))
  }
  
  
  defFT <- xcms::featureDefinitions(y)
  defFT <- defFT[base::row.names(defFT) %in% FT,]
  
  features <- base::character()
  for (f in 1:base::length(FT)) {
    features <- c(features, paste("M", base::round(defFT$mzmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                  "_R", base::round(defFT$rtmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                  "_", FT[f],sep=""))
  }
  
  # temp_rtFT <- as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
  # chrom <- MSnbase::chromatogram(y, aggregationFun = "sum", rt = rtr, msLevel = 1, missing = 0,
  #                                mz = c(base::min(temp_rtFT$mzmin-(ppm/1E6*temp_rtFT$mzmin), na.rm = T),
  #                                       base::max(temp_rtFT$mzmax+(ppm/1E6*temp_rtFT$mzmax), na.rm = T)))
  chrom <- xcms::featureChromatograms(y, features = FT, aggregationFun = "sum", expandRt = 60, include = "any", filled = TRUE, missing = 0)
  
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
  
  #colors <- ntsIUTA::getColors(y, "samples")
  
  if (base::is.null(title)) title <- if (plotBy == "features") base::ifelse(base::length(FT) == 1, features, "")
  
  title <- base::list(text = title, x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))
  
  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                      titlefont = base::list(size = 12, color = "black"),
                      range = c(rtr[1], rtr[2]), autotick = T, ticks = "outside")
  
  yaxis = base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                     titlefont = list(size = 12, color = "black"))
  
  
  if (plotBy == "samples") {
    colors <- ntsIUTA::getColors(y, "samples")
    legG <- y$sample_name
  } else {
    if (plotBy == "features") {
      colors <- ntsIUTA::getColors(base::length(features))
      legG <- features
    } else {
      colors = ntsIUTA::getColors(y, "groups")
      legG <- y$sample_name
    }
  }
  
  
  # if (plotBy == "replicates") legG <- y$sample_group
  # if (plotBy == "samples") legG <- y$sample_name
  
  #By feature the eic followed by the integrated are, legend grouped as function input plotBy
  plot <- plotly::plot_ly(df)
  showlegend = base::rep(0,base::length(features))
  for (s in 1:base::ncol(chrom)) { #base::ncol(feat)
    for (f in 1:base::length(FT)) { #base::nrow(feat)
      showlegend[f] = 1 + showlegend[f]
      plot <- plot %>% plotly::add_trace(df,
                                         x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "rtime"],
                                         y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "intensity"],
                                         type = "scatter", mode = "lines",
                                         line = base::list(width = 0.5, color = base::unname(colors[base::ifelse(plotBy == "features", f, s)])),
                                         connectgaps = TRUE,
                                         name = legG[base::ifelse(plotBy == "features", f, s)],
                                         legendgroup = legG[base::ifelse(plotBy == "features", f, s)],
                                         showlegend = base::ifelse(plotBy == "samples", base::ifelse(f == 1, T, F),
                                                                   base::ifelse(plotBy == "features", base::ifelse(showlegend[f] == 1, T, F),
                                                                                base::ifelse(f == 1, T, F))))
      
      rtFT <- base::as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
      if(unique(TRUE %in% (rtFT$sample == s))) {
        rtFT <- rtFT[rtFT$sample == s, c("rtmin","rtmax")]
        plot <- plot %>%  plotly::add_trace(df,
                                            x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2], "rtime"],
                                            y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2],"intensity"],
                                            type = "scatter", mode =  "lines+markers", fill = 'tozeroy', connectgaps = TRUE,
                                            fillcolor = base::paste(color = base::unname(colors[base::ifelse(plotBy == "features", f, s)]),50, sep = ""),
                                            line = base::list(width = 0.1, color = base::unname(colors[base::ifelse(plotBy == "features", f, s)])),
                                            marker = base::list(size = 3, color = base::unname(colors[base::ifelse(plotBy == "features", f, s)])),
                                            name = legG[base::ifelse(plotBy == "features", f, s)],
                                            legendgroup = legG[base::ifelse(plotBy == "features", f, s)],
                                            showlegend = F,
                                            hoverinfo = 'text', text = base::paste('</br> feature: ', features[f],
                                                                                   '</br> sample: ', y$sample_name[s]))
      }
    }
  }
  # showlegend = 0
  # for (s in 1:base::ncol(chrom)) { #base::ncol(feat)
  #   for (f in 1:base::nrow(chrom)) { #base::nrow(feat)
  #     rtFT <- base::as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
  #     if(base::unique(TRUE %in% (rtFT$sample == s))) {
  #       showlegend <- 1 + showlegend
  #       rtFT <- rtFT[rtFT$sample == s, c("rtmin","rtmax")]
  #       plot <- plot %>%  plotly::add_trace(df,
  #                                           x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2], "rtime"],
  #                                           y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2],"intensity"],
  #                                           type = "scatter", mode = "lines+markers", fill = 'tozeroy', connectgaps = TRUE, fillcolor = paste(color = base::unname(colors[s]),50, sep = ""),
  #                                           line = base::list(width = 0.1, color = base::unname(colors[s])),
  #                                           marker = base::list(size = 3, color = base::unname(colors[s])),
  #                                           name = features[f], legendgroup = FT[f], showlegend = base::ifelse(showlegend == 1, T, F))
  #     }
  #   }
  # }
  
  plot <- plot %>% plotly::layout(xaxis = xaxis,yaxis = yaxis, title = title) #legend = list(title = list(text='<b> Sample: </b>'))
  
  return(plot)
  
}




