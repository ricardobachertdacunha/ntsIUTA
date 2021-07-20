

#' @title plotFeatures
#' @description Method for plotting features from a \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param fileIndex The index or name of the sample/s.
#' The default is \code{NULL} and all samples are used.
#' @param ID The identifier of the features of interest.
#' When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.
#' @param mz Optional target \emph{m/z} to find features using
#' the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the features
#' when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect features.
#' The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}.
#' The default is \code{min}.
#' @param msLevel The MS level to extract the data.
#' For the moment, only 1 is possible.
#' @param colorBy Possible values are \code{"features"},
#' \code{"samples"} or \code{samplegroups} (the default),
#' for colouring by samples or sample replicate groups respectively.
#' @param interactive Logical, set to \code{TRUE} to use
#' the \pkg{plotly} instead of \pkg{ggplot2}. The default is \code{TRUE}.
#'
#' @return A feature/s plot produced through \pkg{gglpot2} for interactive
#' \code{FALSE} or an interactive plot through the \pkg{plotly}.
#' 
#' @export
#'
#' @importFrom dplyr between filter
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
setMethod("plotFeatures", "ntsData", function(obj, fileIndex = NULL,
                                              ID = NULL,
                                              mz = NULL, ppm = 20,
                                              rt = NULL, rtWindow = NULL,
                                              rtUnit = "sec",
                                              msLevel = 1,
                                              colorBy = "features",
                                              interactive = TRUE) {
  
  # obj <- dtxcms2
  # fileIndex <- 2:5
  # ID <- c("M275_R622_2158" , "M210_R586_836")
  # mz <- 213.1869
  # rt <- 15.47
  # rtUnit <- "min"
  # ppm <- NULL
  # rtWindow <- NULL
  # colorBy <- "features"
  # interactive <- FALSE
  
  checkmate::checkSubset(rtUnit, c("sec", "min"))
  checkmate::checkSubset(colorBy, c("features", "samples", "samplegroups"))
  
  if (!is.null(fileIndex)) obj <- filterFileFaster(obj, fileIndex)
  
  rtr <- NULL
  
  if (!is.null(ID)) {
    ft <- obj@features[obj@features$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      ft <- dplyr::filter(obj@features,
                          dplyr::between(mz, mzr[1], mzr[2]),
                          dplyr::between(rt, rtr[1], rtr[2]))
    } else {
        return(cat("One of ID or mz should be given."))
    }
  }
  
  if (nrow(ft) == 0) return(cat("No features found."))
  
  pk <- list()
  for (i in seq_len(nrow(ft))) {
    pk[[ft$ID[i]]] <- obj@peaks[obj@peaks$ID %in% unlist(ft$pIdx[i]), ]
  }
  
  pk <- pk[lapply(pk, nrow) > 0]
  
  if (length(pk) == 0) return(cat("No features found."))
  
  if (is.null(rtr)) rtr <- c(min(ft$rtmin) - 60, max(ft$rtmax) + 60)
  
  if (is.null(ppm)) ppm = 5
  
  EICs <- list()
  EICs <- lapply(pk, function(x) {
    extractEIC(obj,
               mz = c(min(x$mzmin) - ((ppm / 1E6) * min(x$mzmin)),
                      max(x$mzmax) + ((ppm / 1E6) * max(x$mzmax))),
               rtWindow = rtr,
               rtUnit = "sec")
  })
  
  sp <- obj@samples$sample
  rg <- obj@samples$group
  
  if (colorBy == "samples") {
    colors <- getColors(obj, "samples")
    spleg <- lapply(pk, function(x) x$sample)
    spleg2 <- spleg[[1]]
    if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, spleg[[s]])
    legG <- spleg2
    sleg <- !duplicated(spleg2)
  } else {
    if (colorBy == "features") {
      colors <- getColors(nrow(ft))
      legG <- ft$ID
      spleg <- lapply(pk, function(x) x$sample)
      spleg2 <- rep(legG[1], length(spleg[[1]]))
      if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, rep(legG[s], length(spleg[[s]])))
      sleg <- !duplicated(spleg2)
    } else {
      colors <- getColors(obj, "samplegroups")
      spleg <- lapply(pk, function(x) x$group)
      spleg2 <- spleg[[1]]
      if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, spleg[[s]])
      legG <- spleg2
      sleg <- !duplicated(spleg2)
    }
  }
  
  if (interactive) {
    
    plot <- plot_ly()
    
    counter <- 1
    
    for (i in seq_len(nrow(ft))) {
      
      for (z in seq_len(nrow(pk[[i]]))) {
        
        df <- EICs[[i]][EICs[[i]]$file == which(sp == pk[[i]]$sample[z]), ]
        
        plot <- plot %>% add_trace(df,
          x = df$rt,
          y = df$i,
          type = "scatter", mode = "lines",
          line = list(width = 0.5,
                            color = colors[ifelse(colorBy == "features", i, z)]),
          connectgaps = TRUE,
          name = legG[ifelse(colorBy == "features", i, z)],
          legendgroup = legG[ifelse(colorBy == "features", i, z)],
          showlegend = sleg[counter]
        )
        
        df <- df[df$rt >= pk[[i]]$rtmin[z] & df$rt <= pk[[i]]$rtmax[z], ]
        
        plot <- plot %>%  add_trace(df,
          x = df$rt,
          y = df$i,
          type = "scatter", mode =  "lines+markers",
          fill = 'tozeroy', connectgaps = TRUE,
          fillcolor = paste(color = colors[ifelse(colorBy == "features", i, z)],50, sep = ""),
          line = list(width = 0.1, color = colors[ifelse(colorBy == "features", i, z)]),
          marker = list(size = 3, color = colors[ifelse(colorBy == "features", i, z)]),
          name = legG[ifelse(colorBy == "features", i, z)],
          legendgroup = legG[ifelse(colorBy == "features", i, z)],
          showlegend = FALSE,
          hoverinfo = 'text',
          text = paste('</br> feature: ', ft$ID[i],
                       '</br> peak: ', pk[[i]]$ID[z],
                       '</br> sample: ', pk[[i]]$sample[z],
                       '</br> <i>m/z</i>: ', round(df$mz, digits = 4),
                       '</br> rt: ', round(df$rt, digits = 0),
                       '</br> Int: ', round(df$i, digits = 0))
        )
        
        counter <- counter + 1
        
      }
    }
    
    shapes <- list()
    if (colorBy == "features") {
      for (i in seq_len(length(pk))) {
        shapes[[i]] <- list(type = "rect",
                            fillcolor = NULL,
                            opacity = NULL,
                            line = list(color = colors[i],
                                        dash = 'dot',
                                        width = 0.5),
                            x0 = min(pk[[i]]$rtmin),
                            x1 = max(pk[[i]]$rtmax),
                            xref = "x",
                            y0 = 0,
                            y1 = max(pk[[i]]$intensity),
                            yref = "y")
      }
    }
    
    title <- if (colorBy == "features") ifelse(nrow(ft) == 1, ft$ID, "")
    title <- list(text = title, x = 0.1, y = 0.98,
                  font = list(size = 14, color = "black"))
    
    xaxis <- list(linecolor = toRGB("black"),
                  linewidth = 2, title = "Retention Time (sec.)",
                  titlefont = list(size = 12, color = "black"),
                  range = c(rtr[1], rtr[2]),
                  autotick = T, ticks = "outside")
    
    yaxis = list(linecolor = toRGB("black"),
                 linewidth = 2, title = "Intensity",
                 titlefont = list(size = 12, color = "black"))
    
    plot <- plot %>% plotly::layout(xaxis = xaxis, 
                                    yaxis = yaxis,
                                    title = title,
                                    shapes = shapes)
    
    return(plot)
    
  } else {
    
    names(colors) <- legG[1:length(colors)]
    
    plot <- ggplot()
      
      for (i in seq_len(nrow(ft))) {
        
        for (z in seq_len(nrow(pk[[i]]))) {
          
          df <- EICs[[i]][EICs[[i]]$file == which(sp == pk[[i]]$sample[z]), ]
          
          df$col <- legG[ifelse(colorBy == "features", i, z)]
          
          plot <- plot + geom_line(data = df, aes(x = rt, y = i, color = col),
                                   size = 0.5)
          
          df <- df[df$rt >= pk[[i]]$rtmin[z] & df$rt <= pk[[i]]$rtmax[z], ]
          
          plot <- plot + geom_area(data = df, aes(x = rt, y = i),
                                   fill = colors[ifelse(colorBy == "features", i, z)],
                                   alpha = 0.3)
          
        }
      }
    
    if (colorBy == "features") {
      
      rec <- data.frame(xmin = sapply(pk, function(x) min(x$rtmin)),
                        xmax = sapply(pk, function(x) max(x$rtmax)),
                        ymin = rep(0, length(pk)),
                        ymax = sapply(pk, function(x) max(x$intensity)),
                        col = legG)
      
      plot <- plot + geom_rect(dat = rec, aes(xmin = xmin,
                                   xmax = xmax,
                                   ymin = ymin,
                                   ymax = ymax,
                                   color = col),
                               linetype = "dashed",
                               fill = NA)
      
    }
    
    plot <- plot +
      theme_bw() +
      ylab("Intensity") +
      xlab("Retention Time (sec.)") +
      labs(color = colorBy) +
      scale_color_manual(values = colors, aesthetics = "color")
    
    return(plot)
    
  }
  
})
