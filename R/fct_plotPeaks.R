

#' @title mapPeaks
#' @description Method for plotting peaks from a \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param samples The index or name of the sample/s.
#' The default is \code{NULL} and all samples are used.
#' @param ID The identifier of the peaks of interest.
#' When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.
#' @param mz Optional target \emph{m/z} to find peaks using
#' the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the peaks
#' when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect peaks.
#' The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}.
#' The default is \code{min}.
#' @param colorBy Possible values are \code{"peaks"} (the default),
#' \code{"samples"} or \code{sampleGroups},
#' for colouring by peaks, samples or sample replicate groups, respectively.
#'
#' @return A peak/s map plot produced through \pkg{plotly}.
#'
#' @export
#'
#' @importFrom checkmate assertClass assertSubset
#' @importFrom dplyr between filter
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
mapPeaks <- function(obj, samples = NULL,
                      ID = NULL,
                      mz = NULL, ppm = 20,
                      rt = NULL, rtWindow = NULL,
                      rtUnit = "sec",
                      colorBy = "samples") {
  
  assertClass(obj, "ntsData")
  
  assertSubset(rtUnit, c("sec", "min"))
  
  assertSubset(colorBy, c("peaks", "samples", "sampleGroups"))
  
  if (!is.null(samples)) obj <- filterFileFaster(obj, samples)
  
  if (!is.null(ID)) {
    pki <- obj@peaks[obj@peaks$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      if (is.null(rtr)) rtr <- c(min(obj@peaks$rtmin), max(obj@peaks$rtmax))
      pki <- dplyr::filter(obj@peaks,
                          dplyr::between(mz, mzr[1], mzr[2]),
                          dplyr::between(rt, rtr[1], rtr[2]))
    } else {
      return(cat("One of ID or mz should be given."))
    }
  }
  
  if (nrow(pki) == 0) return(cat("No features found."))
  
  if (is.null(rtr)) rtr <- c(min(pki$rtmin) - 60, max(pki$rtmax) + 60)
  
  if (is.null(ppm)) ppm <- 5
  
  sp <- obj@samples$sample
  
  rg <- obj@samples$group
  
  if (colorBy == "samples") {
    colors <- getColors(obj, "samples")
    col_val <- pki$sample
  } else {
    if (colorBy == "peaks") {
      colors <- getColors(nrow(pki))
      col_val <- pki$ID
    } else {
      colors <- getColors(obj, "sampleGroups")
      col_val <- pki$group
    }
  }
  
  plot <- plot_ly()
  
  plot <- plot %>% add_trace(data = pki, x = pki$rt, y = pki$mz, color = pki$sample, 
                             type = 'scatter', mode = 'markers', colors = colors,
                             marker = list(size = 8),
                             hoverinfo = "text",
                             text = paste("</br> peak: ", pki$ID,
                                          "</br> sample: ", pki$sample,
                                          "</br> <i>m/z</i>: ", round(pki$mz, digits = 4),
                                          "</br> dppm: ", round(((pki$mzmax - pki$mzmin) / pki$mz) * 1E6, digits = 0),
                                          "</br> rt: ", round(pki$rt, digits = 0),
                                          "</br> drt: ", round(pki$rtmax - pki$rtmin, digits = 0),
                                          "</br> Int: ", round(pki$intensity, digits = 0),
                                          "</br> Filled: ",
                                          if ("is_filled" %in% colnames(pki)) {
                                            ifelse(pki$is_filled == 1, TRUE, FALSE)
                                          } else {
                                            FALSE
                                          }))
  
  shapes <- list()
  
  for (i in seq_len(nrow(pki))) {  
    shapes[[i]] <- list(type = "rect",
                        fillcolor = colors[names(colors) %in% pki$sample[i]],
                        opacity = 0.2,
                        line = list(color = colors[names(colors) %in% pki$sample[i]]),
                        x0 = pki$rtmin[i],
                        x1 = pki$rtmax[i],
                        xref = "x",
                        y0 = pki$mzmin[i],
                        y1 = pki$mzmax[i],
                        yref = "y")
  }
  
  dppm <- c(round(min(pki$mzmin), digits = 4), round(max(pki$mzmax), digits = 4))
  dppm_val <- round((dppm[2] - dppm[1])/dppm[2] * 1E6, digits = 0)
  drt <- c(round(min(pki$rtmin), digits = 0), round(max(pki$rtmax), digits = 0))
  
  title_text <- paste0("<i>m/z</i>",": ", dppm[1], " - ", dppm[2],
                       " (", dppm_val, " ppm)",
                       "  rt: ", drt[1], " - ", drt[2],
                       " (",  drt[2] - drt[1], " sec.)")
  
  title <- list(text = title_text, x = 0.1, y = 0.98,
                font = list(size = 9, color = "black"))
  
  xaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Retention Time (sec.)",
                titlefont = list(size = 12, color = "black"),
                range = c(rtr[1], rtr[2]),
                autotick = TRUE, ticks = "outside")
  
  yaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "<i>m/z</i>",
                titlefont = list(size = 12, color = "black"))
  
  plot <- plot %>% plotly::layout(xaxis = xaxis,
                                  yaxis = yaxis,
                                  title = title,
                                  shapes = shapes)
  
  return(plot)

}


#' @title plotPeaks
#' @description Method for plotting peaks from a \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param samples The index or name of the sample/s.
#' The default is \code{NULL} and all samples are used.
#' @param ID The identifier of the peaks of interest.
#' When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.
#' @param mz Optional target \emph{m/z} to find peaks using
#' the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the peaks
#' when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect peaks.
#' The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}.
#' The default is \code{min}.
#' @param colorBy Possible values are \code{"peaks"} (the default),
#' \code{"samples"} or \code{sampleGroups},
#' for colouring by peaks, samples or sample replicate groups, respectively.
#'
#' @return A peak/s map plot produced through \pkg{plotly}.
#'
#' @export
#'
#' @importFrom checkmate assertClass assertSubset
#' @importFrom dplyr between filter
#' @importFrom plotly toRGB plot_ly add_trace layout add_segments
#'
plotPeaks <- function(obj, samples = NULL,
                     ID = NULL,
                     mz = NULL, ppm = 20,
                     rt = NULL, rtWindow = NULL,
                     rtUnit = "sec",
                     colorBy = "samples") {
  
  assertClass(obj, "ntsData")
  
  assertSubset(rtUnit, c("sec", "min"))
  
  assertSubset(colorBy, c("peaks", "samples", "sampleGroups"))
  
  if (!is.null(samples)) obj <- filterFileFaster(obj, samples)
  
  if (!is.null(ID)) {
    pki <- obj@peaks[obj@peaks$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      if (is.null(rtr)) rtr <- c(min(obj@peaks$rtmin), max(obj@peaks$rtmax))
      pki <- dplyr::filter(obj@peaks,
                           dplyr::between(mz, mzr[1], mzr[2]),
                           dplyr::between(rt, rtr[1], rtr[2]))
    } else {
      return(cat("One of ID or mz should be given."))
    }
  }
  
  if (nrow(pki) == 0) return(cat("No features found."))
  
  if (is.null(rtr)) rtr <- c(min(pki$rtmin) - 60, max(pki$rtmax) + 60)
  
  if (is.null(ppm)) ppm <- 5
  
  sp <- obj@samples$sample
  
  rg <- obj@samples$group
  
  if (colorBy == "samples") {
    colors <- getColors(obj, "samples")
    col_val <- apply(pki, MARGIN = 1, FUN = function(x) colors[names(colors) %in% x["sample"]])
    leg <- pki$sample
    sleg <- !duplicated(leg)
  } else {
    if (colorBy == "peaks") {
      colors <- getColors(nrow(pki))
      col_val <- colors
      leg <- paste0( pki$ID, "/", round(pki$mz, digits = 4), "/", round(pki$rt, digits = 0))
      sleg <- !duplicated(leg)
    } else {
      colors <- getColors(obj, "sampleGroups")
      col_val <- apply(pki, MARGIN = 1, FUN = function(x) colors[names(colors) %in% x["sample"]])
      leg <- pki$group
      sleg <- !duplicated(leg)
    }
  }
    
  EIC <- extractEIC(obj,
                    mz = c(min(pki$mzmin), max(pki$mzmax)),
                    rtWindow = rtr, rtUnit = "sec")
  
  plot <- plot_ly()
  
  for (i in seq_len(nrow(pki))) {
    
    df <- EIC[EIC$mz >= pki$mzmin[i] &
                EIC$mz <= pki$mzmax[i] &
                EIC$file == which(sp == pki$sample[i]), ]
    
    
    
    plot <- plot %>% add_trace(df,
                               x = df$rt,
                               y = df$i,
                               type = "scatter", mode = "markers",
                               marker = list(size = 0.2,
                                             color = col_val[i]),
                               connectgaps = TRUE,
                               name = leg[i],
                               legendgroup = leg[i],
                               showlegend = sleg[i]
    )
    
    df <- df[df$rt >= pki$rtmin[i] & df$rt <= pki$rtmax[i], ]
    
    plot <- plot %>%  add_trace(df,
                                x = df$rt,
                                y = df$i,
                                type = "scatter", mode =  "markers",
                                fill = "tozeroy", connectgaps = TRUE,
                                fillcolor = paste(color = col_val[i], 50, sep = ""),
                                #line = list(width = 0.1, color = col_val[i]),
                                marker = list(size = 3, color = col_val[i]),
                                name = leg[i],
                                legendgroup = leg[i],
                                showlegend = FALSE,
                                hoverinfo = "text",
                                text = paste("</br> peak: ", pki$ID[i],
                                             "</br> sample: ", pki$sample[i],
                                             "</br> <i>m/z</i>: ", round(pki$mz[i], digits = 4),
                                             "</br> dppm: ", round(((pki$mzmax[i] - pki$mzmin[i]) / pki$mz[i]) * 1E6, digits = 0),
                                             "</br> rt: ", round(pki$rt[i], digits = 0),
                                             "</br> drt: ", round(pki$rtmax[i] - pki$rtmin[i], digits = 0),
                                             "</br> Int: ", round(pki$intensity[i], digits = 0),
                                             "</br> Filled: ",
                                             if ("is_filled" %in% colnames(pki)) {
                                               ifelse(pki$is_filled[i] == 1, TRUE, FALSE)
                                             } else {
                                               FALSE
                                             }))
    
    plot <- plot %>% add_segments(x = pki$rt[i], xend = pki$rt[i], y = 0, yend = pki$intensity[i],
                                  legendgroup = leg[i], showlegend = FALSE, line = list(color = col_val[i], size = 0.5))
    
  }
  
  dppm <- c(round(min(pki$mzmin), digits = 4), round(max(pki$mzmax), digits = 4))
  dppm_val <- round((dppm[2] - dppm[1])/dppm[2] * 1E6, digits = 0)
  drt <- c(round(min(pki$rtmin), digits = 0), round(max(pki$rtmax), digits = 0))
  
  title_text <- paste0("<i>m/z</i>",": ", dppm[1], " - ", dppm[2],
                       " (", dppm_val, " ppm)",
                       "  rt: ", drt[1], " - ", drt[2],
                       " (",  drt[2] - drt[1], " sec.)")
  
  title <- list(text = title_text, x = 0.1, y = 0.98,
                font = list(size = 9, color = "black"))
  
  xaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Retention Time (sec.)",
                titlefont = list(size = 12, color = "black"),
                range = c(rtr[1], rtr[2]),
                autotick = TRUE, ticks = "outside")
  
  yaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Intensity",
                titlefont = list(size = 12, color = "black"))
  
  plot <- plot %>% plotly::layout(xaxis = xaxis,
                                  yaxis = yaxis,
                                  title = title)
  
  return(plot)
  
  # TODO Idea for 3D plotly
  # if (td) {
  #   
  #   df <- EIC
  #   df$file <- sapply(df$file,FUN = function(x) sp[x])
  #   df$file <- factor(df$file, levels = sp, ordered = FALSE)
  #   
  #   plot_ly(df, x = ~rt, y = ~mz, z = ~i, color = ~file, colors = c('#BF382A', '#0C4B8E'),
  #           marker = list(size = 2), mode = "scatter3d")
  # }
  
}
