

#' @title plotTargetCentroids
#' @description Plot centroids or profile data for a target compound
#' using expected \emph{m/z}, retention time and respective deviations.
#'
#' @param obj An \linkS4class{ntsData} object with one or more MS files.
#' @param samples The index or names of the sample/s to extract the centroids or profile data.
#' @param mz The target \emph{m/z}. Note that \code{m/z} is not the expected
#' monoisotopic mass but the expected adduct (\emph{i.e.} \code{[M+H]+}).
#' @param ppm The mass deviation to extract centroids in \code{ppm}.
#' The default is set to 20 ppm.
#' @param rt The expected retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time deviation to collect centroids or profile data.
#' The time unit is the defined by \code{rtUnit}.
#' It is recommended a minimum of 1 minute.
#' @param rtUnit Possible entries are \code{min} or \code{sec}.
#' The default is \code{sec}.
#' @param title Optional title for the plot/s.
#' @param plotTargetMark Logical, set to TRUE, the defaulft, to plot a
#' target mark with +/- 5 ppm and +/- 10 seconds deviation.
#'
#' @return An iterative plot of the centroids for the requested
#' target based on given \code{mz}.
#'
#' @export
#'
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importMethodsFrom MSnbase filterFile filterMz filterMsLevel filterRt
#' @importFrom methods as
#' @importFrom grDevices colorRamp
#' @importFrom plotly plot_ly layout add_annotations toRGB subplot hide_colorbar hide_legend
#'
plotTargetCentroids <- function(obj = NULL,
                                samples = NULL,
                                mz = NULL, ppm = 20,
                                rt = NULL, rtWindow = NULL,
                                rtUnit = "min", title = NULL,
                                plotTargetMark = TRUE) {
  
  if (!is.null(samples)) obj <- filterFileFaster(obj, samples)

  if (is.null(mz)) return(cat("Target mz should be defined!"))

  df <- extractEIC(obj = obj,
                   samples = NULL,
                   mz = mz, ppm = ppm,
                   rt = rt, rtWindow = rtWindow,
                   rtUnit = rtUnit, msLevel = 1,
                   normIntensity = FALSE)

  if (rtUnit == "min") if (!is.null(rt)) rt <- rt * 60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow * 60

  if (is.null(title)) {
    fns <- basename(obj@samples$file)
  } else {
    fns <- rep(title, length(obj@samples$file))
  }

  sNames <- obj@samples$sample
  rtmin <- min(df$rt, na.rm = TRUE)
  rtmax <- max(df$rt, na.rm = TRUE)
  mzmin <- min(df$mz, na.rm = TRUE)
  mzmax <- max(df$mz, na.rm = TRUE)
  maxInt <- max(df$i, na.rm = TRUE) * 1.1
  df <- split(df, df$file)

  if (any(unlist(lapply(df, nrow)) > 20000))
    warning("The MS area to be plotted seems rather large. It is suggested",
            " to restrict the data first using 'filterRt' and 'filterMz'. ",
            "See also ?chromatogram and ?Chromatogram for more efficient ",
            "functions to plot a total ion chromatogram or base peak ",
            "chromatogram.",
            immediate = TRUE, call = FALSE)

  colors <- colorRamp(c("#383E47", "#5E8CAA", "#16B9E5", "#16E5C9", "#16E54C"))

  line <- list(type = "line",
               line = list(color = "red", dash = "dash", width = 0.5),
               xref = "x", yref = "y")

  plotList <- list()
  vline1 <- list()
  vline2 <- list()
  hline <- list()
  rect <- list()

  for (s in seq_len(length(df))) {

    temp <- df[[s]]

    if (plotTargetMark) {
      vline1 <- list(x0 = rt, x1 = rt, y0 = 0, y1 = max(temp$i, na.rm = TRUE))
      vline2 <- list(x0 = rt, x1 = rt, y0 = mzmin, y1 = mzmax)
      hline <- list(x0 = min(temp$rt, na.rm = TRUE), x1 =  max(temp$rt, na.rm = TRUE), y0 = mz, y1 = mz)
      rect <- list(type = "rect", fillcolor = "red", line = list(color = "red"), opacity = 0.1,
                  x0 = rt - 10, x1 = rt + 10, xref = "x",
                  y0 = mz - ((5 / 1E6) * mz), y1 = mz + ((5 / 1E6) * mz), yref = "y")
    }

    p1 <- plot_ly(data = temp, x = temp$rt, y = temp$i,
                  type = "scatter", mode = "markers", color = temp$i, colors = colors,
                  marker = list(size = 8, line = list(color = "white", width = 0.5)), name = paste0(s, "p1"))

    if (plotTargetMark) p1 <- p1 %>% plotly::layout(shapes = c(vline1, line))

    p1 <- p1 %>% add_annotations(text = sNames[s], x = 0.05, y = 1, yref = "paper", xref = "paper",
                                 xanchor = "left", yanchor = "bottom", align = "center",
                                 showarrow = FALSE, font = list(size = 14))

    p2 <- plot_ly(temp, x = temp$rt, y = temp$mz,
                  type = "scatter", mode = "markers", color = temp$i, colors = colors,
                  marker = list(size = 8, line = list(color = "white", width = 0.5)), name = paste0(s, "p2"))

    if (plotTargetMark) p2 <- p2 %>% plotly::layout(shapes = list(c(vline2, line), c(hline, line), rect))

    plotList[[paste0("p1", s)]] <- p1
    plotList[[paste0("p2", s)]] <- p2
  }

  plotList <- plotList[order(names(plotList))]

  xaxis <- list(linecolor = toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                titlefont = list(size = 12, color = "black"))

  yaxis1 <- list(linecolor = toRGB("black"), linewidth = 2, title = "Intensity",
                titlefont = list(size = 12, color = "black"), range = c(0, maxInt))

  yaxis2 <- list(linecolor = toRGB("black"), linewidth = 2, title = paste("<i>m/z</i>"),
                titlefont = list(size = 12, color = "black"), range = c(mzmin, mzmax))

  plot <- subplot(plotList, nrows = 2, margin = 0.04, shareX = TRUE, shareY = TRUE, which_layout = "merge")
  plot <- hide_colorbar(plot)
  plot <- hide_legend(plot)
  plot <- plot %>% plotly::layout(xaxis = xaxis, xaxis2 = xaxis,
                                  xaxis3 = xaxis, xaxis4 = xaxis,
                                  xaxis5 = xaxis, xaxis6 = xaxis,
                                  yaxis = yaxis1, yaxis2 = yaxis2)

  return(plot)

}
