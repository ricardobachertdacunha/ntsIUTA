

#' @title plotFeaturePeaks
#' @description Plots peaks for each feature in an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param samples The index or name of the sample/s.
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
#' @param names A character string with names for each feature given in \code{ID}.
#' Note that length should match between \code{names} and \code{ID}.
#' If length does not match \code{names} are not used.
#' @param interactive Logical, set to \code{TRUE} to use
#' the \pkg{plotly} instead of \pkg{ggplot2}. The default is \code{TRUE}.
#'
#' @return A double plot with peak chromatograms on the top part
#' and feature peak groups below.
#'
#' @export
#'
#' @importFrom checkmate assertClass assertSubset
#' @importFrom plotly plot_ly add_trace layout hide_colorbar subplot toRGB
#' @importFrom dplyr filter between
#' @importFrom stats setNames
#'
plotFeaturePeaks <- function(obj, samples = NULL,
                             ID = NULL,
                             mz = NULL, ppm = 20,
                             rt = NULL, rtWindow = NULL,
                             rtUnit = "sec",
                             msLevel = 1,
                             names = NULL,
                             interactive = TRUE) {

  assertClass(obj, "ntsData")

  assertSubset(rtUnit, c("sec", "min"))

  if (!is.null(names)) if (!(length(names) == length(ID))) names <- NULL

  if (!is.null(samples)) obj <- filterFileFaster(obj, samples)

  rtr <- NULL

  if (!is.null(ID)) {
    ft <- obj@features[obj@features$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      if (is.null(rtr)) rtr <- c(min(obj@features$rtmin), max(obj@features$rtmax))
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

  if (is.null(ppm)) ppm <- 5

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

  colors <- getColors(nrow(ft))

  FlegG <- ft$ID

  if (!is.null(names)) FlegG <- names

  if (interactive) {

    plot <- plot_ly()

    for (i in seq_len(nrow(ft))) {

      for (z in seq_len(nrow(pk[[i]]))) {

        df <- EICs[[i]][EICs[[i]]$file == which(sp == pk[[i]]$sample[z]), ]

        plot <- plot %>% add_trace(df,
                                   x = df$rt,
                                   y = df$i,
                                   type = "scatter", mode = "lines",
                                   line = list(width = 0.5,
                                               color = colors[i]),
                                   connectgaps = TRUE,
                                   name = FlegG[i],
                                   legendgroup = FlegG[i],
                                   showlegend = FALSE
        )

        df <- df[df$rt >= pk[[i]]$rtmin[z] & df$rt <= pk[[i]]$rtmax[z], ]

        plot <- plot %>%  add_trace(df,
                                    x = df$rt,
                                    y = df$i,
                                    type = "scatter", mode =  "lines+markers",
                                    fill = "tozeroy", connectgaps = TRUE,
                                    fillcolor = paste(color = colors[i], 50, sep = ""),
                                    line = list(width = 0.1, color = colors[i]),
                                    marker = list(size = 3, color = colors[i]),
                                    name = FlegG[i],
                                    legendgroup = FlegG[i],
                                    showlegend = FALSE,
                                    hoverinfo = "text",
                                    text = paste(ifelse(!is.null(names), "</br> Name: ", ""),
                                                 ifelse(!is.null(names), names[i], ""),
                                                 "</br> feature: ", ft$ID[i],
                                                 "</br> peak: ", pk[[i]]$ID[z],
                                                 "</br> sample: ", pk[[i]]$sample[z],
                                                 "</br> <i>m/z</i>: ", round(df$mz, digits = 4),
                                                 "</br> rt: ", round(df$rt, digits = 0),
                                                 "</br> Int: ", round(df$i, digits = 0))
        )
      }
    }

    plot2 <- plot_ly()

    rect <- list()

    dotsColor <- c("#000000", "#FDFEFE00")

    dotsColor <- stats::setNames(dotsColor, c("0", "1"))

    for (f in seq_len(nrow(ft))) {

      df2 <- pk[[f]]

      plot2 <- plot2 %>% add_trace(df2,
                                   x = df2$rt,
                                   y = df2$sample, type = "scatter", mode = "markers",
                                   error_x = list(type = "data", symmetric = FALSE,
                                                  array = df2$rtmax - df2$rt,
                                                  arrayminus = df2$rt - df2$rtmin,
                                                  color = colors[f],
                                                  width = 5),
                                   color = as.character(df2$is_filled), colors = dotsColor, size = 6,
                                   marker = list(line = list(color = colors[f], width = 3)),
                                   name = FlegG[f],
                                   legendgroup = FlegG[f],
                                   showlegend = TRUE,
                                   hoverinfo = "text", text = paste(ifelse(!is.null(names), "</br> Name: ", ""),
                                                                    ifelse(!is.null(names), names[i], ""),
                                                                    "</br> feature: ", names(pk[f]),
                                                                    "</br> sample: ", df2$sample,
                                                                    "</br> height: ", round(df2$intensity, digits = 0),
                                                                    "</br> width: ", round(df2$rtmax - df2$rtmin, digits = 0),
                                                                    "</br> dppm: ", round(((df2$mzmax - df2$mzmin) / df2$mz) * 1E6, digits = 1),
                                                                    "</br> filled: ", ifelse(df2$is_filled == 1, "TRUE", "FALSE")))
    }

    plot2 <- hide_colorbar(plot2)

    plotList <- list()

    plotList[["plot"]] <- plot

    plotList[["plot2"]] <- plot2

    xaxis <- list(linecolor = toRGB("black"), linewidth = 2,
                  title = "Retention Time (sec.)",
                  titlefont = list(size = 12, color = "black"),
                  range = rtr, autotick = TRUE, ticks = "outside")

    yaxis1 <- list(linecolor = toRGB("black"), linewidth = 2,
                  title = "Intensity",
                  titlefont = list(size = 12, color = "black"))

    yaxis2 <- list(linecolor = toRGB("black"), linewidth = 2,
                  title = "Sample",
                  titlefont = list(size = 12, color = "black"),
                  tick0 = 0, dtick = 1)

    plotf <- subplot(plotList, nrows = 2, margin = 0.04,
                     shareX = TRUE, which_layout = "merge")

    plotf <- plotf %>% layout(xaxis = xaxis, yaxis = yaxis1, yaxis2 = yaxis2)

    plotf

    return(plotf)

  } else {

    # TODO non-interactive plotFeaturePeaks

  }
}
