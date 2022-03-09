

#' @title plotFeaturePeaks
#' 
#' @description Plots peaks for each feature in an \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
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
plotFeaturePeaks <- function(object,
                             samples = NULL,
                             targets = NULL,
                             mz = NULL, ppm = 20,
                             rt = NULL, sec = 30,
                             legendNames = NULL) {

  checkmate::assertClass(object, "ntsData")

  fts <- features(
    object,
    samples,
    targets,
    mz, ppm,
    rt, sec
  )

  pks <- peaks(
    object,
    samples,
    targets = fts$id
  )

  eic <- lapply(split(pks, pks$sample), function(x, object) {
    return(
      extractEICs(
        object,
        samples = unique(x$sample),
        mz = x
      )
    )
  }, object = object)

  eic <- rbindlist(eic)

  if (nrow(eic) < 1) return(cat("Data was not found for any of the targets!"))

  eic$var <- sapply(eic$id, function(x, pks) {
    pks[id == x, feature]
  }, pks = pks)

  if (!is.null(legendNames) & length(legendNames) == length(unique(eic$var))) {
    leg <- legendNames
    names(leg) <- unique(eic$var)
    eic$var <- sapply(eic$var, function(x) leg[x])
  } else {
    leg <- unique(eic$var)
    names(leg) <- unique(eic$var)
  }

  colors <- getColors(leg)

  showleg <- rep(TRUE, length(leg))
  names(showleg) <- names(leg)

  plot <- plot_ly()

  for (i in fts$id) {

    pk_temp <- pks[feature == i, ]

    for (z in pk_temp$id) {

      df <- eic[id == z, ]

      plot <- plot %>% add_trace(df,
        x = df$rt,
        y = df$intensity,
        type = "scatter", mode = "lines",
        line = list(width = 0.5,
                    color = colors[i]),
        connectgaps = TRUE,
        name = leg[i],
        legendgroup = leg[i],
        showlegend = FALSE
      )

      df <- df[rt >= pk_temp[id == z, rtmin] & rt <= pk_temp[id == z, rtmax], ]

      plot <- plot %>%  add_trace(
        df,
        x = df$rt,
        y = df$intensity,
        type = "scatter", mode =  "lines+markers",
        fill = "tozeroy", connectgaps = TRUE,
        fillcolor = paste(color = colors[i], 50, sep = ""),
        line = list(width = 0.1, color = colors[i]),
        marker = list(size = 3, color = colors[i]),
        name = leg[i],
        legendgroup = leg[i],
        showlegend = showleg[i],
        hoverinfo = "text",
        text = paste(
          "</br> name: ", leg[i],
          "</br> feature: ", i,
          "</br> peak: ", z,
          "</br> sample: ", pk_temp[id == z, sample],
          #"</br> <i>m/z</i>: ", round(df$mz, digits = 4),
          "</br> rt: ", round(df$rt, digits = 0),
          "</br> Int: ", round(df$intensity, digits = 0)
        )
      )

      showleg[i] <- FALSE
    }
  }

  plot2 <- plot_ly()

  for (i in fts$id) {

    df2 <- pks[feature == i, ]

    df_p <- df2[is_filled == 0, ]

    plot2 <- plot2 %>% add_trace(
      x = df_p$rt,
      y = df_p$sample,
      type = "scatter",
      mode = "markers",
      marker = list(
        line = list(color = colors[i], width = 3),
        color = "#000000", size = 10
      ),
      error_x = list(
        type = "data",
        symmetric = FALSE,
        arrayminus = df_p$rt - df_p$rtmin,
        array = df_p$rtmax - df_p$rt,
        color = colors[i],
        width = 5
      ),
      name = leg[i],
      legendgroup = leg[i],
      showlegend = FALSE,
      hoverinfo = "text",
      text = paste(
        "</br> name: ", leg[i],
        "</br> feature: ", i,
        "</br> peak: ", df_p$id,
        "</br> sample: ", df_p$sample,
        "</br> height: ", round(df_p$intensity, digits = 0),
        "</br> width: ", round(df_p$rtmax - df_p$rtmin, digits = 0),
        "</br> dppm: ", round(((df_p$mzmax - df_p$mzmin) / df_p$mz) * 1E6, digits = 1),
        "</br> filled: ", ifelse(df_p$is_filled == 1, "TRUE", "FALSE")
      )
    )

    df_f <- df2[is_filled == 1, ]

    if (nrow(df_f) > 0) {
      plot2 <- plot2 %>% add_trace(
        x = df_f$rt,
        y = df_f$sample,
        type = "scatter",
        mode = "markers",
        marker = list(
          line = list(color = colors[i], width = 3),
          color = "#f8f8f8",
          size = 10
        ),
        error_x = list(
          type = "data",
          symmetric = FALSE,
          arrayminus = df_f$rt - df_f$rtmin,
          array = df_f$rtmax - df_f$rt,
          color = colors[i],
          width = 5
        ),
        name = leg[i],
        legendgroup = leg[i],
        showlegend = FALSE,
        hoverinfo = "text",
        text = paste(
          "</br> name: ", leg[i],
          "</br> feature: ", i,
          "</br> peak: ", df_f$id,
          "</br> sample: ", df_f$sample,
          "</br> height: ", round(df_f$intensity, digits = 0),
          "</br> width: ", round(df_f$rtmax - df_f$rtmin, digits = 0),
          "</br> dppm: ", round(((df_f$mzmax - df_f$mzmin) / df_f$mz) * 1E6, digits = 1),
          "</br> filled: ", ifelse(df_f$is_filled == 1, "TRUE", "FALSE")
        )
      )
    }
  }

  plot2 <- hide_colorbar(plot2)

  plotList <- list()

  plotList[["plot"]] <- plot

  plotList[["plot2"]] <- plot2

  xaxis <- list(linecolor = toRGB("black"), linewidth = 2,
                title = "Retention Time (sec.)",
                titlefont = list(size = 12, color = "black"),
                range = c(min(eic$rt), max(eic$rt)), autotick = TRUE, ticks = "outside")

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

  return(plotf)
}
