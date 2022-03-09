

#' @title plotAlignment
#'
#' @description Plots the results from the retention time alignment across samples.
#'
#' @param object An \linkS4class{ntsData} object with adjusted retention time.
#' @param samples A numeric or character vector with the sample to plot the
#' retention time adjustment.
#' @param interactive Logical, set to \code{TRUE} for plotting with \pkg{plotly} package.
#'
#' @return A plot with the retention time alignment differences for each sample.
#'
#' @export
#'
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
plotAlignment <- function(object, samples = NULL, interactive = TRUE) {

  if (hasAdjustedRetentionTime(object)) {

    if (is.character(samples)) {
      if (FALSE %in% (samples %in% samples(object))) {
        warning("Given sample names not found in the ntsData object!")
        return(data.table())
      }
      samples <- which(samples(object) %in% samples)
    }

    if (is.null(samples)) samples <- seq_len(length(samples(object)))

    scans <- object@scans[samples]

    colors <- getColors(samples(object)[samples])

    if (interactive) {

      xaxis <- list(
        linecolor = toRGB("black"),
        linewidth = 2,
        title = "Retention time (seconds)",
        titlefont = list(size = 12, color = "black")
      )

      yaxis <- list(
        linecolor = toRGB("black"),
        linewidth = 2,
        title = "RT<sub>Adjusted</sub> - RT<sub>Raw</sub> (seconds)",
        titlefont = list(size = 12, color = "black")
      )

      plot <- plot_ly()

      for (i in samples) {

        plot  <- plot %>% add_trace(
          x = scans[[i]]$retentionTime,
          y = scans[[i]]$adjustment,
          type = "scatter",
          mode = "lines",
          line = list(
            shape = "spline", width = 0.5,
            color = colors[i]
          ),
          name = samples(object)[i],
          legendgroup = samples(object)[i],
          showlegend = TRUE
        )

        df <- scans[[i]][!is.na(adjPoints), ]

        plot <- plot %>% add_trace(
          x = df$adjPoints,
          y = df$adjustment,
          type = "scatter",
          mode = "markers",
          marker = list(
            size = 5,
            color = colors[i]
          ),
          name = samples(object)[i],
          legendgroup = samples(object)[i],
          showlegend = FALSE
        )
      }

      plot <- plot %>% layout(
        legend = list(title = list(text = "<b> Sample: </b>")),
        xaxis = xaxis, yaxis = yaxis
      )

      return(plot)
    }

  } else {

    warning("Alignment information not found!")
    return(NULL)
  }
}
