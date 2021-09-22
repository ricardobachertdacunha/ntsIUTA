

#' @title plotAlignment
#' @description Plots the results from the retention time alignment across samples.
#'
#' @param x The \linkS4class{XCMSnExp} object obtained from the \code{\link{makeFeatures}} function.
#'
#' @return A plot with the retention time alignment differences for each samples.
#'
#' @export
#'
#' @import magrittr
#' @importFrom xcms hasAdjustedRtime
#' @importMethodsFrom MSnbase rtime
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
plotAlignment <- function(x = features) {

  if (x@algorithm == "xcms3") {

    if (xcms::hasAdjustedRtime(x@xdata)) {

      x <- x@xdata

      colors <- ntsIUTA::getColors(x, "samples")

      main <- base::paste0("Time Alignment")

      diffRt <- MSnbase::rtime(x, adjusted = TRUE) - MSnbase::rtime(x, adjusted = FALSE)
      diffRt <- base::split(diffRt, MSnbase::fromFile(x))

      xRt <- MSnbase::rtime(x, adjusted = TRUE, bySample = TRUE)

      title <- base::list(text = main, x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))

      xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                          titlefont = list(size = 12, color = "black"),
                          range = base::range(xRt, na.rm = TRUE))

      yaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "RT<sub>Adjusted</sub> - RT (sec.)",
                          titlefont = base::list(size = 12, color = "black"),
                          range = base::range(diffRt, na.rm = TRUE) * 1.15)

      plot <- plotly::plot_ly(x = xRt[[1]],
                              y = diffRt[[1]],
                              type = "scatter", mode = "lines",
                              line = base::list(shape = "spline", width = 0.5,
                                          color = base::unname(ntsIUTA::getColors(x, "samples")[1])),
                              name = x$sample_name[1])

      if (base::length(x$sample_name) > 1) {
        for (i in 2:base::length(x$sample_name)) {
          plot  <- plot %>% plotly::add_trace(x = xRt[[i]],
                                              y = diffRt[[i]],
                                              type = "scatter", mode = "lines",
                                              line = base::list(shape = "spline", width = 0.5,
                                                          color = base::unname(ntsIUTA::getColors(x, "samples")[i])),
                                              name = x$sample_name[i])
        }
      }

      plot <- plot %>% plotly::layout(legend = base::list(title = base::list(text = "<b> Sample: </b>")),
                                      xaxis = xaxis, yaxis = yaxis, title = title)

      return(plot)

    } else {
      print("Alignment information not found.")
    }

  } else {
    print("Alignment information not found.")
  }

}
