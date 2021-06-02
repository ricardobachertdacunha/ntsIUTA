

#' @title plotRawChrom
#' @description Plots total, base and extracted ion chromatograms (TIC, BPC and EIC, respectively) of an \linkS4class{OnDiskMSnExp} object,
#' using the \pkg{plotly} package to procude an iteractive plot.
#'
#' @param raw An \linkS4class{OnDiskMSnExp} object with one or more files.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param mz Optional target \emph{m/z} to obtain an extracted ion chromatogram (EIC).
#' When not \code{NULL} (the default), EIC is always returned.
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{sec}.
#' @param msLevel The MS level to extract the data. For the momment, only 1 is possible.
#' @param type The type of chromatogram. Possible entries are "bpc" for base peak chromatogram or "tic" for total ion chromatogram.
#' The default is "tic". If \code{mz} is specified (not \code{NULL}), the type is set automatically to EIC.
#'
#' @return An iterative plot for inspection of the raw data in the given \linkS4class{OnDiskMSnExp} object.
#'
#' @export
#'
#' @importMethodsFrom MSnbase filterFile
#' @importFrom plotly toRGB plot_ly add_trace layout
#' @importFrom dplyr group_by arrange top_n summarize
#'
#' @examples
#'
plotRawChrom <- function(raw = rawData, fileIndex = NULL,
                         mz = NULL, ppm = 20,
                         rt = NULL, rtWindow = NULL,
                         rtUnit = "sec",
                         msLevel = 1, type = "tic") {

  raw = rawData
  fileIndex = NULL
  mz <- c(233.0243)
  rt <- NULL
  rtUnit = "min"
  ppm <- 20
  rtWindow = NULL
  msLevel = 1

  if (!is.null(mz) && length(mz) == 1) {
    if (!is.null(ppm)) ppm <- 20
    main <- paste0("EIC of ", round(mz, digits = 4), " +/- ", round(ppm, digits = 0), " ppm")
    type <- "eic"
  } else {
    main <- toupper(type)
  }

  if (!is.null(fileIndex)) {
    raw <- filterFile(raw, fileIndex)
  }

  colors <- getColors(raw, "samples")

  df <- extractEIC(raw = raw,
                   fileIndex = NULL,
                   mz = mz, ppm = ppm,
                   rt = rt, rtWindow = rtWindow,
                   rtUnit = rtUnit, msLevel = 1,
                   normIntensity = FALSE)

  if (type == "tic") {
    mz <- df %>% dplyr::group_by(rt) %>% dplyr::top_n(1, i)
    mz <- base::as.data.frame(mz)
    mz <- dplyr::arrange(mz, rt)
    y <- df %>% dplyr::group_by(rt) %>% dplyr::summarize(i = sum(i))
    y <- base::as.data.frame(y)
    y <- dplyr::arrange(y, rt)
    mz$i <- y[, "i", drop = TRUE]
    df <- mz
  }

  if (type == "bpc") {
    y <- df %>% dplyr::group_by(rt) %>% dplyr::top_n(1, i)
    y <- base::as.data.frame(y)
    df <- dplyr::arrange(y, rt)
  }

  for (i in base::seq_len(base::length(base::unique(df$file)))) {
    df[df$file == i, "file"] <- raw$sample_name[i]
  }

  title <- base::list(text = main, x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))

  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                      titlefont = base::list(size = 12, color = "black"))

  yaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                     titlefont = base::list(size = 12, color = "black"))

  plot <- plotly::plot_ly(df,
                          x = df[df$file == raw$sample_name[1], "rt"],
                          y = df[df$file == raw$sample_name[1], "i"],
                          type = "scatter", mode = "lines+markers",
                          line = base::list(width = 0.5, color = base::unname(ntsIUTA::getColors(raw, "samples")[1])),
                          marker = base::list(size = 2, color = base::unname(ntsIUTA::getColors(raw, "samples")[1])),
                          name = raw$sample_name[1])

  for (i in 2:base::length(base::unique(df$file))) {
    plot  <- plot %>% plotly::add_trace(df,
                                        x = df[df$file == raw$sample_name[i], "rt"],
                                        y = df[df$file == raw$sample_name[i], "i"],
                                        type = "scatter", mode = "lines+markers",
                                        line = base::list(width = 0.5, color = base::unname(ntsIUTA::getColors(raw, "samples")[i])),
                                        marker = base::list(size = 2, color = base::unname(ntsIUTA::getColors(raw, "samples")[i])),
                                        name = raw$sample_name[i])
  }

  plot <- plot %>% plotly::layout(legend = base::list(title = base::list(text = "<b> Sample: </b>")),
                                  xaxis = xaxis, yaxis = yaxis, title = title)

  return(plot)

}
