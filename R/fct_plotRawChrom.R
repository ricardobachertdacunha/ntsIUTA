

#' @title plotRawChrom
#' @description A method for plotting total, base and extracted ion chromatograms
#' (TIC, BPC and EIC, respectively) of an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object with one or more files.
#' @param fileIndex The index of the file/s to extract the data.
#' @param mz Optional target \emph{m/z} to obtain an EIC.
#' Note that when not \code{NULL} (the default), EIC is always returned.
#' @param ppm The mass deviation to extract the data for the EIC, in \code{ppm}.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time window or deviation to collect the data.
#' The time unit is defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' A vector of length 1 is assumed as a deviation of a given \code{rt}.
#' @param rtUnit Possible entries are \code{sec} (the default) or \code{min}.
#' @param msLevel The MS level to extract the data.
#' For the moment, only 1 is possible.
#' @param type The type of chromatogram.
#' Possible entries are "bpc" for base peak chromatogram
#' or "tic" for total ion chromatogram.
#' The default is "tic". If \code{mz} is specified (not \code{NULL}),
#' the type is set automatically to EIC.
#' @param colorBy Possible values are \code{"samples"} or \code{samplegroups}
#' (the default), for colouring by samples or sample replicate groups respectively.
#' @param interactive Logical, set to \code{TRUE} to use
#' the \pkg{plotly} instead of \pkg{ggplot2}. The default is \code{FALSE}.
#'
#' @return A plot for inspection of the raw data
#' in the given \linkS4class{ntsData} object.
#'
#' @export
#'
#' @importFrom plotly toRGB plot_ly add_trace layout
#' @importFrom dplyr group_by arrange top_n summarize
#' @importFrom checkmate checkSubset
#'
setMethod("plotRawChrom", "ntsData", function(obj, fileIndex = NULL,
                                              mz = NULL, ppm = 20,
                                              rt = NULL, rtWindow = NULL,
                                              rtUnit = "sec",
                                              msLevel = 1,
                                              type = "tic",
                                              colorBy = "samplegroups",
                                              interactive = FALSE) {
  
  # obj = dtcent
  # fileIndex = 1:2
  # mz <- 213.1869
  # rt <- NULL
  # rtUnit = "min"
  # ppm <- 20
  # rtWindow = NULL
  # msLevel = 1
  # type = "tic"
  # colorBy = "samplegroups"
  # interactive = FALSE
  
  checkmate::checkSubset(rtUnit, c("sec", "min"))
  checkmate::checkSubset(type, c("tic", "bpc"))
  checkmate::checkSubset(colorBy, c("samples", "samplegroups"))
  
  if (!is.null(mz) && length(mz) == 1) {
    if (!is.null(ppm)) ppm <- 20
    main <- paste0("EIC of ", round(mz, digits = 4), " +/- ", round(ppm, digits = 0), " ppm")
    type <- "eic"
  } else {
    if (!is.null(mz) && length(mz) == 2) {
      main <- paste0("EIC for ", round(mz[1], digits = 4), "to ", round(mz[2], digits = 4))
      type <- "eic"
    } else {
      main <- toupper(type)
    }
  }
  
  if (!is.null(fileIndex)) obj <- filterFileFaster(obj, fileIndex)
  
  df <- extractEIC(obj = obj,
                   fileIndex = NULL,
                   mz = mz, ppm = ppm,
                   rt = rt, rtWindow = rtWindow,
                   rtUnit = rtUnit, msLevel = 1,
                   normIntensity = FALSE)
  
  if (type == "tic") {
    mz <- df %>% group_by(rt) %>% top_n(1, i)
    mz <- as.data.frame(mz)
    mz <- arrange(mz, rt)
    y <- df %>% group_by(rt) %>% summarize(i = sum(i))
    y <- as.data.frame(y)
    y <- arrange(y, rt)
    mz$i <- y[, "i", drop = TRUE]
    df <- mz
  }
  
  if (type == "bpc") {
    y <- df %>% group_by(rt) %>% top_n(1, i)
    y <- as.data.frame(y)
    df <- arrange(y, rt)
  }
  
  for (i in seq_len(length(unique(df$file)))) {
    df[df$file == i, "file"] <- samples(obj)[i]
  }
  
  cl <- getColors(obj, which = colorBy)
  
  if (!interactive) {
    
    plot <- ggplot(data = df, aes(x = rt, y = i, color = file)) +
      geom_line(size = 0.5) +
      scale_color_manual(values = cl) +
      ggtitle(main) +
      theme_bw() +
      ylab("Intensity") +
      xlab("Retention Time") +
      theme(legend.title = element_blank())
    
  } else {
    
    title <- list(text = main, x = 0.1, y = 0.98, font = list(size = 14, color = "black"))
    
    xaxis <- list(linecolor = toRGB("black"),
                  linewidth = 2, title = "Retention Time (sec.)",
                  titlefont = list(size = 12, color = "black"))
    
    yaxis <- list(linecolor = toRGB("black"),
                  linewidth = 2, title = "Intensity",
                  titlefont = list(size = 12, color = "black"))
    
    plot <- plot_ly(df,
                    x = df[df$file == samples(obj)[1], "rt"],
                    y = df[df$file == samples(obj)[1], "i"],
                    type = "scatter", mode = "lines+markers",
                    line = list(width = 0.5, color = unname(cl[1])),
                    marker = list(size = 2, color = unname(cl[1])),
                    name = samples(obj)[1],
                    text = round(df[df$file == samples(obj)[1], "mz"], digits = 4),
                    hovertemplate = paste('<i>m/z</i>: %{text}',
                                          '<br>rt: %{x}<br>',
                                          'Int: %{y}'))
    
    if (length(unique(df$file)) > 1) {
      for (i in 2:length(unique(df$file))) {
        plot  <- plot %>% add_trace(df,
                                    x = df[df$file == samples(obj)[i], "rt"],
                                    y = df[df$file == samples(obj)[i], "i"],
                                    type = "scatter", mode = "lines+markers",
                                    line = list(width = 0.5, color = unname(cl[i])),
                                    marker = list(size = 2, color = unname(cl[i])),
                                    name = samples(obj)[i],
                                    text = round(df[df$file == samples(obj)[i], "mz"], digits = 4),
                                    hovertemplate = paste('<i>m/z</i>: %{text}',
                                                          '<br>rt: %{x}<br>',
                                                          'Int: %{y}'))
      }
    }
    
    plot <- plot %>% layout(legend = list(title = list(text = "<b> Sample: </b>")),
                            xaxis = xaxis, yaxis = yaxis, title = title)
    
  }
  
  return(plot)
  
})
