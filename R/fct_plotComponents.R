

#' @title plotComponents
#' @description Plots annotation details for features selected by \code{ID},
#' \emph{m/z} and retention time or component number (\code{comp}) of all
#' or selected sample replicate groups in an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object containing annotated features.
#' @param samples The samples (name or index) to extract the features.
#' The default is \code{NULL} which considers all existing sample replicate groups,
#' excluding any assigned blank replicate groups as listed in the \code{obj}.
#' Note that the features are taken from sample replicate groups, meaning that
#' features will be extracted from sample replicate groups
#' that contain the specified samples.
#' @param ID The identifier/s of selected features. When specified,
#' it overwrites any given \code{mz} or \code{comp} values.
#' @param mz The \emph{m/z} to find features. Not used if \code{ID} is specified.
#' Can be a vector of length 2, defining the mass range to find features.
#' @param ppm The expected mass deviation (in ppm) to search
#' for features of a given \code{mz}.
#' @param rt The expected retention time of the \emph{m/z} of interest,
#' only used if \code{ID} is not specified.
#' @param rtWindow The expected retention time deviation for searching.
#' Can be a vector of length 2, giving the time range to find features.
#' @param rtUnit The time unit used.
#' Possible values are \code{sec} and \code{min}. Default is \code{sec}.
#' @param comp The component number/s to extract features.
#' Only used if both \code{ID} and \code{mz} are \code{NULL}.
#' @param entireComponents Logical, set to \code{TRUE} (The default) to give the
#' all the features in the components represented by the selected features.
#' @param onlyAnnotated Logical, set to \code{TRUE} to return only annotated features.
#' @param onlyRelated Logical, set to \code{TRUE} to return only features that are related,
#' meaning features annotated with the same molecular ion.
#' @param log Logical, set to \code{TRUE} (the default) to plot the log of the peak heights.
#' @param colorBy Possible values are \code{samples} or \code{groups} (the default)
#' for plotting by samples or sample replicate groups, respectively. Note that for
#' \code{groups} the height of the average peaks is averaged.
#'
#' @return A spectrum plot of annotation details for
#' the selected/found features.
#'
#' @note If all (\code{ID}, \code{mz} and \code{comp}) are \code{NULL},
#' the returned \code{data.frame} will contain all the features
#' for each sample replicate group. Additionally, the three logical arguments
#' are applied in the follwoing order: (1) \code{entireComponents},
#' (2) \code{onlyAnnotated} and (3) \code{onlyRelated}.
#'
#' @export
#'
#' @importFrom checkmate assertSubset assertClass
#' @importFrom stringr str_replace_all
#' @importFrom plotly plot_ly add_bars toRGB layout
#'
plotComponents <- function(obj = NULL,
                           samples = NULL,
                           ID = NULL,
                           mz = NULL, ppm = 5,
                           rt = NULL, rtWindow = 1, rtUnit = "sec",
                           comp = NULL,
                           entireComponents = TRUE,
                           onlyAnnotated = FALSE,
                           onlyRelated = TRUE,
                           log = TRUE,
                           colorBy = c("groups")) {

  assertClass(obj, "ntsData")

  assertSubset(colorBy, c("groups", "samples"))

  assertSubset(rtUnit, c("sec", "min"))

  ft <- components(object = obj,
                   samples = samples,
                   ID = ID,
                   mz = mz, ppm = ppm,
                   rt = rt, rtWindow = rtWindow, rtUnit = rtUnit,
                   comp = comp,
                   entireComponents = entireComponents,
                   onlyAnnotated = entireComponents,
                   onlyRelated = entireComponents)


  #plotting table
  ft$text <- paste0(round(ft$mz, digits = 2), ft$isotopes, ft$adduct)

  mzrange <- c(min(ft$mz) - 0.5, max(ft$mz) + 0.5)

  plot <- plot_ly(type = "bar")

  if (colorBy == "groups") {

    intft <- ft
    colnames(intft) <- str_replace_all(colnames(intft), "[[:punct:]]", "")
    spNames <- str_replace_all(obj@samples$sample[obj@samples$group %in% unique(ft$group)], "[[:punct:]]", "")
    names(spNames) <- obj@samples$group[obj@samples$group %in% unique(ft$group)]
    intft <- intft[, colnames(intft) %in% spNames, drop = FALSE]
    colnames(intft) <- names(spNames)
    intft$group <- ft$group

    ft$int <- apply(intft, MARGIN = 1, FUN = function(x) {
      x <- mean(as.numeric(x[names(x) %in% x[names(x) %in% "group"]]), na.rm = TRUE)
      ifelse(log, log(x), x)
    })

    intrange <- c(0, max(ft$int, na.rm = TRUE) * 2)

    rg <- unique(ft$group)

    colors <- getColors(length(rg))

    for (g in seq_len(length(rg))) {

      plot <- plot %>% add_bars(x = ft$mz[ft$group %in% rg[g]],
                                y = ft$int[ft$group %in% rg[g]],
                                marker = list(color = colors[g]),  width = 0.05,
                                text = ft$text[ft$group %in% rg[g]],
                                textposition = "outside", textangle = -90, insidetextanchor = "middle",
                                textfont = list(size = 10, color = colors[g]),
                                legendgroup = rg[g],
                                name = rg[g])
    }
  } else {

    colnames(ft) <- str_replace_all(colnames(ft), "[[:punct:]]", "")
    spNames <- str_replace_all(obj@samples$sample[obj@samples$group %in% unique(ft$group)], "[[:punct:]]", "")
    names(spNames) <- obj@samples$group[obj@samples$group %in% unique(ft$group)]

    intrange <- c(0, max(ft[, colnames(ft) %in% spNames], na.rm = TRUE) * 1.5)

    colors <- getColors(length(spNames))

    for (s in seq_len(length(spNames))) {
      plot <- plot %>% add_bars(x = ft$mz[ft$group %in% names(spNames)[s]],
                                y = ft[ft$group %in% names(spNames)[s], colnames(ft) %in% spNames[s], drop = TRUE],
                                marker = list(color = colors[s]),  width = 0.05,
                                text = ft$text[ft$group %in% names(spNames)[s]],
                                textposition = "outside", textangle = -90, insidetextanchor = "middle",
                                textfont = list(size = 10, color = colors[s]),
                                name = spNames[s])
    }
  }

  xaxis <- list(linecolor = toRGB("black"), linewidth = 2, title = "m/z",
                      titlefont = list(size = 12, color = "black"),
                      range = mzrange)

  yaxis <- list(linecolor = toRGB("black"), linewidth = 2, title = ifelse(log, "log(Intensity)", "Intensity"),
               titlefont = list(size = 12, color = "black"), range = intrange)

  plot <- plot %>% plotly::layout(xaxis = xaxis, yaxis = yaxis,
                                  barmode = "overlay", uniformtext = list(minsize = 4, mode = "show"))

  return(plot)

}
