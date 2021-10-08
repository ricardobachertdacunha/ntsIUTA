

#' catPlot
#'
#' @param Step A data frame with the categorized monitor data for a given step.
#' @param StepN The step number.
#' @param StepName The step name.
#' @param yUnit The unit to show in the y axis.
#' Possible values are "mz" and "rt".
#' @param size The size to plot the dots.
#' @param constantLevel The instensity variation threshold for constant features between steps.
#' @param ID Optional character vector with features ID to filter the process data.
#' 
#' @import ggplot2
#'
#' @return A ggplot with the features shown by categories according to the quantiles.
#'
catPlot <- function(Step, StepN, StepName, yUnit = "mz", size = NULL, constantLevel = 5000, ID = NULL) {
  
  if (!is.null(ID)) Step <- Step[Step$ID %in% ID, ]
  
  if (is.null(size)) size <- 0.5
  
  if (nrow(Step) == 0) return(list())
  
  if (yUnit == "rt") Step$mz <- Step$rt
    
  plot <- ggplot(data = Step) +
    geom_point(mapping = aes(x = qtlplot, y = mz, colour = cat), size = size, alpha = 1/1.2) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_line(colour = "#C0C0C0", size = 0.1),
                   legend.position =  c(0.73, 1.09),
                   legend.background = element_blank(),
                   legend.direction = "horizontal") +
    labs(title = StepName,
         subtitle = paste0("Constant level threshold set at ", constantLevel, " counts"),
         colour = "Category") +
    xlab(expression("Intensity ( In - Out )  x10"^"5")) +
    ylab(ifelse(yUnit == "mz", expression(italic("m/z")), expression("Retention time (sec)"))) +
    scale_color_manual(values = c("R" = "#006400", "L" = "#90EE90", "C" = "#87CEFA", "H" = "#FFA500", "N" = "#8B0000")) +
    scale_x_continuous(breaks = seq(-3, 3, 1)) +
    expand_limits(x = c(-3.5, 3.5)) +
    geom_vline(xintercept = -3, size = 1) +   
    geom_vline(xintercept = -0.05, size = 0.5) +
    geom_vline(xintercept = 0.05, size = 0.5) +
    geom_vline(xintercept = 3, size = 1)
  
  return(plot)
  
}




#' catPlotInteractive
#'
#' @param Step A data frame with the categorized monitor data for a given step.
#' @param StepN The step number.
#' @param StepName The step name.
#' @param yUnit The unit to show in the y axis.
#' Possible values are "mz" and "rt".
#' @param size The size to plot the dots.
#' @param showL A vector with TRUE and FALSE to hide legends in subplots.
#' @param constantLevel The instensity variation threshold for constant features between steps.
#' @param ID Optional character vector with features ID to filter the process data.
#' 
#' @importFrom plotly plot_ly add_trace toRGB layout
#'
#' @return A plotly with the features shown by categories according to the quantiles.
#'
catPlotInteractive <- function(Step, StepN, StepName, yUnit = "mz", size = NULL,
                               showL = TRUE, constantLevel = 5000, ID = NULL) {
  
  if (!is.null(ID)) Step <- Step[Step$ID %in% ID, ]
  
  if (is.null(size)) size <- 3
  
  if (nrow(Step) == 0) return(list())
  
  plot <- plot_ly()
  
  if (yUnit == "mz") {
    Step$yunit <- Step$mz
  } else {
    Step$yunit <- Step$rt
  }
  
  plot <- plot %>% add_trace(Step,
          x = Step$qtlplot,
          y = Step$yunit,
          type = "scatter",
          mode = 'markers',
          marker = list(size = size),
          color = Step$cat,
          colors = c("R" = "#006400", "L" = "#90EE90", "C" = "#87CEFA", "H" = "#FFA500", "N" = "#8B0000"),
          legendgroup = Step$cat, 
          showlegend = showL,
          hoverinfo = "text",
          text = paste("</br> feature: ", Step$ID,
                       "</br> <i>m/z</i>: ", round(Step$mz, digits = 4),
                       "</br> rt: ", round(Step$rt, digits = 0),
                       "</br> intIN: ", round(Step$intIN, digits = 0),
                       "</br> intOUT: ", round(Step$intOUT, digits = 0)))
  
  
  vline1 <- list(x0 = -3, x1 = -3, y0 = min(Step$yunit, na.rm = TRUE), y1 = max(Step$yunit, na.rm = TRUE), line = list(width = 0.8))
  vline2 <- list(x0 = -0.05, x1 = -0.05, y0 = min(Step$yunit, na.rm = TRUE), y1 = max(Step$yunit, na.rm = TRUE), line = list(width = 0.5))
  vline3 <- list(x0 = 0.05, x1 = 0.05, y0 = min(Step$yunit, na.rm = TRUE), y1 = max(Step$yunit, na.rm = TRUE), line = list(width = 0.5))
  vline4 <- list(x0 = 3, x1 = 3, y0 = min(Step$yunit, na.rm = TRUE), y1 = max(Step$yunit, na.rm = TRUE), line = list(width = 0.8))

  title <- list(text = paste0("Constant level threshold set at ", constantLevel, " counts"),
                x = 0.1, y = 1.07, font = list(size = 10, color = "black"))
  
  subtitle <- list(x = 0.02, y = 1.05,
                  text = StepName,
                  font = list(size = 13, color = "black"),
                  showarrow = FALSE, xref = 'paper', yref = 'paper')
  
  xaxis <- list(linecolor = toRGB("black"),
               linewidth = 2, title = "Intensity ( In - Out )  x 10E5",
               titlefont = list(size = 12, color = "black"),
               range = c(-3.5,3.5),
               autotick = TRUE, ticks = "outside")
  
  yaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = ifelse(yUnit == "mz", paste("<i>m/z</i>"), "Retetnion time (sec.)"),
                titlefont = list(size = 12, color = "black"))
  
  plot <- plot %>% plotly::layout(xaxis = xaxis,
                                  yaxis = yaxis,
                                  title = title,
                                  annotations = list(subtitle),
                                  shapes = list(vline1, vline2, vline3, vline4))
  
  plot
  
  return(plot)
  
}




#' plotProcessCategories
#'
#' @param obj An \linkS4class{ntsMonitoringData} object with monitoring categories to plot.
#' @param sequences Optional vector with sequence names from the \linkS4class{ntsMonitoringData} object to plot.
#' When \code{NULL}, all sequences are plotted.
#' @param ID Optional character vector with features ID to filter the process data. 
#' @param yUnit The unit to use for the y axis. Possible values are "mz" or "rt.
#' @param size The size for ploting the dots.
#' @param interactive Logical, set to \code{TRUE} to plot an interactive plot.
#'
#' @return A list with plots for each given sequence in the \linkS4class{ntsMonitoringData} object.
#' 
#' @importFrom checkmate assertClass
#' @importFrom plotly subplot
#' @importFrom gridExtra grid.arrange
#' 
#' @export
#'
plotProcessCategories <- function(obj, sequences = NULL, ID = NULL,
                                  yUnit = "mz", size = NULL, interactive = TRUE) {
  
  assertClass(obj, "ntsMonitoringData")
  
  obj2 <- obj
  
  if (!is.null(sequences)) {
    if (all(sequences %in% names(obj2@sequences))) {
      obj2@sequences <- obj@sequences[sequences]
      obj2@categories <- obj@categories[sequences]
    } else {
      warning("One or more sequences not found!")
      return(list())
    }
  }
  
  monitorPlot <- list()
  
  for (i in seq_len(length(obj2@sequences))) {
    
    monitorPlot[[names(obj2@sequences[i])]] <- list()
    
    showLegend <- rep(FALSE, nrow(obj2@sequences[[i]]))
    showLegend[1] <- TRUE 
    
    for (j in seq_len(nrow(obj2@sequences[[i]]))) {
      
      
      if (interactive) {
        
        monitorPlot[[names(obj2@sequences[i])]][[obj2@sequences[[i]]$step[j]]] <- catPlotInteractive(
          Step =  obj2@categories[[names(obj2@sequences[i])]][[obj2@sequences[[i]]$step[j]]],
          StepN = j,
          StepName = paste0("Step ", j, ": ", obj2@sequences[[i]]$IN[j], " to ", obj2@sequences[[i]]$OUT[j]),
          yUnit = yUnit,
          showL = showLegend[j],
          constantLevel = obj2@parameters$constantLevel,
          ID = ID,
          size = size)
        
      } else {
        
        monitorPlot[[names(obj2@sequences[i])]][[obj2@sequences[[i]]$step[j]]] <- catPlot(
          Step =  obj2@categories[[names(obj2@sequences[i])]][[obj2@sequences[[i]]$step[j]]],
          StepN = j,
          StepName = paste0("Step ", j, ": ", obj2@sequences[[i]]$IN[j], " to ", obj2@sequences[[i]]$OUT[j]),
          yUnit = yUnit,
          constantLevel = obj2@parameters$constantLevel,
          ID = ID,
          size = size)
      
      }
      
      
    }
    
    if (interactive) {
      
      monitorPlot[[names(obj2@sequences[i])]] <- plotly::subplot(monitorPlot[[names(obj2@sequences[i])]],
                                                                 nrows = length(monitorPlot[[names(obj2@sequences[i])]]),
                                                                 margin = 0.04, shareX = TRUE, shareY = TRUE)
      
    } else {
      
      monitorPlot[[names(obj2@sequences[i])]] <- gridExtra::grid.arrange(grobs = monitorPlot[[names(obj2@sequences[i])]], ncol = 1)
      
    }
    
  }
  
  return(monitorPlot)
  
}




#' plotProcessEfficiency
#'
#' @param obj An \linkS4class{ntsMonitoringData} object with categorized features
#' to evaluate the process efficiency.
#' @param sequences Optional vector with sequence names from the \linkS4class{ntsMonitoringData} object to plot.
#' When \code{NULL}, all sequences are plotted.
#'
#' @return A bar plot with the number of features
#' for each category for each step in the respective process sequence.
#' 
#' @importFrom reshape2 melt
#' @import ggplot2
#' 
#' @export
#'
plotProcessEfficiency <- function(obj, sequences = NULL) {
  
  assertClass(obj, "ntsMonitoringData")
  
  obj2 <- obj
  
  if (!is.null(sequences)) {
    if (all(sequences %in% names(obj2@sequences))) {
      obj2@sequences <- obj@sequences[sequences]
      obj2@efficiency <- obj@efficiency[sequences]
    } else {
      warning("One or more sequences not found!")
      return(list())
    }
  }
  
  efficiencyPlot <- list()

  for (i in seq_len(length(obj2@efficiency))) {
    
    tempEff <- reshape2::melt(obj2@efficiency[[i]], id.vars = "cat")
    
    efficiencyPlot[[names(obj2@efficiency[i])]] <- ggplot(data = tempEff) +
          theme_bw() +
          geom_bar(mapping = aes(x = variable, y = value, fill = cat), stat = "identity", width = 0.4, alpha = 1/1.3) +
          theme(axis.title.x = element_blank()) +
          scale_fill_manual(name = "Category", drop = FALSE,
            values = c("New" = "#8B0000", "Higher" = "#FFA500", "Constant" = "#87CEFA", "Lower" = "#90EE90", "Removed" = "#006400"))+
          ylab("Number of Features") +
          labs(title = paste0("Sequence: ", names(obj2@efficiency[i])),
               subtitle = paste0("Constant level threshold set at ", obj2@parameters$constantLevel, " counts"))
  
  }
  
  return(efficiencyPlot)
  
}



#' plotProcessClusters
#'
#' @param obj An \linkS4class{ntsMonitoringData} object with cluster analysis for process evaluation.
#' @param sequences Optional vector with sequence names from the \linkS4class{ntsMonitoringData} object to plot.
#' When \code{NULL}, all sequences are plotted.
#' @param clusters The number of clusters to plot.
#' @param ID The ID of features to plot.
#' @param interactive Logical, set to \code{TRUE} for using plotly instead of ggplot2.
#' @param ncol The number of columns to arrange the clusters grid.
#'
#' @return A grid plot list with the intensity trends of features in each cluster for each sequence.
#' 
#' @export
#'
plotProcessClusters <- function(obj, sequences = NULL,
                                clusters = NULL, ID = NULL,
                                interactive = FALSE, ncol = 3) {
  
  assertClass(obj, "ntsMonitoringData")
  
  obj2 <- obj
  
  if (!is.null(sequences)) {
    if (all(sequences %in% names(obj2@sequences))) {
      obj2@sequences <- obj@sequences[sequences]
      obj2@clusters <- obj@clusters[sequences]
    } else {
      warning("One or more sequences not found!")
      return(list())
    }
  }
  
  clusterPlot <- list()

  for (i in seq_len(length(obj2@clusters))) {
    
    cl <- obj2@clusters[[i]]
    
    clusterPlot[[names(obj2@clusters[i])]] <- list()
    
    if (!is.null(ID)) cl <- cl[cl$ID %in% ID, ]
    
    if (!is.null(clusters)) cl <- cl[cl$cluster %in% clusters, ]
    
    for (cN in seq_len(max(cl$cluster))) {
    
      cls <- cl[cl$cluster == cN, ]
      
      if (nrow(cls) > 0) {
      
      cls <- cls %>% tidyr::pivot_longer(!c("cluster", "ID", "mz", "rt"), names_to = "step", values_to = "int")
      
      cls$step <- factor(cls$step, levels = unique(c(obj2@sequences[[i]]$IN, obj2@sequences[[i]]$OUT)),
                         labels = unique(c(obj2@sequences[[i]]$IN, obj2@sequences[[i]]$OUT)))
      
      plot <- ggplot(cls, aes(x = step, y = int, group = ID)) +
                geom_line(aes(colour = mz), size = 0.2) +
                #scale_colour_gradient2() +
                theme_bw() +
                scale_x_discrete(expand = c(0.01, 0.01), drop=FALSE) +
                theme(panel.grid.major.x = element_blank(),
                      panel.grid.minor.x = element_blank(),
                      panel.grid.minor.y = element_blank(),
                      axis.title.x = element_blank(),
                      plot.title = element_text(size = 10),
                      axis.text.y = element_text(size = 5),
                      axis.title.y = element_text(size = 8),
                      legend.position = "none") +
                labs(title = paste0("Cluster ", unique(cls$cluster), " (", nrow(cls), " feature/s)")) +
                ylab("Intensity (relative)")
      
      clusterPlot[[names(obj2@clusters[i])]][[cN]] <- plot
      
      }
      
    }
    
    clusterPlot[[names(obj2@clusters[i])]] <- gridExtra::grid.arrange(grobs = clusterPlot[[names(obj2@clusters[i])]],
                                                                      ncol = ncol)
    
  }
  
  return(clusterPlot)
  
}
