

#' catStep
#'
#' @param PL A features data frame.
#' @param IN The name of the sample replicate group as influent/input.
#' @param OUT The name of the sample replicate group as effluent/output.
#'
#' @return A data frame with features categorized as
#' removed (R), Lower (L), Constant (C), Higher (H) and New (N).
#' 
#' @importFrom dplyr mutate
#'
catStep <- function(PL, IN, OUT, constantLevel = 5000) {
    
  St <- PL
  St <- St[, 1:which(colnames(St) == "width")]
  St$intIN <- PL[, IN, drop = TRUE]
  St$intOUT <- PL[, OUT, drop = TRUE]
  St$qtl <-  St$intOUT - St$intIN
  St$qtlplot <- 0
  St <- St[!(St$intIN == 0 & St$intOUT == 0), ] #removes features not in replicate, thus 0
  
  St <- mutate(St, cat = ifelse(qtl < -constantLevel, "L", "wait"))
  St <- mutate(St, cat = ifelse(qtl > constantLevel, "H", cat))
  St$cat[St$intOUT <= 1000] <- "R"
  St$cat[St$intIN <= 1000] <- "N"
  St <- mutate(St, cat = ifelse(qtl >= -constantLevel & qtl <= constantLevel, "C", cat))

  St <- mutate(St, qtlplot = qtl/10^5)
  St <- mutate(St, qtlplot = ifelse(qtlplot > 3, 3.25, qtlplot))
  St <- mutate(St, qtlplot = ifelse(qtlplot < -3, -3.25, qtlplot))

  return(St)

}




### ntsMonitoringData -----

#' @title ntsMonitoringData
#'
#' @slot sequences A list of sequences reflecting the treatment/process chain.
#' The treatments/processes names in each sequence must match sample replicate group names
#' in an associated \linkS4class{ntsData} object.
#' @slot clusters A list of results from the clustering analysis for each sequence. 
#' @slot categories A list of results from the categorization of each step in each sequence.
#' @slot efficiency A list of data.frames with summarized categories results per sequence.
#'
#' @return An \linkS4class{ntsFragmentsData} object to be added to
#' the workflows slot of an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("ntsMonitoringData",
  slots = c(
    sequences = "list",
    parameters = "list",
    clusters_raw = "list",
    clusters = "list",
    categories = "list",
    efficiency = "list"
  ),
  prototype = list(
    sequences = list(),
    parameters = list(constantLevel = 5000,
                      minClust = 8,
                      sizeClus = 2:100),
    clusters_raw = list(),
    clusters = list(),
    categories = list(),
    efficiency = data.frame()
  )
)




#' @describeIn ntsMonitoringData Method to export plots.
#'
#' @param object An \linkS4class{ntsMonitoringData} object.
#' @param path A character vector with the path to export the monitoring results.
#' @param sequences Other parameters used in the production of process monitoring results.
#' @param yUnit The unit to show in the y axis.
#' Possible values are "mz" and "rt".
#' @param size The size for plotting the dots.
#' @param ID Optional character vector with features ID to filter the process data.
#' @param interactive Logical, set to \code{TRUE} to plot an interactive plot.
#' @param clusters The number of clusters to plot.
#' @param ncol The number of columns to arrange the clusters grid.
#'
#' @export
#'
setMethod("exportPlots", c("ntsMonitoringData", "character"), function(object, path = NULL,
                                                                       sequences = NULL,
                                                                       yUnit = "mz",
                                                                       size = NULL,
                                                                       ID = NULL,
                                                                       interactive = FALSE,
                                                                       clusters = NULL,
                                                                       ncol = 3) {
  
  
  
  if (missing(sequences)) sequences <- NULL
  if (missing(yUnit)) yUnit <- "mz"
  if (missing(ID)) ID <- NULL
  if (missing(interactive)) interactive <- FALSE
  if (missing(clusters)) clusters <- NULL
  if (missing(ncol)) ncol <- 3
  if (missing(size)) size <- NULL
  
  if (is.null(path)) path <- paste0(getwd(),"/results/monitoring")
  
  catplot <- plotProcessCategories(object, sequences = sequences, yUnit = yUnit,
                                   ID = ID, interactive = interactive, size = size)

  effplot <- plotProcessEfficiency(object, sequences = sequences)

  if (length(object@clusters) > 0) {
    clusplot <- plotProcessClusters(object, sequences = sequences,
                                    clusters = clusters, ID = ID, ncol = ncol,
                                    interactive = interactive)
  }
  
  
  for (i in seq_len(length(object@sequences))) {
    
    plot_path <- paste0(path, "/", names(object@sequences[i]))
    
    if (!dir.exists(plot_path)) dir.create(plot_path)
    
    if (interactive) {
      
      saveWidget(partial_bundle(catplot[[i]]),
                 file = paste0(plot_path,"/catPlot_", names(object@sequences[i]), ".html"))
      
      if (length(object@clusters) > 0) {
      
        k <- length(unique(object@clusters[[i]]$cluster))
        
        ggsave(paste0(plot_path,"/clusterPlot_", names(object@sequences[i]), ".pdf"),
               plot = clusplot[[i]],
               device = "pdf", path = NULL, scale = 1, width = 20,
               height = 6 * ceiling(k/ncol),
               units = "cm", dpi = 600, limitsize = TRUE)
      }
      
    } else {
    
      ggsave(paste0(plot_path,"/catPlot_", names(object@sequences[i]), ".tiff"),
             plot = catplot[[i]],
             device = "tiff", path = NULL, scale = 1, width = 20,
             height = 9 * ceiling(nrow(object@sequences[[i]])/2),
             units = "cm", dpi = 600, limitsize = TRUE)
      
      if (length(object@clusters) > 0) {
      
        k <- length(unique(object@clusters[[i]]$cluster))
        
        ggsave(paste0(plot_path,"/clusterPlot_", names(object@sequences[i]), ".pdf"),
               plot = clusplot[[i]],
               device = "pdf", path = NULL, scale = 1, width = 20,
               height = 6 * ceiling(k/ncol),
               units = "cm", dpi = 600, limitsize = TRUE)
      }
      
    }
    
    ggsave(paste0(plot_path,"/efficiencyPlot_", names(object@sequences[i]), ".tiff"),
           plot = effplot[[i]], device = "tiff",
           path = NULL, scale = 1, width = 17,
           height = 10, units = "cm", dpi = 300, limitsize = TRUE)
    
  }

})




#' processMonitoring
#'
#' @param obj An \linkS4class{ntsData} object with features for process monitoring.
#' @param sequences A list of sequences reflecting the treatment/process chain.
#' The treatments/processes names in each sequence must match sample replicate group names
#' in an associated \linkS4class{ntsData} object.
#' @param title Optional title for the workflow object entry in the \linkS4class{ntsData} object.
#' @param constantLevel The level, in counts, for assuming constant features through steps.
#' @param doClustering Logical, set to \code{TRUE} to perform clustering analysis for each sequence.
#' @param minClust The minimum number of clusters.
#' @param sizeClus The size range for clustering. The default is \code{2:100}.
#' @param exportResults Logical, set to \code{TRUE} to export
#' the standard monitoring results to the results forder.
#'
#' @return An \linkS4class{ntsMonitoringData} object added to the workflows slot
#' of the \linkS4class{ntsData} object.
#' 
#' @importFrom checkmate assertClass
#' @importFrom dplyr full_join inner_join rename select count
#' @importMethodsFrom  patRoon generateComponentsIntClust
#' @importFrom cluster silhouette
#' @importFrom stats cutree
#' @importFrom data.table setDT
#' @importFrom patRoon treeCut
#' @importFrom stringr str_extract
#' 
#' @export
#'
processMonitoring <- function(obj,
                              sequences = NULL,
                              title = NULL,
                              constantLevel = 5000,
                              doClustering = FALSE,
                              minClust = 10, sizeClus = 2:100, exportResults = FALSE) {
  
  assertClass(obj, "ntsData")
  
  ft <- obj@features
  
  if (nrow(ft) == 0) {
    warning("Featrues not found in the given ntsData.")
    return(obj)
  }
  
  if (is.null(sequences) | length(sequences) < 1) {
    warning("Treatment sequences not found.")
    return(obj)
  }
  
  #check if all sequences are larger than 2
  check <- rep(FALSE, length(sequences))
  check <- sapply(sequences, function(x) length(x) > 2)
  
  if (!all(check)) {
    warning("Sequences with less than 2 sample replicate groups not possible but found.")
    return(obj)
  }
  
  data <- new("ntsMonitoringData")
  
  monitorSeq <- list()
  
  for (s in seq_len(length(sequences))) {
    if (length(sequences[[s]]) > 2) {
    monitorSeq[[names(sequences)[s]]] <- data.frame(step = c(), IN = c(), OUT = c())
      for (step in seq_len((length(sequences[[s]]) - 1))) {
        monitorSeq[[names(sequences)[s]]] <- rbind(monitorSeq[[names(sequences)[s]]],
                                                   data.frame(step = step,
                                                   IN = sequences[[s]][step],
                                                   OUT = sequences[[s]][step + 1]))
      }
    }
  }
  
  data@sequences <- monitorSeq
  
  if (is.null(constantLevel)) constantLevel <- data@parameters$constantLevel
  
  monitorCat <- list()
  
  monitorEff <- list()
  
  for (s in seq_len(length(sequences))) {
    monitorCat[[names(sequences)[s]]] <- list()
    for (step in seq_len(nrow(monitorSeq[[s]]))) {
      monitorCat[[names(sequences)[s]]][[monitorSeq[[s]]$step[step]]] <- catStep(PL = ft,
        IN = monitorSeq[[s]]$IN[step],
        OUT = monitorSeq[[s]]$OUT[step],
        constantLevel = constantLevel)
      
      m1 <- monitorCat[[names(sequences)[s]]][[monitorSeq[[s]]$step[step]]]
      m1 <- m1 %>% dplyr::count(cat)
      colnames(m1) <- c("cat", paste0("Step", step, "_", monitorSeq[[s]]$IN[step], "_", monitorSeq[[s]]$OUT[step]))
      
      if (step == 1) { m2 <- m1 } else {
        m2 <- full_join(m2, m1, by = c("cat"))
      }
    }
    
    m2$cat <- factor(m2$cat, levels = c("N","H","C","L","R"),
                     labels = c("New", "Higher","Constant", "Lower","Removed"))
    
    monitorEff[[names(sequences)[s]]] <- m2
  }
  
  data@categories <- monitorCat
  
  data@efficiency <- monitorEff
  
  if (doClustering) {
    
    if (is.null(minClust)) minClust <- data@parameters$minClust
    if (is.null(sizeClus)) sizeClus <- data@parameters$sizeClus
  
    monitorClustRaw <- list()
    
    monitorClust <- list()
    
    for (s in seq_len(length(sequences))) {
      if (length(sequences[[s]]) > 2) {
        
        monitorClustRaw[[names(sequences)[s]]] <- generateComponentsIntClust(
          obj@patdata[which(sampleGroups(obj) %in% sequences[[s]]), ft$ID],
          method = "complete")
        
        meanws <- sapply(sizeClus, function(k) {
          sil <- cluster::silhouette(cutree(monitorClustRaw[[names(sequences)[s]]]@clust, k = k), 
                                     monitorClustRaw[[names(sequences)[s]]]@distm)
          summary(sil)$avg.width
        })
      
        k <- 1 + minClust + match(max(meanws[-c(1:minClust)]), meanws[-c(1:minClust)])
      
        clus0 <- patRoon::treeCut(monitorClustRaw[[names(sequences)[s]]], k = k)
        
        #make summary for components
        clus1 <- clus0@components
        clus1 <- rbindlist(clus1, idcol = "cluster")
        clus1 <- rename(clus1, ID = group)
        clus2 <- as.data.frame(clus0@clusterm)
        clus2 <- data.table::setDT(clus2, keep.rownames = TRUE)
        clus2 <- rename(clus2, ID = rn)
        clus2 <- inner_join(clus1, clus2, by = "ID")
        clus2 <- rename(clus2, rt = ret)
        clus2 <- select(clus2, -intensity)
        clus2$cluster <- as.numeric(str_extract(clus2$cluster, "[:digit:]"))
        
        monitorClust[[names(sequences)[s]]] <- clus2
        
      }
    }
    
    data@parameters$minClust <- minClust
    
    data@parameters$sizeClus <- sizeClus
    
    data@clusters_raw <- monitorClustRaw
    
    data@clusters <- monitorClust
    
  }
  
  if (exportResults) {
    
    exportPath <- paste0(obj@path,"/results")
    if (!dir.exists(exportPath)) dir.create(exportPath)
    
    exportPath <- paste0(exportPath,"/monitoring")
    if (!dir.exists(exportPath)) dir.create(exportPath)
    
    exportPlots(data, path = exportPath)
    
  }
  
  if (is.null(title)) title <- "processMonitoring"

  obj@workflows[[title]] <- data
  
  return(obj)
  
}
