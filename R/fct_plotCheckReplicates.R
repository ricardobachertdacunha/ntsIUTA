

#' @title plotCorReplicates
#' @description Evaluates the quality of replicate sample groups using a pearson
#' correlation of binned raw data.
#'
#' @param obj An \linkS4class{ntsData} object with more than one file.
#' @param binsize The bin size to apply to the raw data. The default is 2.
#'
#' @return A \emph{heatmap} plot to evaluate the correlation.
#'
#' @export
#'
#' @importFrom pheatmap pheatmap
#' @importClassesFrom MSnbase MSnExp
#' @importMethodsFrom MSnbase chromatogram bin
#'
plotCorReplicates <- function(obj, binsize = 2) {

  bpc <- MSnbase::chromatogram(obj@MSnExp, aggregationFun = "max")
  bpc <- MSnbase::bin(bpc, binSize = binsize)
  bpc <- log2(do.call(cbind, lapply(bpc, intensity)))
  bpc <- replace(bpc, bpc == "-Inf", 0)
  bpc <- cor(bpc)
  colnames(bpc) <- rownames(bpc) <- obj@MSnExp$sample_name

  ann <- data.frame(group = obj@MSnExp$sample_group) #for annotation of the heat map
  rownames(ann) <- obj@MSnExp$sample_name

  cl <- getColors(obj, which = "groups")

  clr <- colorRampPalette(brewer.pal(n = 7, name = "PRGn"))(100)

  pheatmap::pheatmap(bpc,
                     cluster_cols = FALSE,
                     cluster_rows = FALSE,
                     color = clr,
                     annotation = ann,
                     annotation_color = list(group = cl))

}




#' getCorReplicates
#'
#' @param obj An \linkS4class{ntsData} object with more than one file.
#' @param binsize The bin size to apply to the raw data. The default is 2.
#' @param exportResults Logical, set to \code{TRUE} for exporting the results to the project folder.
#'
#' @return A data frame with the minimum, average and respective sd of the correlation between
#' sample replicate groups.
#' 
#' @export
#' 
#' @importClassesFrom MSnbase MSnExp
#' @importMethodsFrom MSnbase chromatogram bin
#' @importFrom utils read.csv write.csv
#'
getCorReplicates <- function(obj, binsize = 2, exportResults = FALSE) {
  
  bpc <- MSnbase::chromatogram(obj@MSnExp, aggregationFun = "max")
  bpc <- MSnbase::bin(bpc, binSize = binsize)
  bpc <- log2(do.call(cbind, lapply(bpc, intensity)))
  bpc <- replace(bpc, bpc == "-Inf", 0)
  bpc <- cor(bpc)
  colnames(bpc) <- rownames(bpc) <- obj@MSnExp$sample_name
  
  df <- data.frame(replicate = unique(sampleGroups(obj)), cor_min = 0, cor_av = 0, cor_sd = 0)
  
  for (i in seq_len(nrow(df))) {
    sp <- samples(obj)[which(sampleGroups(obj) == df$replicate[i])]
    df$cor_min[i] <- min(bpc[rownames(bpc) %in% sp, colnames(bpc) %in% sp])
    df$cor_av[i] <- mean(bpc[rownames(bpc) %in% sp, colnames(bpc) %in% sp])
    df$cor_sd[i] <- sd(bpc[rownames(bpc) %in% sp, colnames(bpc) %in% sp])
  }
  
  if (exportResults) {
    
    results <- paste0(obj@path, "/results")
    if (!dir.exists(results)) dir.create(results)

    write.csv(df, file = paste0(results, "/RepliCor_df.csv"))
    
  }
  
  return(df)
  
}
