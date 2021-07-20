

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
setMethod("plotCorReplicates", "ntsData", function(obj, binsize = 2) {
  
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
  
})
