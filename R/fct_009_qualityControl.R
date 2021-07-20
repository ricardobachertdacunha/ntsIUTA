

#' @title checkQC
#' @description Function to check QC samples in an \linkS4class{ntsData} object.
#' 
#' @param obj An \linkS4class{ntsData} object, containing QC samples.
#' @param targets A \code{data.frame} or a location for a csv
#' with details for each standard in the QC samples.
#' See details for more information about the required data.frame structure.
#' @param algorithmPeakPicking The algorithm to use for peak picking.
#' One of "xcms3" (the default), "xcms", "openms" or "envipick".
#' @param algorithmMakeFeatures The algorithm to use for peak alignment and grouping.
#' One of "xcms3" (the default), "xcms" or "openms".
#' @param paramPeakPicking  List of arguments for the specified algorithm.
#' See \code{\link[patRoon]{findFeatures}} for more information.
#' For instance, for "xcms3" as \code{algorithm}, \code{param} should be given
#' with the method parameters for the peak finding.
#' See \href{https://rdrr.io/bioc/xcms/man/chromatographic-peak-detection.html}{\code{xcms::findChromPeaks}}
#' for more information.
#' @param rtAlignment Logical, set to \code{TRUE} (The default) for performing peak alignment before grouping.
#' @param paramAlignment Applicable for algorithm "xcms3" only,
#' the parameters for the chosen retention time alignment method.
#' See documentation of \code{\link[xcms]{adjustRtime}} for more information.
#' For alignment with the method \code{PeakGroups}, a list of length two
#' should be given, where the first element
#' is the pre-grouping parameters using \code{\link[xcms]{groupChromPeaks}}
#' and the second element the actual alignment parameters.
#' @param paramGrouping he parameters for the chosen grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} or
#' \code{\link[patRoon]{groupFeatures}} for more information.
#' @param recurvive Logical, set to \code{TRUE} for applying recursive
#' integration for samples with missing peaks.
#' @param paramFill An object of class \code{FillChromPeaksParam}
#' or \code{ChromPeakAreaParam} containing the parameters
#' to apply the recursive integration.
#' #' See \code{?\link[xcms]{fillChromPeaks}} for more information.
#' @param rtWindow The retention time deviation, in seconds,
#' to screen for the QC target standards. The default is 30 seconds.
#' @param ppm The mass deviation, in ppm, to screen for the QC target standards.
#' The default is 15 ppm.
#' @param exportResults Logical, set to \code{TRUE} for exporting the results to the project folder.
#' @param save  Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @details The \code{screeningList} template can be obtained as csv via \code{ntsIUTA::getScreeningListTemplate()}.
#' Add other details of the template.
#'
#' @return The \linkS4class{ntsData} with quality control data and
#' evaluation results added to the slot \code{QC}.
#' 
#' @export
#'
#' @importFrom checkmate checkChoice
#' @importFrom utils read.csv write.csv
#' @importClassesFrom patRoon featureGroupsScreening
#' @importMethodsFrom patRoon screenSuspects screenInfo generateMSPeakLists as.data.table
#' @importFrom patRoon getDefAvgPListParams
#' @importFrom dplyr select arrange left_join everything filter mutate top_n
#' @importFrom fuzzyjoin difference_inner_join
#' @importFrom stats cor
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly partial_bundle
#' @importFrom xcms applyAdjustedRtime
#'
#'
#' @examples
#' 
#' 
#' 
checkQC <- function(obj = NULL,
                    targets = utils::choose.files(base::getwd(), "Select the QC screening list"),
                    algorithmPeakPicking = NULL,
                    algorithmMakeFeatures = NULL,
                    paramPeakPicking = NULL,
                    rtAlignment = TRUE,
                    paramAlignment = NULL,
                    paramGrouping = NULL,
                    recursive = TRUE,
                    paramFill = NULL,
                    rtWindow = 30,
                    ppm = 15,
                    exportResults = FALSE,
                    save = FALSE) {
  
  
  if (is.null(obj)) return(cat("An ntsData object with QC sample replicate group/s assigned must be given!"))
  
  if (nrow(obj@QC$samples) == 0) return(cat("No sample replicate group/s assinged for QC in the ntsData object!"))
  
  if (is.null(algorithmPeakPicking)) algorithmPeakPicking <- obj@algorithms$peakPicking
  
  checkmate::checkChoice(algorithmPeakPicking, c("xcms3", "xcms", "openms", "envipick"))
  
  if (is.null(algorithmMakeFeatures)) algorithmMakeFeatures <- obj@algorithms$makeFeatures
  
  checkmate::checkChoice(algorithmMakeFeatures, c("xcms3", "xcms", "openms"))
  
  if (is.null(paramPeakPicking)) paramPeakPicking <- obj@parameters$peakPicking
  if (is.null(paramAlignment)) paramAlignment <- obj@parameters$peakAlignment
  if (is.null(paramGrouping)) paramGrouping <- obj@parameters$peakGrouping
  if (is.null(paramFill)) paramFill <- obj@parameters$fillMissing
  
  data <- new("ntsData")
  
  data@samples <- obj@QC$samples
  
  data@polarity <- obj@polarity
  
  data <- peakPicking(data, algorithm = algorithmPeakPicking,
                      param = paramPeakPicking, save = FALSE)
  
  data <- makeFeatures(data,
                       algorithm = algorithmMakeFeatures,
                       rtAlignment = rtAlignment,
                       paramAlignment = paramAlignment,
                       paramGrouping = paramGrouping,
                       recursive = recursive,
                       paramFill = paramFill,
                       save = FALSE)
  
  
  targets <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),"/QC_ScreeningList_ntsIUTA_MS2_pos.csv")
  
  if (!is.data.frame(targets)) targets <- utils::read.csv(targets)
  
  if (max(targets$rt) < 120) targets$rt <- targets$rt * 60
  
  adduct <- ifelse(data@polarity == "positive", "[M+H]+",
                   ifelse(data@polarity == "negative", "[M-H]-",
                          stop("polarity argument must be 'positive' or 'negative'")))
  
  
  screen <- screenSuspects(data@patdata, select(targets, -mz),
                           rtWindow = rtWindow, mzWindow = 0.02,
                           adduct = adduct, onlyHits = TRUE)
  
  df <- arrange(patRoon::as.data.frame(screen, average = FALSE), group)
  df <- left_join(df, targets[,c("name", "formula", "int10", "hasMS2", "mzMS2", "intMS2", "preMS2")], by = "name")
  df <- left_join(df, select(arrange(screenInfo(screen), group), group, d_mz, d_rt), by = "group")
  df$av_into <- rowMeans(select(df, data@samples$sample))
  df <- df %>% mutate(sd_into = apply(select(., data@samples$sample), 1, sd))
  df <- mutate(df, sd_intop = sd_into / av_into * 100)
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- select(df, name, formula, d_ppm, d_rt, mz, ret, group, everything(), -d_mz)
  df <- filter(df, d_ppm <= ppm)
  
  
  control_avgPListParams <- getDefAvgPListParams(
    clusterMzWindow = 0.003,
    topMost = 50,
    minIntensityPre = 10,
    minIntensityPost = 10
  )
  
  MS2 <- suppressWarnings(generateMSPeakLists(
    screen, "mzr",
    maxMSRtWindow = 10,
    precursorMzWindow = 1.5,
    avgFeatParams = control_avgPListParams, 
    avgFGroupParams = control_avgPListParams
  ))
  
  df$nfrag <- 0
  df$pfrag <- 0
  df$intoscore <- 0
  
  for (i in 1:nrow(df)) {
    xgroup <- df$group[i]
    xname <- df$name[i]
    if (!is.na(xgroup) & df$hasMS2[i]) {
      xMS2 <- MS2[[xgroup]]$MSMS
      if (!is.null(xMS2)) {
        dbMS2 <- data.frame(mz = as.numeric(unlist(strsplit(df$mzMS2[i], split = ";"))),
                            intensity = as.numeric(unlist(strsplit(df$intMS2[i], split = ";"))),
                            precursor = as.logical(unlist(strsplit(df$preMS2[i], split = ";"))))
        
        xMS2 <- top_n(xMS2, 10, intensity)
        xMS2 <- mutate(xMS2, into_ind = intensity / max(xMS2$intensity))
        
        dbMS2 <- top_n(dbMS2, 10, intensity)
        dbMS2 <- mutate(dbMS2, into_ind = intensity / max(dbMS2$intensity))
        
        combi <- fuzzyjoin::difference_inner_join(xMS2, dbMS2,
                                                  by = c("mz"),
                                                  max_dist = 0.005,
                                                  distance_col = "diff")
        
        df$nfrag[i] <- nrow(combi)
        df$pfrag[i] <- nrow(combi) / nrow(dbMS2)
        
        if (df$nfrag[i] == 1) { df$intoscore[i] <- 1 }
        else {
          df$intoscore[i] <- stats::cor(combi$into_ind.x,
                                        combi$into_ind.y,
                                        use = "everything",
                                        method = "pearson")
        }
      } else {
        df$nfrag[i] <- 0
        df$pfrag[i] <- 0
        df$intoscore[i] <- 0
      }
    }
  }
  
  df$pfrag <- round(df$pfrag, digits = 2)
  
  df$intoscore <- round(as.numeric(df$intoscore), digits = 4)
  
  obj@QC$targets <- targets
  
  obj@QC$data <- data
  
  obj@algorithms$peakPicking <- algorithmPeakPicking
  
  obj@algorithms$makeFeatures <- algorithmMakeFeatures
  
  obj@parameters$peakPicking <- paramPeakPicking
  
  obj@parameters$peakAlignment <- paramAlignment
  
  obj@parameters$peakGrouping <- paramGrouping
  
  obj@parameters$fillMissing <- paramFilling
  
  
  
  #TODO Implement isotope check for QC targets?
  
  QC <- list()
  QC[["df"]] <- dplyr::select(df, -hasMS2, -mzMS2, -intMS2, -preMS2)
  
  
  
  if (plot) {
    evalPlot <- gridExtra::arrangeGrob(
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = -10, ymax = 10, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -15, ymax = -10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = 15, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = -15, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 15, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = d_rt)) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(size = 7),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("RT diff (sec)") +
        ggplot2::ylim(-rtWindow,rtWindow) +
        ggplot2::coord_flip(),
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 5, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 5, ymax = 10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = d_ppm)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("m/z diff (ppm)") +
        ggplot2::ylim(0,ppmWindow) +
        ggplot2::coord_flip(), 
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = 2, ymax = 5, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 5, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 2, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = nfrag)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("Nr Shared Top10 MS2") +
        #ylim(0,10) +
        ggplot2::scale_y_continuous(limits = c(0, 10), breaks = base::seq(0, 10, by = 2)) +
        ggplot2::coord_flip(),
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = 0.95, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 0.9, ymax = 0.95, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 0.9, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = intoscore)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("MS2 Intensity Corr.") +
        ggplot2::ylim(0.7,1) +
        ggplot2::coord_flip(),
      ncol = 4, widths = base::c(8,4,4,4))
    
    
    featPlot <- ntsIUTA::plotFeaturePeaks(qcPat, fileIndex = NULL,
                                          features = qcdf$group,
                                          mz = NULL, ppm = NULL,
                                          rt = NULL, rtWindow = NULL,
                                          rtUnit = "min",
                                          plotBy = "features",
                                          names = qcdf$name)
    
    if (save) {
      results <- base::paste0(projPath,"\\results")
      if (!base::dir.exists(results)) base::dir.create(results)
      ggplot2::ggsave(base::paste0(projPath,"/results/QC_Deviations.tiff"),
                      plot = evalPlot, device = "tiff", path = NULL, scale = 1,
                      width = 17, height = 10, units = "cm", dpi = 300, limitsize = TRUE)
      htmlwidgets::saveWidget(plotly::partial_bundle(alignPlot), file = base::paste0(projPath,"/results/QC_Alignment.html"))
      htmlwidgets::saveWidget(plotly::partial_bundle(featPlot), file = base::paste0(projPath,"/results/QC_Features.html"))
    }
    
    QC[["alignPlot"]] <- alignPlot
    QC[["featPlot"]] <- featPlot
    QC[["evalPlot"]] <- evalPlot
    
  }
  
  if (save) {
    results <- base::paste0(projPath,"\\results")
    if (!base::dir.exists(results)) base::dir.create(results)
    utils::write.csv(dplyr::select(qcdf, -hasMS2, -mzMS2, -intMS2, -preMS2), file = base::paste0(projPath,"/results/QC_CheckResults.csv"))
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(QC, file = base::paste0(rData,"\\QC.rds"))
  }
  
  return(QC)
  
}

#add MS2 to screening list, adds intensity at 10ng/ml to respective column and adds mz corresponding to the polarity for MS2 matching
# for (i in 1:base::nrow(qcdf)) {
#   xgroup <- qcdf$group[i]
#   xfrag <- MS2[[xgroup]]$MSMS
# 
#   if (!base::is.null(xfrag)) {
#     sl$hasMS2[sl$name %in% qcdf$name[i]] <- TRUE
#     sl$mzMS2[sl$name %in% qcdf$name[i]] <- base::paste(xfrag$mz, collapse = ";")
#     sl$intMS2[sl$name %in% qcdf$name[i]] <- base::paste(xfrag$intensity, collapse = ";")
#     sl$preMS2[sl$name %in% qcdf$name[i]] <- base::paste(xfrag$precursor, collapse = ";")
#   }
# 
#   sl$int10[sl$name %in% qcdf$name[i]] <- qcdf$av_into[i]
#   sl$rt[sl$name %in% qcdf$name[i]] <- qcdf$ret[i]/60
#   sl$mz[sl$name %in% qcdf$name[i]] <- sl$neutralMass[sl$name %in% qcdf$name[i]] +  base::ifelse(polarity == "positive", 1.007276, -1.007276)
# }
# 
# utils::write.csv(sl, file = base::paste0(projPath,"/ScreeningList_QC_ntsIUTA_MS2_pos.csv"))
