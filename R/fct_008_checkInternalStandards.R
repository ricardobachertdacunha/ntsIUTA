

#' @title checkIS
#' @description Evaluation of internal standards for samples within a given \linkS4class{XCMSnExp} or \linkS4class{featureGroups} object.
#'
#' @param x A \linkS4class{XCMSnExp} object.
#' @param xPat Alternatively, a \linkS4class{featureGroups} object.
#' @param polarity The polarity for the search. Possible values are \code{positive} or \code{negative}.
#' @param screeningList A \code{data.frame} with details from each internal standard in the samples.
#' See details for more information about the required data.frame structure.
#' @param ppmWindow The mass deviation, in ppm, allowed to screen for internal standards.
#' @param rtWindow The retention time deviation, in seconds, allowed to screen for internal standards.
#' @param plot Logical, set to \code{TRUE} for plotting the results.
#' @param save Logical, set to \code{TRUE} for storing  the generated IS object in the rData folder.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#'
#' @return A \code{list} containing the plots and a data.frame with the summary of the
#' retention time, mass and intensity deviations as well as the quality of the MS2 data for each internal standard found.
#' 
#' @details The \code{screeningList} template can be obtained as csv via \code{ntsIUTA::getScreeningListTemplate()}.
#' Add other details of the template.
#' 
#' @export
#' 
#' @import magrittr
#' @importFrom utils read.csv write.csv
#' @importFrom patRoon screenSuspects screenInfo getDefAvgPListParams generateMSPeakLists as.data.table importFeatureGroupsXCMS3 analysisInfo analyses
#' @importFrom dplyr select arrange left_join everything filter mutate top_n
#' @importFrom fuzzyjoin difference_inner_join
#' @importFrom stats cor sd
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly partial_bundle
#' @importMethodsFrom MSnbase fileNames
#'
#' @examples
#' 
#' 
#' 
checkIS <- function(x = NULL,
                    xPat = NULL,
                    polarity = "positive",
                    screeningList = utils::choose.files(base::getwd(), "Select the QC screening list"),
                    ppmWindow = 15,
                    rtWindow = 30,
                    plot = TRUE,
                    save = TRUE, projPath = setup$projPath) {
  
  # TODO Adapt to ntsData framework
  
  # x = featDataSample
  # xPat = patDataRef
  # polarity = "positive"
  # screeningList = utils::choose.files(base::getwd(), "Select the QC screening list")
  # ppmWindow = 15
  # rtWindow = 30
  
  if (base::is.null(x) & base::is.null(xPat)) stop("group data should be provided using either x or xPat arguments. See ?checkIS for more information.")
  
  #getting patData if not given
  if (base::is.null(xPat)) {
    patSampleInfo <- base::data.frame(path = base::dirname(MSnbase::fileNames(x)),
                                      analysis = x$sample_name,
                                      group = x$sample_group,
                                      blank = "")
    xPat <- patRoon::importFeatureGroupsXCMS3(x, patSampleInfo)
  }
  
  #load or read screeningList
  if (!base::is.data.frame(screeningList)) {
    sl <- utils::read.csv(screeningList)
  } else {
    sl <- screeningList
  }
  if (base::max(sl$rt) < 120) sl$rt <- sl$rt * 60
  
  #filterIS
  #sl <- sl[sl$comment == "IS",]
  
  #check polarity
  if (base::length(polarity) == 1 & !base::is.null(polarity)) {
    adduct <- base::ifelse(polarity == "positive", "[M+H]+", base::ifelse(polarity == "negative", "[M-H]-", stop("polarity argument must be 'positive' or 'negative'")))
  } else {
    stop("polarity argument must be 'positive' or 'negative'")
  }
  
  isID <- patRoon::screenSuspects(xPat, sl, rtWindow = rtWindow, mzWindow = 0.01, adduct = adduct, onlyHits = TRUE)
  isdf <- dplyr::arrange(patRoon::as.data.table(isID, average = TRUE), group)
  isdf <- dplyr::left_join(isdf, sl[,base::c("name", "int10", "hasMS2", "mzMS2", "intMS2", "preMS2")], by = "name")
  isdf  <- dplyr::left_join(isdf, dplyr::select(dplyr::arrange(patRoon::screenInfo(isID), group), group, d_mz, d_rt), by = "group")
  isdf$av_into <- base::rowMeans(dplyr::select(isdf, base::unique(patRoon::analysisInfo(xPat)$group)))
  isdf <- isdf %>% dplyr::mutate(sd_into = base::apply(dplyr::select(., base::unique(patRoon::analysisInfo(xPat)$group)), 1, sd))
  isdf <- dplyr::mutate(isdf, sd_intop = sd_into/av_into*100)
  isdf$d_ppm <- (base::abs(isdf$d_mz)/isdf$mz)*1E6
  isdf <- dplyr::select(isdf, name, formula, d_ppm, d_rt, mz, ret, group, dplyr::everything(), -d_mz)
  isdf<- dplyr::filter(isdf, d_ppm <= ppmWindow)
  isdf <- base::as.data.frame(isdf)
  
  control_avgPListParams <- patRoon::getDefAvgPListParams(
    clusterMzWindow = 0.005,
    topMost = 50,
    minIntensityPre = 10,
    minIntensityPost = 10
  )
  
  MS2 <- base::suppressWarnings(patRoon::generateMSPeakLists(
    isID, "mzr",
    maxMSRtWindow = 5,
    precursorMzWindow = 3,
    avgFeatParams = control_avgPListParams, 
    avgFGroupParams = control_avgPListParams
  ))
  
  #add MS2 to screening list, adds intensity at 10ng/ml to respective column and adds mz corresponding to the polarity for MS2 matching
  # for (i in 1:base::nrow(isdf)) {
  #   xgroup <- isdf$group[i]
  #   xfrag <- MS2[[xgroup]]$MSMS
  # 
  #   if (!base::is.null(xfrag)) {
  #     sl$hasMS2[sl$name %in% isdf$name[i]] <- TRUE
  #     sl$mzMS2[sl$name %in% isdf$name[i]] <- base::paste(xfrag$mz, collapse = ";")
  #     sl$intMS2[sl$name %in% isdf$name[i]] <- base::paste(xfrag$intensity, collapse = ";")
  #     sl$preMS2[sl$name %in% isdf$name[i]] <- base::paste(xfrag$precursor, collapse = ";")
  #   } else {
  #      sl$hasMS2[sl$name %in% isdf$name[i]] <- FALSE
  #   }
  # 
  #   sl$int10[sl$name %in% isdf$name[i]] <- isdf$av_into[i]
  #   sl$rt[sl$name %in% isdf$name[i]] <- isdf$ret[i]/60
  #   # sl$mz[sl$name %in% isdf$name[i]] <- sl$neutralMass[sl$name %in% isdf$name[i]] +  base::ifelse(polarity == "positive", 1.007276, -1.007276)
  # }
  # utils::write.csv(sl, file = base::paste0(setup$projPath,"/IS_ScreeningList_ntsIUTA_MS2_pos.csv"))
  
  # ntsIUTA::plotRawChrom(MSnbase::filterFile(rawData, 7:11) , type = "tic", mz = 751.5045, ppm = 15)
  # 
  # ntsIUTA::plotFeaturePeaks(patDataRef, fileIndex = NULL, 
  #                           features = c("M752_R850_5641","M752_R866_5642"),
  #                           mz = NULL, ppm = 15,
  #                           rt = NULL, rtWindow = 1,
  #                           rtUnit = "min",
  #                           plotBy = "samples")
  
  
  isdf$nfrag <- 0
  isdf$pfrag <- 0
  isdf$intoscore <- 0
  
  for (i in 1:base::nrow(isdf)) {
    xgroup <- isdf$group[i]
    xname <- isdf$name[i]
    if (!base::is.na(xgroup) & isdf$hasMS2[i]) {
      xMS2 <- MS2[[xgroup]]$MSMS
      if (!base::is.null(xMS2)) {
        dbMS2 <- base::data.frame(mz = base::as.numeric(base::unlist(base::strsplit(isdf$mzMS2[i], split=";"))),
                                  intensity = base::as.numeric(base::unlist(base::strsplit(isdf$intMS2[i], split=";"))),
                                  precursor = base::as.logical(base::unlist(base::strsplit(isdf$preMS2[i], split=";"))))
        
        xMS2 <- dplyr::top_n(xMS2, 10, intensity)
        xMS2 <- dplyr::mutate(xMS2, into_ind = intensity/base::max(xMS2$intensity))
        
        dbMS2 <- dplyr::top_n(dbMS2, 10, intensity)
        dbMS2 <- dplyr::mutate(dbMS2, into_ind = intensity/base::max(dbMS2$intensity))
        
        combi <- fuzzyjoin::difference_inner_join(xMS2, dbMS2, by = c("mz"), max_dist = 0.005, distance_col = "diff")
        
        isdf$nfrag[i] <- base::nrow(combi)
        isdf$pfrag[i] <- base::nrow(combi)/base::nrow(dbMS2)
        
        if (isdf$nfrag[i] == 1) { isdf$intoscore[i] <- 1 }
        else {
          isdf$intoscore[i] <- stats::cor(combi$into_ind.x, combi$into_ind.y, use = "everything", method = "pearson")
        }
      } else {
        isdf$nfrag[i] <- 0
        isdf$pfrag[i] <- 0
        isdf$intoscore[i] <- 0
      }
    }
  }
  
  isdf$pfrag <- base::round(isdf$pfrag, digits = 2)
  isdf$intoscore <- base::round(base::as.numeric(isdf$intoscore), digits = 4)
  
  IS <- base::list()
  IS[["df"]] <- dplyr::select(isdf, -hasMS2, -mzMS2, -intMS2, -preMS2)

  # if(plot) {
  # 
  #   evalPlot <- base::list()
  # 
  #   evalPlot[[1]] <- ggplot2::ggplot(isdf) +
  #     ggplot2::theme_bw() +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = -10, ymax = 10, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = -15, ymax = -10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = 15, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = -15, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = 15, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
  #     ggplot2::geom_point(stat = "identity",
  #                         ggplot2::aes(x = name,
  #                                      y = d_rt)) +
  #     ggplot2::theme(axis.title.y = ggplot2::element_blank(),
  #                    axis.text.y = ggplot2::element_text(size = 7),
  #                    axis.text.x = ggplot2::element_text(size = 7),
  #                    axis.title.x = ggplot2::element_text(size = 7),
  #                    plot.title = ggplot2::element_text(hjust = 0.5, size = 9)) +
  #     ggplot2::labs(title = paste("")) +
  #     ggplot2::ylab("RT diff (sec)") +
  #     ggplot2::ylim(-rtWindow,rtWindow) +
  #     ggplot2::coord_flip()
  # 
  #   evalPlot[[2]] <- ggplot2::ggplot(isdf) +
  #     ggplot2::theme_bw() +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 5, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = 5, ymax = 10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
  #     ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
  #     ggplot2::geom_point(stat = "identity",
  #                         ggplot2::aes(x = name,
  #                                      y = d_ppm)) +
  #     ggplot2::theme(axis.text.y = ggplot2::element_blank(),
  #                    axis.title.y = ggplot2::element_blank(),
  #                    axis.text.x = ggplot2::element_text(size = 7),
  #                    axis.title.x = ggplot2::element_text(size = 7),
  #                    plot.title = ggplot2::element_text(hjust = 0.5, size = 9)) +
  #     ggplot2::labs(title = paste("")) +
  #     ggplot2::ylab("m/z diff (ppm)") +
  #     ggplot2::ylim(0,ppmWindow) +
  #     ggplot2::coord_flip()
  # 
  #   samplesInt <- base::as.data.frame(patRoon::as.data.table(isID[,isdf$group], average = F))
  #   rGroups <- base::unique(patRoon::analysisInfo(xPat)$group)
  # 
  #   for (i in 1:base::length(rGroups)) {
  #     sNames <- patRoon::analyses(isID)[patRoon::analysisInfo(xPat)$group == rGroups[i]]
  #     tempdf <- base::data.frame(A = isdf$name,
  #                                B = base::apply(samplesInt[,base::colnames(samplesInt) %in% sNames], 1, function(x) base::mean(x, na.rm = T)),
  #                                C = base::apply(samplesInt[,base::colnames(samplesInt) %in% sNames], 1, function(x) stats::sd(x, na.rm = T)),
  #                                D = isdf$int10)
  #     evalPlot[[2 + i]] <- ggplot2::ggplot(data = tempdf) +
  #       ggplot2::theme_bw() +
  #       ggplot2::geom_rect(ggplot2::aes(ymin = 50, ymax = 150, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
  #       ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = 50, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
  #       ggplot2::geom_rect(ggplot2::aes(ymin = 150, ymax = 190, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
  #       ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 10, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
  #       ggplot2::geom_rect(ggplot2::aes(ymin = 190, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
  #       ggplot2::geom_errorbar(ggplot2::aes(x = A,
  #                                           ymin = base::ifelse((B/D*100) == 0, 0, base::ifelse(B/D*100 > 200, 200, B/D*100) - (base::ifelse(B/D*100 > 200, 200, B/D*100)*C/B)),
  #                                           ymax = base::ifelse((B/D*100) == 0, 0, base::ifelse(B/D*100 > 200, 200, B/D*100) + (base::ifelse(B/D*100 > 200, 200, B/D*100)*C/B))),
  #                                           width = 0.2, position = ggplot2::position_dodge(0.9),
  #                                           colour = "gray45") +
  #       ggplot2::geom_point(stat = "identity", ggplot2::aes(x = A, y = base::ifelse(B/D*100 > 200, 200, B/D*100)), size = 2) +
  #       ggplot2::theme(axis.title.y = ggplot2::element_blank(),
  #             axis.text.y = ggplot2::element_blank(),
  #             axis.text.x = ggplot2::element_text(size = 7),
  #             axis.title.x = ggplot2::element_text(size = 7),
  #             plot.title = ggplot2::element_text(hjust = 0.5, size = 9),
  #             panel.grid.minor.x = ggplot2::element_blank()) +
  #       ggplot2::labs(title = paste(rGroups[i])) +
  #       ggplot2::ylab("Recovery (%)") +
  #       ggplot2::scale_y_continuous(limits = c(0, 210), breaks = base::seq(0, 200, by = 50)) +
  #       ggplot2::coord_flip()
  #   }
  #   
  #   widths <- c(8,5, base::rep(5, base::length(rGroups)))
  #   
  #   evalPlot <- gridExtra::grid.arrange(grobs = evalPlot, nrow = 1, heights = 10, widths = widths)
  #   
  #   featPlot <- ntsIUTA::plotFeaturePeaks(xPat, fileIndex = NULL,
  #                                         features = isdf$group,
  #                                         mz = NULL, ppm = NULL,
  #                                         rt = NULL, rtWindow = NULL,
  #                                         rtUnit = "min",
  #                                         plotBy = "features",
  #                                         names = isdf$name)
  #   
  #   if (save) {
  #     results <- base::paste0(projPath,"\\results")
  #     if (!base::dir.exists(results)) base::dir.create(results)
  #     ggplot2::ggsave(base::paste0(projPath,"/results/IS_Deviations.tiff"),
  #                     plot = evalPlot, device = "tiff", path = NULL, scale = 1,
  #                     width = base::sum(widths), height = 10, units = "cm", dpi = 300, limitsize = TRUE)
  #     htmlwidgets::saveWidget(plotly::partial_bundle(featPlot), file = base::paste0(projPath,"/results/IS_Features.html"))
  #   }
  #   
  #   IS[["featPlot"]] <- featPlot
  #   IS[["evalPlot"]] <- evalPlot
  #   
  # }

  if (save) {
    results <- base::paste0(projPath,"\\results")
    if (!base::dir.exists(results)) base::dir.create(results)
    utils::write.csv(dplyr::select(isdf, -hasMS2, -mzMS2, -intMS2, -preMS2), file = base::paste0(projPath,"/results/IS_CheckResults.csv"))
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(IS, file = base::paste0(rData,"\\IS.rds"))
  }
  
  return(IS)
  
}
