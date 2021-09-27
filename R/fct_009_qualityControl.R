

#' @title checkQC
#' @description Function to check QC samples in an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object, containing QC samples.
#' @param targets The \linkS4class{suspectList} object with the target compounds
#' for quality control.
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
#' @importFrom checkmate assertClass checkChoice
#' @importClassesFrom xcms XCMSnExp
#' @importMethodsFrom xcms chromPeaks
#' @importFrom utils read.csv write.csv
#' @importClassesFrom patRoon featureGroupsScreening
#' @importMethodsFrom patRoon screenSuspects screenInfo generateMSPeakLists as.data.table
#' @importFrom patRoon getDefAvgPListParams
#' @importFrom dplyr select arrange left_join everything mutate filter top_n rename
#' @importFrom fuzzyjoin difference_inner_join
#' @importFrom stats cor
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly partial_bundle
#'
#'
checkQC <- function(
  obj = NULL,
  targets = NULL,
  rtWindow = NULL,
  ppm = NULL,
  exportResults = FALSE,
  save = FALSE) {

  assertClass(obj, "ntsData")

  assertClass(targets, "suspectList")

  if (targets@length == 0) {
    warning("The targets list is empty!")
    return(obj)
  }

  if (is.null(rtWindow)) rtWindow <- obj@QC@rtWindow
  
  if (is.null(ppm)) ppm <- obj@QC@ppm
    
  #creates a temporary ntsData to produce results
  data <- new("ntsData")

  data@samples <- obj@QC@samples

  data@samples$blank <- ""

  data@polarity <- obj@polarity
  
  data@parameters <- obj@parameters

  data <- importRawData(data,
                        rtFilter = NULL,
                        rtUnit = "min",
                        centroidedData = NA,
                        removeEmptySpectra = TRUE,
                        save = FALSE)

  data <- peakPicking(data, save = FALSE)

  data <- makeFeatures(data, save = FALSE)

  data <- annotateFeatures(data, excludeBlanks = FALSE, save = FALSE)

  
  if (!("adduct" %in% colnames(targets@data))) {
    adduct <- ifelse(data@polarity == "positive", "[M+H]+",
                     ifelse(data@polarity == "negative", "[M-H]-",
                            stop("polarity argument must be 'positive' or 'negative'")))
  } else {
    adduct <- NULL
  }

  require(xcms)

  screen <- screenSuspects(data@patdata, select(targets@data, -mz),
                           rtWindow = rtWindow, mzWindow = 0.02,
                           adduct = adduct, onlyHits = TRUE)

  df <- arrange(patRoon::as.data.frame(screen, average = FALSE), group)
  df <- rename(df, name = susp_name)
  df <- left_join(df, targets@data[, colnames(targets@data) %in% c("name", "formula", "hasFragments", "intControl")], by = "name")
  df <- left_join(df, select(arrange(screenInfo(screen), group), group, d_mz, d_rt), by = "group")
  df$av_into <- rowMeans(select(df, data@samples$sample))
  df <- df %>% mutate(sd_into = apply(select(., data@samples$sample), 1, sd))
  df$sd_into <- round(df$sd_into, digits = 1)
  df <- mutate(df, sd_intop = round(sd_into / av_into * 100, digits = 1))
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- dplyr::rename(df, rt = ret, ID = group)
  df$isoN <- 0
  df <- select(df, name, formula, d_ppm, d_rt, ID, mz, rt,
               av_into, sd_into, sd_intop, isoN, everything(), -d_mz)

  df <- dplyr::filter(df, d_ppm <= ppm)
  df$d_ppm <- round(df$d_ppm, digits = 1)
  df$d_rt <- round(df$d_rt, digits = 1)

  screen <- screen[, df$ID]

  annot <- data@annotation$comp[data@annotation$comp$ID %in% df$ID, ]
  for (i in seq_len(nrow(annot))) {
    df$isoN[i] <- nrow(data@annotation$comp[data@annotation$comp$isogroup %in% annot$isogroup[i], ]) - 1
  }

  if ("hasFragments" %in% colnames(df)) {

    ms2df <- left_join(df[, c("name", "ID", "hasFragments")],
                       targets@data[, colnames(targets@data) %in% c("name", "fragments_mz", "fragments_int", "fragments_pre")],
                       by = "name")

    ms2df <- mutate(ms2df,
                    hasFragments_Exp = FALSE,
                    fragments_mz_Exp = NA,
                    fragments_int_Exp = NA,
                    fragments_pre_Exp = NA)

    MS2 <- extractMS2(screen, param = data@parameters@MS2)

    df$nfrag <- 0

    df$pfrag <- 0

    df$intoscore <- 0

    for (i in seq_len(nrow(ms2df))) {

      xgroup <- ms2df$ID[i]

      xname <- ms2df$name[i]

      if (!is.na(xgroup)) {

        xMS2 <- MS2[[xgroup]]$MSMS

        if (!is.null(xMS2)) {

          ms2df$hasFragments_Exp[i] <- TRUE
          ms2df$fragments_mz_Exp[i] <- paste(xMS2$mz, collapse = ";")
          ms2df$fragments_int_Exp <- paste(xMS2$intensity, collapse = ";")
          ms2df$fragments_pre_Exp <- paste(xMS2$precursor, collapse = ";")

          if (ms2df$hasFragments[i]) {

            dbMS2 <- data.frame(mz = as.numeric(unlist(strsplit(ms2df$fragments_mz[i], split = ";"))),
                                intensity = as.numeric(unlist(strsplit(ms2df$fragments_int[i], split = ";"))),
                                precursor = as.logical(unlist(strsplit(ms2df$fragments_pre[i], split = ";"))))

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

            if (df$nfrag[i] == 1) {
              df$intoscore[i] <- 1
            } else {
              df$intoscore[i] <- stats::cor(combi$into_ind.x,
                                            combi$into_ind.y,
                                            use = "everything",
                                            method = "pearson")
            }
          }
        }
      }
    }

    df$pfrag <- round(df$pfrag, digits = 2)

    df$intoscore <- round(as.numeric(df$intoscore), digits = 4)

    df <- left_join(df, select(ms2df, -name, -hasFragments), by = "ID")

  }

  obj@QC@targets <- targets

  obj@QC@rtWindow <- rtWindow
  
  obj@QC@ppm <- ppm
  
  obj@QC@MSnExp <- data@MSnExp

  obj@QC@patdata <- data@patdata
  
  obj@QC@peaks <- data@peaks

  obj@QC@features <- data@features
  
  obj@QC@annotation <- data@annotation
  
  obj@QC@results <- df

  if (exportResults) {

    results <- paste0(obj@path, "/results")
    if (!dir.exists(results)) dir.create(results)

    write.csv(df, file = paste0(results, "/QC01_df.csv"))

    plotqc <- plotCheckQC(obj, rtWindow = rtWindow, ppm = ppm)

    ggsave(paste0(results, "/QC02_deviations.tiff"),
           plot = plotqc, device = "tiff", path = NULL, scale = 1,
           width = 17, height = 10, units = "cm", dpi = 300, limitsize = TRUE)

    plotfp <- plotFeaturePeaks(obj = data,
                               ID = obj@QC$results$ID,
                               mz = NULL,
                               rt = NULL,
                               rtUnit = "min",
                               ppm = NULL,
                               rtWindow = NULL,
                               names = obj@QC$results$name,
                               interactive = TRUE)

    saveWidget(partial_bundle(plotfp), file = paste0(results, "/QC03_featurePeaks.html"))

  }

  return(obj)

}


#' @title plotCheckQC
#'
#' @param obj obj An \linkS4class{ntsData} object, containing QC samples.
#' @param rtWindow The time window to show the time deviation.
#' @param ppm The ppm window to show the mass deviation.
#'
#' @return A plot of the QC check results.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#'
plotCheckQC <- function(obj, rtWindow = NULL, ppm = NULL) {
  
  qcdf <- obj@QC@results

  if (nrow(qcdf) == 0) return(cat("QC results not found or empty."))

  if (is.null(rtWindow)) rtWindow <- obj@QC@rtWindow
  
  if (is.null(ppm)) ppm <- obj@QC@ppm
  
  grobs <- list()

  grobs[["rt"]] <- ggplot(qcdf) +
    theme_bw() +
    geom_rect(aes(ymin = -10, ymax = 10, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
    geom_rect(aes(ymin = -15, ymax = -10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
    geom_rect(aes(ymin = 10, ymax = 15, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
    geom_rect(aes(ymin = -Inf, ymax = -15, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
    geom_rect(aes(ymin = 15, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
    geom_point(stat = "identity", aes(x = name, y = d_rt)) +
    theme(axis.title.y = element_blank(),
          axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 7),
          axis.title.x = element_text(size = 7)) +
    ylab("RT diff (sec)") +
    ylim(-rtWindow, rtWindow) +
    coord_flip()

  grobs[["ppm"]] <- ggplot(qcdf) +
    theme_bw() +
    geom_rect(aes(ymin = -Inf, ymax = 5, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
    geom_rect(aes(ymin = 5, ymax = 10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
    geom_rect(aes(ymin = 10, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
    geom_point(stat = "identity", aes(x = name, y = d_ppm)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.title.x = element_text(size = 7)) +
    ylab("m/z diff (ppm)") +
    ylim(0, ppm) +
    coord_flip()

  if ("intControl" %in% colnames(qcdf)) {

    qcdf$intRecovery <- (qcdf$av_into / qcdf$intControl) * 100

    grobs[["int"]] <- ggplot(data = qcdf) +
      theme_bw() +
      geom_rect(aes(ymin = 50, ymax = 150, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 10, ymax = 50, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 150, ymax = 190, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
      geom_rect(aes(ymin = -Inf, ymax = 10, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_rect(aes(ymin = 190, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_point(stat = "identity", aes(x = name, y = ifelse(intRecovery > 200, 200, intRecovery))) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 7),
            axis.title.x = element_text(size = 7),
            panel.grid.minor.x = element_blank()) +
      ylab("Recovery (%)") +
      scale_y_continuous(limits = c(0, 210), breaks = seq(0, 200, by = 50)) +
      coord_flip()
  }

  if ("isoN" %in% colnames(qcdf)) {

    grobs[["isoN"]] <- ggplot(qcdf) +
      theme_bw() +
      geom_rect(aes(ymin = -Inf, ymax = 2, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.05) +
      geom_rect(aes(ymin = 2, ymax = 5, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 5, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_point(stat = "identity", aes(x = name, y = isoN)) +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),
            axis.title.x = element_text(size = 7)) +
      ylab("Nr Isotopes") +
      ylim(0, 6) +
      coord_flip()

  }

  if ("hasFragments" %in% colnames(qcdf)) {

    grobs[["frag"]] <- ggplot(qcdf) +
      theme_bw() +
      geom_rect(aes(ymin = 2, ymax = 5, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 5, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
      geom_rect(aes(ymin = -Inf, ymax = 2, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_point(stat = "identity", aes(x = name, y = nfrag)) +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),
            axis.title.x = element_text(size = 7)) +
      ylab("Nr Shared Top10 MS2") +
      scale_y_continuous(limits = c(0, 10), breaks = base::seq(0, 10, by = 2)) +
      coord_flip()

    grobs[["intoscore"]] <- ggplot(qcdf) +
      theme_bw() +
      geom_rect(aes(ymin = 0.95, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 0.9, ymax = 0.95, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
      geom_rect(aes(ymin = -Inf, ymax = 0.9, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_point(stat = "identity", aes(x = name, y = intoscore)) +
      theme(axis.text.y = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 7),
            axis.title.x = element_text(size = 7)) +
      ylab("MS2 Intensity Corr.") +
      ylim(0.7, 1) +
      coord_flip()

  }

  return(grid.arrange(grobs = grobs,
                      ncol = length(grobs),
                      widths = c(8, rep(4, length(grobs) - 1))))

}




#' qc2ntsData
#' 
#' @description Converts the QC data in an \linkS4class{ntsData} object to a new
#' \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object containing QC results.
#'
#' @return An \linkS4class{ntsData} object with QC data added to the main data slots.
#' 
#' @export
#'
qc2ntsData <- function(obj) {
  
  assertClass(obj, "ntsData")
  
  if (nrow(obj@QC@features) == 0) {
    warning("QC data not found in the given ntsData!")
    return(NULL)
  }
  
  newObj <- new("ntsData")
  
  newObj@title <- "qcData"
  
  newObj@path <- obj@path
  
  newObj@date <- obj@date
  
  newObj@polarity <- obj@polarity
  
  newObj@samples <- obj@QC@samples
  
  newObj@parameters <- obj@parameters
  
  newObj@QC@rtWindow <- obj@QC@rtWindow
  
  newObj@QC@ppm <- obj@QC@ppm
  
  newObj@MSnExp <- obj@QC@MSnExp
  
  newObj@patdata <- obj@QC@patdata
  
  newObj@peaks <- obj@QC@peaks
  
  newObj@features <- obj@QC@features
  
  newObj@annotation <- obj@QC@annotation
  
  newObj@QC@results <- obj@QC@results
  
  return(newObj)
  
}
