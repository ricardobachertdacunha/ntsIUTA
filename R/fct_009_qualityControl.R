

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
#' @param recursive Logical, set to \code{TRUE} for applying recursive
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
#' @importFrom checkmate assertClass checkChoice
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
#' @examples
#'
checkQC <- function(
  obj = NULL,
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


  assertClass(obj, "ntsData")

  if (nrow(obj@QC$samples) == 0) {
    warning("No sample replicate group/s assinged for QC in the ntsData object!")
    return(obj)
  }

  if (is.null(algorithmPeakPicking)) algorithmPeakPicking <- obj@algorithms$peakPicking

  if (!testChoice(algorithmPeakPicking, c("xcms3", "xcms", "openms", "envipick"))) {
    warning("Algorithm not recognized. See ?peakPicking for more information.")
    return(obj)
  }

  if (is.null(algorithmMakeFeatures)) algorithmMakeFeatures <- obj@algorithms$makeFeatures

  if (!testChoice(algorithmMakeFeatures, c("xcms3", "xcms", "openms"))) {
    warning("Algorithm not recognized. See ?makeFeatures for more information.")
    return(obj)
  }

  if (is.null(paramPeakPicking)) paramPeakPicking <- obj@parameters$peakPicking
  if (is.null(paramAlignment)) paramAlignment <- obj@parameters$peakAlignment
  if (is.null(paramGrouping)) paramGrouping <- obj@parameters$peakGrouping
  if (is.null(paramFill)) paramFill <- obj@parameters$fillMissing
  # TODO Add check for parameters

  data <- new("ntsData")

  data@samples <- obj@QC$samples

  data@polarity <- obj@polarity

  data <- importRawData(data,
                        rtFilter = NULL,
                        rtUnit = "min",
                        centroidedData = NA,
                        removeEmptySpectra = TRUE,
                        save = FALSE)

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

  if (!is.data.frame(targets)) targets <- utils::read.csv(targets)

  if (max(targets$rt) < 120) targets$rt <- targets$rt * 60

# TODO change check from supects to targets
  if (!("adduct" %in% colnames(targets))) {
    adduct <- ifelse(data@polarity == "positive", "[M+H]+",
                     ifelse(data@polarity == "negative", "[M-H]-",
                            stop("polarity argument must be 'positive' or 'negative'")))
  } else {
    adduct <- NULL
  }

  screen <- screenSuspects(data@patdata, select(targets, -mz),
                           rtWindow = rtWindow, mzWindow = 0.02,
                           adduct = adduct, onlyHits = TRUE)

  df <- arrange(patRoon::as.data.frame(screen, average = FALSE), group)
  df <- left_join(df, targets[, colnames(targets) %in% c("name", "formula", "hasFragments", "intControl")], by = "name")
  df <- left_join(df, select(arrange(screenInfo(screen), group), group, d_mz, d_rt), by = "group")
  df$av_into <- rowMeans(select(df, data@samples$sample))
  df <- df %>% mutate(sd_into = apply(select(., data@samples$sample), 1, sd))
  df$sd_into <- round(df$sd_into, digits = 1)
  df <- mutate(df, sd_intop = round(sd_into / av_into * 100, digits = 1))
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- dplyr::rename(df, rt = ret, ID = group)
  df <- select(df, name, formula, d_ppm, d_rt, ID, mz, rt,
               av_into, sd_into, sd_intop, everything(), -d_mz)
  df <- dplyr::filter(df, d_ppm <= ppm)
  df$d_ppm <- round(df$d_ppm, digits = 1)
  df$d_rt <- round(df$d_rt, digits = 1)

  if ("hasFragments" %in% colnames(df)) {

    ms2df <- left_join(df[, c("name", "ID", "hasFragments")],
                       targets[, colnames(targets) %in% c("name", "fragments_mz", "fragments_int", "fragments_pre")],
                       by = "name")

    ms2df <- mutate(ms2df, hasFragments_Exp = FALSE,
                    fragments_mz_Exp = NA, fragments_int_Exp = NA, fragments_pre_Exp = NA)

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

    for (i in 1:nrow(ms2df)) {

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

  obj@QC$targets <- targets

  obj@QC$data <- data

  obj@QC$results <- df

  obj@algorithms$peakPicking <- algorithmPeakPicking

  obj@algorithms$makeFeatures <- algorithmMakeFeatures

  obj@parameters$peakPicking <- paramPeakPicking

  obj@parameters$peakAlignment <- paramAlignment

  obj@parameters$peakGrouping <- paramGrouping

  obj@parameters$fillMissing <- paramFill

  #TODO Implement isotope check for QC targets? By using annotation results.

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
#' @examples
#'
plotCheckQC <- function(obj, rtWindow = 30, ppm = 15) {

  qcdf <- obj@QC$results

  if (nrow(qcdf) == 0) return(cat("QC results not found or empty."))

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
