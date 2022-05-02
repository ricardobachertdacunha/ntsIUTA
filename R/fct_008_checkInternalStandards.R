

#' @title checkIS
#' @description Evaluation of internal standards for each sample replicate group
#' within the \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param targets The \linkS4class{suspectList} object with the information of
#' internal standards to be found.
#' @param ppm The mass deviation, in ppm, to screen for the internal standards.
#' The default is 10 ppm.
#' @param rtWindow The retention time deviation, in seconds,
#' to screen for the internal standards. The default is 30 seconds.
#' @param MS2param The parameters to extract MS2 information of each internal standard.
#' @param recoveryFrom A character vector with the name of the sample replicate group
#' to use for calculating the recovery.
#' The default is \code{NULL} to use the intControl in the \linkS4class{suspectList}.
#' @param exportResults Logical, set to \code{TRUE} for exporting the results
#' to the project folder.
#' @param save  Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object with results from evaluation of
#' internal standards in each sample replicate group added to the IS slot.
#'
#' @export
#'
#' @importFrom checkmate assertClass checkChoice
#' @importFrom utils read.csv write.csv
#' @importClassesFrom patRoon featureGroupsScreening
#' @importMethodsFrom patRoon screenSuspects screenInfo generateMSPeakLists as.data.table
#' @importFrom patRoon getDefAvgPListParams
#' @importFrom dplyr select arrange left_join everything mutate filter top_n rename all_of
#' @importFrom fuzzyjoin difference_inner_join
#' @importFrom stats cor
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly partial_bundle
#'
checkIS <- function(obj = NULL,
                    targets = NULL,
                    ppm = NULL,
                    rtWindow = NULL,
                    MS2param = NULL,
                    recoveryFrom = NULL,
                    exportResults = FALSE,
                    save = FALSE) {

  assertClass(obj, "ntsData")

  if (is.null(targets)) targets <- obj@IS@targets
  
  if (nrow(targets@data) == 0) {
    warning("Targets not found in the given ntsData. They should be provided!")
    return(obj)
  }
  
  if (is.null(ppm)) ppm <- obj@IS@ppm
  
  if (is.null(rtWindow)) rtWindow <- obj@IS@rtWindow
  
  if (is.null(MS2param)) MS2param <- obj@parameters@MS2
  
  if (is.null(recoveryFrom)) recoveryFrom <- obj@IS@recoveryFrom
  
  if (!("adduct" %in% colnames(targets@data))) {
    adduct <- ifelse(obj@polarity == "positive", "[M+H]+",
                     ifelse(obj@polarity == "negative", "[M-H]-",
                            stop("polarity argument must be 'positive' or 'negative'")))
  } else {
    adduct <- NULL
  }

  screen <- screenSuspects(obj@patdata, select(targets@data, -mz),
                           rtWindow = rtWindow, mzWindow = 0.02,
                           adduct = adduct, onlyHits = FALSE)

  df <- arrange(patRoon::as.data.table(screen, average = FALSE, onlyHits = TRUE), group)
  df <- base::as.data.frame(df)
  df <- rename(df, name = susp_name)
  df <- left_join(df, select(targets@data, -mz, -rt), by = "name")
  df <- left_join(df, select(arrange(screenInfo(screen), group), group, d_mz, d_rt), by = "group")
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- rename(df, rt = ret, ID = group)
  df <- select(df, name, formula, d_ppm, d_rt, ID, mz, rt, everything(), -d_mz)

  ## TODO Idea to add infor from rt and mz variance from peaks in the set

  # pksrt <- obj@peaks[obj@peaks$feature %in% df$ID, ]
  # pksrt <- pksrt %>% select(feature, rt) %>% group_by(feature) %>% summarize(sd(rt))
  # colnames(pksrt) <- c("ID", "rt_sd")
  # df <- left_join(df, pksrt, by = c("ID"))
  # 
  # pksmz <- obj@peaks[obj@peaks$feature %in% df$ID, ]
  # pksmz <- pksmz %>% select(feature, mz) %>% group_by(feature) %>% summarize(sd(rt))
  # colnames(pksrt) <- c("ID", "rt_sd")
  # df <- left_join(df, pksrt, by = c("ID")

  df <- dplyr::filter(df, d_ppm <= ppm)
  df$d_ppm <- round(df$d_ppm, digits = 1)
  df$d_rt <- round(df$d_rt, digits = 1)

  screen <- screen[, df$ID]

  #evaluate per replicate group
  ls <- list()

  rg <- unique(sampleGroups(obj))

  for (g in seq_len(length(rg))) {

    sp <- samples(obj)[sampleGroups(obj) %in% rg[g]]

    inv_sp <- samples(obj)[!(sampleGroups(obj) %in% rg[g])]

    temp <- df[, !(colnames(df) %in% inv_sp)]

    temp$av_into <- rowMeans(select(temp, all_of(sp)))

    temp <- temp %>% mutate(sd_into = apply(select(., all_of(sp)), 1, sd))

    temp$sd_into <- round(temp$sd_into, digits = 1)

    temp <- mutate(temp, sd_intop = round(sd_into / av_into * 100, digits = 1))

    if (is.null(recoveryFrom) | is.na(recoveryFrom)) {

      temp$recovery <- temp$av_into / temp$intControl
      temp$recovery_sd <- (temp$sd_into / temp$av_into) * temp$recovery

    } else {

      intControl <- df[, colnames(df) %in% samples(obj)[which(sampleGroups(obj) == recoveryFrom)]]
      intControl_av <- apply(intControl, MARGIN = 1, mean)
      intControl_sd <- apply(intControl, MARGIN = 1, sd)

      temp$recovery <- temp$av_into / intControl_av
      temp$recovery_sd <- ((temp$sd_into / temp$av_into) + (intControl_sd / intControl_av)) * temp$recovery

    }

    temp$isoN <- 0

    annot <- obj@annotation$comp[(obj@annotation$comp$ID %in% df$ID) & obj@annotation$comp$group %in% rg[g], ]

    IDs <- temp$ID

    for (f in seq_len(length(IDs))) {
      tempf <- annot$isogroup[annot$ID %in% IDs[f]]
      tempf <- obj@annotation$comp[obj@annotation$comp$isogroup %in% tempf, ]
      temp$isoN[f] <- length(unique(tempf$isoclass)) - 1
    }

    MS2 <- extractMS2(screen[which(analyses(screen) == sp), ], param = MS2param)

    temp <- mutate(temp,
                   hasFragments_Exp = FALSE,
                   fragments_mz_Exp = NA,
                   fragments_int_Exp = NA,
                   fragments_pre_Exp = NA)

    temp$nfrag <- 0

    temp$pfrag <- 0

    temp$intoscore <- 0

    for (i in seq_len(nrow(temp))) {

      xgroup <- temp$ID[i]

      xname <- temp$name[i]

      if (!is.na(xgroup)) {

        xMS2 <- MS2[[xgroup]]$MSMS

        if (!is.null(xMS2)) {

          temp$hasFragments_Exp[i] <- TRUE
          temp$fragments_mz_Exp[i] <- paste(xMS2$mz, collapse = ";")
          temp$fragments_int_Exp <- paste(xMS2$intensity, collapse = ";")
          temp$fragments_pre_Exp <- paste(xMS2$precursor, collapse = ";")
#TODO 220308 change back to temp$hasFragments ???
          if (temp$hasFragments_Exp[i]) {

            dbMS2 <- data.frame(mz = as.numeric(unlist(strsplit(temp$fragments_mz[i], split = ";"))),
                                intensity = as.numeric(unlist(strsplit(temp$fragments_int[i], split = ";"))),
                                precursor = as.logical(unlist(strsplit(temp$fragments_pre[i], split = ";"))))

            xMS2 <- top_n(xMS2, 10, intensity)
            xMS2 <- mutate(xMS2, into_ind = intensity / max(xMS2$intensity))

            dbMS2 <- top_n(dbMS2, 10, intensity)
            dbMS2 <- mutate(dbMS2, into_ind = intensity / max(dbMS2$intensity))

            combi <- fuzzyjoin::difference_inner_join(xMS2, dbMS2, by = c("mz"), max_dist = 0.02, distance_col = "diff")
            combi$diff <- (combi$diff/combi$mz.x) * 1E6
            combi <- combi[combi$diff < 10, ] #remove entries with max diff of 10 ppm

            temp$nfrag[i] <- nrow(combi)
            temp$pfrag[i] <- nrow(combi) / nrow(dbMS2)

            if (temp$nfrag[i] == 1) {
              temp$intoscore[i] <- 1
            } else {
              temp$intoscore[i] <- stats::cor(combi$into_ind.x,
                                              combi$into_ind.y,
                                              use = "everything",
                                              method = "pearson")
            }
          }
        }
      }
    }

    temp$pfrag <- round(temp$pfrag, digits = 2)

    temp$intoscore <- round(as.numeric(temp$intoscore), digits = 4)

    ls[[rg[g]]] <- temp

  }

  obj@IS@targets <- targets
  
  obj@IS@ppm <- ppm
  
  obj@IS@rtWindow <- rtWindow
  
  obj@IS@recoveryFrom <- recoveryFrom

  obj@IS@results <- ls

  if (exportResults) {

    results <- paste0(obj@path, "/results")
    if (!dir.exists(results)) dir.create(results)

    ls2 <- ls
    for (i in seq_len(length(ls))) ls2[[i]] <- ls2[[i]][, !colnames(ls2[[i]]) %in% samples(obj)] 
    ls2 <- rbindlist(ls2, idcol = "group")

    write.csv(ls2, file = paste0(results, "/IS01_df.csv"))

    plotis <- plotCheckIS(obj, rtWindow = rtWindow, ppm = ppm)

    widths <- c(8,5, rep(5, length(ls)))

    ggsave(paste0(results, "/IS02_DeviationsAndRecovery.tiff"),
           plot = plotis, device = "tiff", path = NULL, scale = 1,
           width = sum(widths), height = 10, units = "cm", dpi = 300, limitsize = TRUE)

    plotfp <- plotFeaturePeaks(obj = obj,
                               ID = df$ID,
                               mz = NULL,
                               rt = NULL,
                               rtUnit = "min",
                               ppm = NULL,
                               rtWindow = NULL,
                               names = df$name,
                               interactive = TRUE)

    saveWidget(partial_bundle(plotfp), file = paste0(results, "/IS03_featurePeaks.html"))

  }

  if (save) saveObject(obj = obj)
  
  return(obj)

}




#' @title plotCheckIS
#'
#' @param obj obj An \linkS4class{ntsData} object.
#' @param rtWindow The time window to show the time deviation.
#' @param ppm The ppm window to show the mass deviation.
#'
#' @return A plot of the IS check results.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#'
plotCheckIS <- function(obj, rtWindow = NULL, ppm = NULL) {

  if (length(obj@IS@results) < 1) {
    warning("No IS results found in the ntsdata object.")
    return(obj)
  }

  if (is.null(ppm)) ppm <- obj@IS@ppm
  
  if (is.null(rtWindow)) rtWindow <- obj@IS@rtWindow
  
  ls <- obj@IS@results

  maindf <- ls[[1]][, 1:7]

  if (nrow(maindf) == 0) {
    warning("IS results empty.")
    return(obj)
  }

  grobs <- list()

  grobs[["rt"]] <- ggplot(maindf) +
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
          axis.title.x = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 9)) +
    labs(title = "") +
    ylab("RT diff (sec)") +
    ylim(-rtWindow, rtWindow) +
    coord_flip()

  grobs[["ppm"]] <- ggplot(maindf) +
    theme_bw() +
    geom_rect(aes(ymin = -Inf, ymax = 5, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
    geom_rect(aes(ymin = 5, ymax = 10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
    geom_rect(aes(ymin = 10, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
    geom_point(stat = "identity", aes(x = name, y = d_ppm)) +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 7),
          axis.title.x = element_text(size = 7),
          plot.title = element_text(hjust = 0.5, size = 9)) +
    labs(title = "") +
    ylab("m/z diff (ppm)") +
    ylim(0, ppm) +
    coord_flip()

  for (i in seq_len(length(ls))) {

    max <- max((ls[[i]]$recovery * 100 + ls[[i]]$recovery_sd * 100) * 1.05)

    grobs[[names(ls)[i]]] <- ggplot(data = ls[[i]]) +
      theme_bw() +
      geom_rect(aes(ymin = 50, ymax = 150, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 10, ymax = 50, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
      geom_rect(aes(ymin = 150, ymax = 190, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
      geom_rect(aes(ymin = -Inf, ymax = 10, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_rect(aes(ymin = 190, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
      geom_point(stat = "identity", aes(x = name, y = recovery * 100)) +
      geom_errorbar(aes(x = name, y = recovery * 100,
                        ymin = recovery * 100 - recovery_sd * 100, ymax = recovery * 100 + recovery_sd * 100),
                    width = 0.2,
                    position = position_dodge(0.05)) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.text.x = element_text(size = 7),
            axis.title.x = element_text(size = 7),
            panel.grid.minor.x = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 9)) +
      labs(title = names(ls)[i]) +
      ylab("Recovery (%)") +
      scale_y_continuous(limits = c(0, ifelse(max < 200, 200, max))) +
      coord_flip()
  }

  return(grid.arrange(grobs = grobs,
                      ncol = length(grobs),
                      widths = c(8, rep(5, length(grobs) - 1))))

}
