

# library(ntsIUTA)

# projPath <- system.file(package = "ntsIUTA", dir = "extdata")
# sampleInfo <- setupProject(projPath, save = FALSE, makeNewProject = FALSE)
# sampleInfo <- sampleInfo[1:6, ]

# rawData <- importRawData(sampleInfo,
#                          rtFilter = c(13, 17),
#                          timeUnit = "min",
#                          centroidedData = NA,
#                          removeEmptySpectra = TRUE,
#                          save = FALSE)

# EIC_diuron <- ntsIUTA::extractEIC(rawData, fileIndex = 4, mz = 233.0243, ppm = 20)

# plotRawChrom(rawData, fileIndex = 4, type = "tic", mz = 233.0243, ppm = 20, rt = 15, rtWindow = 2, rtUnit = "min")

# projPath <- system.file(package = "ntsIUTA", dir = "extdata")

# sI <- ntsIUTA::setupProject(projPath, save = FALSE, makeNewProject = FALSE)

# sI <- sI[7:9, ]

# rD <- ntsIUTA::importRawData(sI,
#                              rtFilter = c(c(13, 17)),
#                              timeUnit = "min",
#                              centroidedData = FALSE,
#                              removeEmptySpectra = TRUE,
#                              save = FALSE)

# ntsIUTA::plotTargetCentroids(rD, fileIndex = 1,
#                              mz = 242.1434, ppm = 90,
#                              rt = 14.65, rtWindow = 0.5,
#                              rtUnit = "min", plotTargetMark = TRUE)

# rawData_prof2 <- centroidProfileData(MSnbase::filterFile(rawData_prof,1),
#                                      smoothing = FALSE, refineMZ = TRUE,
#                                      methodRefineMz = "kNeighbors", k = 1,
#                                      save = FALSE)


# x <- rD %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 2, polynomialOrder = 3)
# x <- MSnbase::pickPeaks(x)

# x <- MSnbase::pickPeaks(rD)

# x <- MSnbase::pickPeaks(rD,
#                         halfWindowSize = 1,
#                         SNR = 0,
#                         method = "SuperSmoother",
#                         refineMz = "none")

# x <- MSnbase::pickPeaks(rD,
#                         halfWindowSize = 0.8,
#                         SNR = 0,
#                         method = "SuperSmoother",
#                         refineMz = "kNeighbors",
#                         k = 2)

# x <- x %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 2, polynomialOrder = 3)


# x <- MSnbase::pickPeaks(rD,
#                         halfWindowSize = 20,
#                         SNR = 0,
#                         method = "SuperSmoother",
#                         refineMz = "descendPeak",
#                         signalPercentage = 0.3,
#                         stopAtTwo = TRUE)

# ntsIUTA::plotTargetCentroids(x, fileIndex = 1,
#                              mz = 242.1434, ppm = 30,
#                              rt = 14.65, rtWindow = 0.5,
#                              rtUnit = "min",
#                              plotTargetMark = TRUE)

























### Code -----

# TODO make function to add MS2 data to a screening list, based on suspect screening workflow.

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









#' For \pkg{patRoon} the rawData should be the \code{sampleInfo} returning the respectives    and returns a \code{list} of \linkS4class{XCMSnExp}
#' objects for each replicate sample group as defined in the \code{\link{setupProject}} function.
#' The \linkS4class{XCMSnExp} objects can be conconated via \code{c()} function of the \pkg{xcms} package.
#' We separate the peaks from each replicate group to facilitate workflows that include multi-project cross-analysis.
#' The peak picking uses the function \code{\link[xcms]{chromatographic-peak-detection}} from the \pkg{xcms} package.
#' Different methods for peak picking can be used. For more information, see documentation of \code{\link[xcms]{chromatographic-peak-detection}}.

# sampleInfo <- ntsIUTA::setupProject(save = FALSE)
# View(sampleInfo)
# sampleInfo <- sampleInfo[19:21, ]
# View(sampleInfo)
# rawData <- ntsIUTA::importRawData(sampleInfo = sampleInfo, save = FALSE, centroidedData = FALSE, removeEmptySpectra = FALSE)
# library(magrittr)
# rawData2 <- rawData %>% MSnbase::filterRt(rt = c(13 * 60, 17 * 60), msLevel. = c(1, 2)) %>% MSnbase::filterMz(mz = c(230, 300), msLevel. = c(1))
# MSnbase::writeMSData(rawData2, file = c("Sample04_prof.mzML", "Sample05_prof.mzML", "Sample06_prof.mzML"), copy = TRUE)


#Options to save information from raw data
  #raw[["rawchrom"]] <- chromatogram(rawData, aggregationFun = "max")
  #raw[["emptyspec"]] <- length(which(MSnbase::peaksCount(rawData) == 0))

  #For examples afterwards
  # file <- dir(system.file(package = "MSnbase", dir = "extdata"),
  #             full.name = TRUE,
  #             pattern = "mzXML$")

  #Verify if mzML files are centroided ot not
  #TestDataTypeCent <- table(MSnbase::isCentroided(MSnbase::filterFile(rawData_cent, 1)), useNA = "always")
  #TestDataTypeCent <- table(MSnbase::isCentroidedFromFile(MSnbase::filterFile(rawData_cent, 1)), useNA = "always")



# @title paramListExample
#
# @description A \code{list} containing the parameters needed for \code{\link{peakPicking}} and \code{\link{makeFeatures}}.
#
# @format A \code{list} with following entries:
# #' \describe{
#   \item{instName}{The name of the list. Recommended to reflect the instrumentation and conditions used for data acquisition.}
#   \item{PP}{A \code{list} with the method specific parameters used for peak picking.
#   See \code{\link[xcms]{chromatographic-peak-detection}} for more information and which parameters can be added.}
#   \item{preGrouping}{A \code{list} with the method specific parameters used for pre-grouping, which is necessary
#   for \code{\link[xcms]{adjustRtime}} function when using the method \code{PeakGroups} as it requires peaks groups (\emph{i.e.} features)
#   for alignment. See \code{\link[xcms]{groupChromPeaks-density}} for more information}
#   \item{alignment}{A \code{list} with the method specific parameters used for retention time alignment.
#   See \code{\link[xcms]{adjustRtime}} for more information on the possible methods and parameters that can be used.}
#  \item{grouping}{A \code{list} with the method specific parameters used for grouping peaks.
#   See \code{\link[xcms]{groupChromPeaks}} for more information on the possible methods and parameters that can be used.}
# }
#
# Note: For all lists a \code{funcName} entry should be given matching the expected class object (\emph{e.g.} \code{MassifquantParam} from
# \code{\link[xcms]{findChromPeaks-massifquant}}).
#
# @source
#"paramListExample"


# ft <- "FT001"
# chrom <- xcms::featureChromatograms(featData, feature = ft, include = "feature_only", filled = TRUE)
# chrom
# hasFilledChromPeaks(chrom)
# chromP <- chromPeaks(chrom)
# chromP
# colors <- RColorBrewer::brewer.pal.info[RColorBrewer::brewer.pal.info$category == 'qual',]
# colors <- unlist(mapply(RColorBrewer::brewer.pal, colors$maxcolors, rownames(colors)))
# colors <- paste0(colors[1:nrow(setup$sampleInfo)],"60")
# names(colors) <- setup$sampleInfo$analysis
# plot(chrom, col = colors, peakBg = colors[chromP[, "sample"]])



#defFT <- base::as.data.frame(xcms::featureDefinitions(featDataExample))
#FT <- "FT106"
#defFT <- defFT[base::row.names(defFT) %in% FT,]
#chromPeaks <- as.data.frame(xcms::chromPeaks(featData)[unlist(defFT[base::row.names(defFT) %in%  FT,"peakidx"]), ])
#chromPeaksX <- as.data.frame(xcms::peaks(xSet)[unlist(defFT[base::row.names(defFT) %in%  FT,"peakidx"]), ])

#test <- CAMERA::getPeaklist(xA1)
#testX <-base::row.names(xcms::featureDefinitions(featData))
#test <- dplyr::filter(test, pcgroup == 15)


# CAMERA::getPeaklist(xA[[rIndex]])
# test <- dplyr::filter(test, pcgroup == 7)
# testx <- as.data.frame(xcms::featureDefinitions(featData))
# #testX <-base::row.names(xcms::featureDefinitions(featData))





#setup <- ntsIUTA::makeSetup(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)

#rawData <- ntsIUTA::importRawData(setup$sampleInfo[1:3, ], save = FALSE, centData = T)






# rawData_cent <- ntsIUTA::importRawData(setup$sampleInfo[3:4, ], save = FALSE, centData = T)
# rawData <- ntsIUTA::importRawData(setup$sampleInfo[5:6, ], save = FALSE, centData = F)

# mz <- 407.1966
# rt <- 22.25

# #ntsIUTA::plotTargetCentroids(rawData = rawData_cent, fileIndex = c(1,2), mz = mz, ppmWindow = 20, rt = rt, rtWindow = 1, rtUnit = "min")
# #rawData_cent_2 <- rawData_cent %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 6, polynomialOrder = 3)
# #rawData_cent_2 <- rawData_cent %>% MSnbase::smooth(method = "MovingAverage", halfWindowSize = 2, weighted = T)
# #ntsIUTA::plotTargetCentroids(rawData = rawData_cent_2, fileIndex = c(1,2), mz = 287.0582, ppmWindow = 20, rt = 22.8, rtWindow = 1, rtUnit = "min")
#
# #Test Smoothing
#
# ntsIUTA::plotTargetCentroids(rawData = rawData_cent, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Proteo")
#
# rawData0 <- MSnbase::pickPeaks(rawData)
# ntsIUTA::plotTargetCentroids(rawData = rawData0, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Control")
#
# rawData1 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 2, polynomialOrder = 3)
# rawData1 <- MSnbase::pickPeaks(rawData1)
# ntsIUTA::plotTargetCentroids(rawData = rawData1, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Test 1")
#
# rawData2 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 2, polynomialOrder = 4)
# rawData2 <- MSnbase::pickPeaks(rawData2)
# ntsIUTA::plotTargetCentroids(rawData = rawData2, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Test 2")
#
# rawData3 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 4, polynomialOrder = 3)
# rawData3 <- MSnbase::pickPeaks(rawData3)
# ntsIUTA::plotTargetCentroids(rawData = rawData3, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Test 3")
#
# rawData4 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 4, polynomialOrder = 4)
# rawData4 <- MSnbase::pickPeaks(rawData4)
# ntsIUTA::plotTargetCentroids(rawData = rawData4, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Test 4")
#
# rawData5 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 4, polynomialOrder = 4)
# rawData5 <- MSnbase::pickPeaks(rawData5, halfWindowSize = 4, refineMz = "descendPeak")
# ntsIUTA::plotTargetCentroids(rawData = rawData5, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min", title = "Test 5")
#
#
#
# rawData3 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 2, polynomialOrder = 3)
# #rawData3 <- MSnbase::pickPeaks(rawData3, halfWindowSize = 3, refineMz = "kNeighbors", k = 1)
# rawData3 <- MSnbase::pickPeaks(rawData3, halfWindowSize = 10, refineMz = "descendPeak")
# ntsIUTA::plotTargetCentroids(rawData = rawData3, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min")
