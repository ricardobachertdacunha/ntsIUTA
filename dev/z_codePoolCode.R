



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
# 
# mz <- 407.1966
# rt <- 22.25
# 
# #ntsIUTA::plotTargetCentroids(rawData = rawData_cent, fileIndex = c(1,2), mz = mz, ppmWindow = 20, rt = rt, rtWindow = 1, rtUnit = "min")
# #rawData_cent_2 <- rawData_cent %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 6, polynomialOrder = 3)
# #rawData_cent_2 <- rawData_cent %>% MSnbase::smooth(method = "MovingAverage", halfWindowSize = 2, weighted = T)
# #ntsIUTA::plotTargetCentroids(rawData = rawData_cent_2, fileIndex = c(1,2), mz = 287.0582, ppmWindow = 20, rt = 22.8, rtWindow = 1, rtUnit = "min")
# 
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
# 
# 
# rawData3 <- rawData %>% MSnbase::smooth(method = "SavitzkyGolay", halfWindowSize = 2, polynomialOrder = 3)
# #rawData3 <- MSnbase::pickPeaks(rawData3, halfWindowSize = 3, refineMz = "kNeighbors", k = 1)
# rawData3 <- MSnbase::pickPeaks(rawData3, halfWindowSize = 10, refineMz = "descendPeak")
# ntsIUTA::plotTargetCentroids(rawData = rawData3, fileIndex = c(1,2), mz = mz, ppmWindow = 30, rt = rt, rtWindow = 1, rtUnit = "min")
