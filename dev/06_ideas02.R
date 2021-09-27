#Copy template from template.R using x and/or use ?ntsIUTA for a tutorial.

#Run the following code to load the project setup:

### Setup -----

library(ntsIUTA)

setup <- ntsIUTA::makeSetup(projPath = getwd(), makeNewProject = F) # choose.dir()

#Getting sample list from setup
sampleInfo <- setup$sampleInfo
View(sampleInfo)

#Assign groups to replicate samples
sampleInfo$group <- c(rep("Blank_Plate",3),
                      rep("Blank_Plate_2",3),
                      rep("Blank",3),
                      rep("QC",5),
                      rep("Ref1_01",2),
                      rep("Ref2_025",2),
                      rep("Ref3_050",2),
                      rep("Ref4_1",2),
                      rep("Ref5_10",2),
                      rep("Ref6_100",2),
                      rep("Sample",3),
                      rep("Sample_Plate",3))

#Correcting blank assignment
sampleInfo$blank[base::grepl("Plate",sampleInfo$group)] <- "Blank_Plate"




### Import Raw Data -----

rawData <- ntsIUTA::importRawData(sampleInfo,
                                  centroidedData = TRUE,
                                  rtFilter = c(0.5,25), timeUnit = "min",
                                  removeEmptySpectra = TRUE,
                                  save = TRUE, projPath = setup$projPath)



### Parameters -----

#Parameters for peak picking
param = xcms::CentWaveParam(ppm = 15, peakwidth = c(6, 60),
                            snthresh = 5, prefilter = c(6, 300),
                            mzCenterFun = "wMean", integrate = 2,
                            mzdiff = -0.0001, fitgauss = TRUE,
                            noise = 0, verboseColumns = TRUE,
                            firstBaselineCheck = FALSE,
                            extendLengthMSW = TRUE)

#Only necessary if method for alignment is via PeaksGroups
param1 <- xcms::PeakDensityParam(sampleGroups = "holder",
                                 bw = 20,
                                 minFraction = 1,
                                 minSamples = 2,
                                 binSize = 0.005,
                                 maxFeatures = 100)

#Parameters for alignment of retention time across samples
#Not used if only one samples is given in peaksData
param2 <- xcms::PeakGroupsParam(minFraction = 1,
                                extraPeaks = 0,
                                smooth = "loess",
                                span = 0.2,
                                family = "gaussian")

#Parameters for final grouping of peaks across samples
param3 <- xcms::PeakDensityParam(sampleGroups = "holder",
                                 bw = 10,
                                 minFraction = 0.5,
                                 minSamples = 1,
                                 binSize = 0.005,
                                 maxFeatures = 100)





### Check QC replicate Group -----

#rawQC <- MSnbase::filterFile(rawData, base::which(rawData$sample_group == "QC"))

# rawQCNew <- MSnbase::filterFile(rawQC, 1:2)
# 
# rawQCOld <- MSnbase::filterFile(rawQC, 3:5)

#ran is script

#slQC <- utils::read.csv(utils::choose.files())

# QC_both <- ntsIUTA::checkQC(rawQC = rawQC,
#                             sampleInfo = sampleInfo,
#                             screeningList = slQC,
#                             paramPeaks = param,
#                             paramPreGrouping = param1,
#                             paramAlignment = param2,
#                             paramGrouping = param3,
#                             ppmForFillingGroups = 5,
#                             rtWindow = 30,
#                             ppmWindow = 15,
#                             polarity = "positive",
#                             plot = TRUE,
#                             save = TRUE,
#                             projPath = setup$projPath)

# QC_both$alignPlot
# plot(QC_both$evalPlot)
# View(QC_both$df)
# QC_both$featPlot


### PeakFinding -----

peaksData <- ntsIUTA::peakPicking(rawData,
                                  param = param,
                                  removeQC = TRUE,
                                  save = TRUE)





### Alignment and Grouping -----

#Filter for Blank2 and Sample
peaksDataSample <- peaksData[base::c("Blank_Plate", "Blank", "Sample", "Sample_Plate")]
featDataSample <- ntsIUTA::makeFeatures(peaksData = peaksDataSample,
                                        paramPreGrouping = param1,
                                        paramAlignment = param2,
                                        paramGrouping = param3,
                                        save = F)

#Filter for Reference Standards
peaksDataRef <- peaksData[base::unique(sampleInfo$group[base::grepl("Ref",sampleInfo$group)])]
featDataRef <- ntsIUTA::makeFeatures(peaksData = peaksDataRef,
                                     paramPreGrouping = param1,
                                     paramAlignment = param2,
                                     paramGrouping = param3,
                                     save = F)





### Apply Retention Time -----

plotAlignSample <- ntsIUTA::plotAlignment(featDataSample)
featDataSample <- xcms::applyAdjustedRtime(featDataSample)

plotAlignRef <- ntsIUTA::plotAlignment(featDataRef)
featDataRef <- xcms::applyAdjustedRtime(featDataRef)

#plotAlignSample

#plotAlignRef


### Get patRoon Data -----

patDataSample <- ntsIUTA::getPatData(featDataSample, sampleInfo[sampleInfo$group %in% base::c("Blank_Plate", "Blank", "Sample", "Sample_Plate"),])

patDataRef <- ntsIUTA::getPatData(featDataRef, sampleInfo[base::grepl("Ref",sampleInfo$group),])






### Check Internal Standards -----

#slIS <- utils::choose.files(base::getwd(), "Select the IS screening list")

ISSample <- ntsIUTA::checkIS(xPat = patDataSample,
                             screeningList = slIS)

# plot(ISSample$evalPlot)
# View(ISSample$df)
# ISSample$featPlot


ISRef <- ntsIUTA::checkIS(xPat = patDataRef,
                          screeningList = slIS,
                          save = FALSE)

# plot(ISRef$evalPlot)
# View(ISRef$df)
# ISRef$featPlot

# #only the samples from plate to confirm internal standards
# peaksDataSamplePlate <- peaksData[base::c("Blank_Plate", "Sample_Plate")]
# featDataSamplePlate <- ntsIUTA::makeFeatures(peaksData = peaksDataSamplePlate,
#                                              paramPreGrouping = param1,
#                                              paramAlignment = param2,
#                                              paramGrouping = param3,
#                                              save = F)
# ISPlate <- ntsIUTA::checkIS(x = featDataSamplePlate, save = FALSE)
# plot(ISPlate$evalPlot)
# View(ISPlate$df)






### Make Components -----

featCompSample <- ntsIUTA::makeFeatureComponents(featData =  featDataSample,
                                                 polarity = "positive",
                                                 sigma = 5, perfwhm = 0.45,
                                                 cor_eic_th = 0.85,
                                                 cor_exp_th = 0.85,
                                                 pval = 0.05,
                                                 calcCaS = TRUE,
                                                 calcIso = TRUE,
                                                 validateIsotopePatterns = TRUE,
                                                 ppmIsotopes = 40, mzabs = 0.01,
                                                 searchAdducts = TRUE,
                                                 ppmAdducts = 5, extendedList = FALSE,
                                                 excludeBlanks = TRUE,
                                                 blankGroups = c("Blank_Plate", "Blank"))

featCompRef <- ntsIUTA::makeFeatureComponents(featData =  featDataRef,
                                              polarity = "positive",
                                              sigma = 5, perfwhm = 0.45,
                                              cor_eic_th = 0.85,
                                              cor_exp_th = 0.85,
                                              pval = 0.05,
                                              calcCaS = TRUE,
                                              calcIso = TRUE,
                                              validateIsotopePatterns = TRUE,
                                              ppmIsotopes = 40, mzabs = 0.01,
                                              searchAdducts = TRUE,
                                              ppmAdducts = 5, extendedList = FALSE,
                                              excludeBlanks = F,
                                              blankGroups = "Blank")







### Make Feature List -----

flSample <- ntsIUTA::buildFeatureList(x = featDataSample, xA = featCompSample, xPat = patDataSample)

flRef <- ntsIUTA::buildFeatureList(x = featDataRef, xA = featCompRef, xPat = patDataRef)








### Filtering Feature List -----




#### Blank subtraction -----

#Features with intensity higher than the defined threshold remain in the feature list, below the threshold are removed.
blankThreshold <- 3

#Remove blanks for each sample, considering a threshold
removeFeatures <- flSample[,base::unique(featDataSample$sample_group)]
removeFeatures$remove <-  FALSE

for (r in 1:base::nrow(nSamples)) {
  for (i in 1:base::nrow(removeFeatures)) {
    if (removeFeatures[i,nSamples$group[r]] < blankThreshold * removeFeatures[i,nSamples$blank[r]]) {
      removeFeatures[i,nSamples$group[r]] <- 0
      removeFeatures[i,nSamples$blank[r]] <- 0
    }
  }
}
for (i in 1:base::nrow(removeFeatures)) {
  checkBlank <- base::unique(removeFeatures[i,base::unique(featDataSample$sample_group), drop = T])
  if (base::length(checkBlank) == 1) {
    if (checkBlank == 0) removeFeatures$remove[i] <- TRUE
  }
}

flSampleFiltered <- flSample
flSampleFiltered[,base::unique(featDataSample$sample_group)] <- dplyr::select(removeFeatures, -remove) 
flSampleFiltered <- flSampleFiltered[!removeFeatures$remove,] 







#### Replicates deviation -----
#zeros the intensity for features with high sd intensity is each replicate group

sdReplicateThreshold <- 0.4

rGroups <- base::unique(featDataSample$sample_group[!(featDataSample$sample_group %in% sampleInfo$blank)])

sdRepli <- flSampleFiltered[, base::colnames(flSampleFiltered) == base::paste(rGroups, "_sd", sep = "")]
for (i in 1:base::ncol(sdRepli)) {
  check <- sdRepli[,i] >= sdReplicateThreshold
  check[base::is.na(check)] <- FALSE
  flSampleFiltered[check, base::colnames(flSampleFiltered)[base::colnames(flSampleFiltered) == rGroups[i]]] <- 0
}

remove <- base::rep(FALSE,base::nrow(flSampleFiltered))
for (i in 1:base::nrow(flSampleFiltered)) {
  checkZeros <- base::unique(flSampleFiltered[i,base::unique(featDataSample$sample_group), drop = T])
  if (base::length(checkZeros) == 1) {
    if (checkZeros == 0) remove[i] <- TRUE
  }
}

flSampleFiltered <- flSampleFiltered[!remove,] 






#### Minimum intensity -----

minIntensityThreshold <- 600

rGroups <- base::unique(featDataSample$sample_group[!(featDataSample$sample_group %in% sampleInfo$blank)])
minIntRepli <- base::apply(flSampleFiltered[, base::colnames(flSampleFiltered) == rGroups], 1, function(x) max(x))
minIntRepli <- minIntRepli < minIntensityThreshold

flSampleFiltered <- flSampleFiltered[!minIntRepli,] 






#### Minimum SignalToNoise -----

minSN <- 3

flSampleFiltered <- flSampleFiltered[flSampleFiltered$sn_max >= minSN,] 






#### Visual Inspection -----


# ntsIUTA::plotFeatures(x = patDataSample,
#                       features = "M304_R659_11734",
#                       rtWindow = NULL,
#                       plotBy = "samples")



#### Count Features -----

nSamples <- sampleInfo[sampleInfo$group %in% base::unique(featDataSample$sample_group),]
nSamples <- nSamples[nSamples$group != nSamples$blank,]
nSamples <- nSamples[,c("group","blank")]
nSamples <- base::unique(nSamples)

#Count features in each sample replicate group
nSamples$nFeat <- 0
for (r in 1:base::nrow(nSamples)) {
  counter <- flSampleFiltered[,nSamples$group[r],drop=T] > 0
  nSamples$nFeat[r] <- base::length(counter[counter == TRUE])
}

View(nSamples)

### Suspect Screening -----


#### Reference Standards -----

#For reference standards
#slRef <- utils::read.csv(base::file.choose())
slRefStandards <- slRef[slRef$comment == "Reference",]

suspectsRef <- ntsIUTA::screenSuspectsFromCSV(patData = patDataRef,
                                              screeningList = slRefStandards,
                                              polarity = "positive", 
                                              ppmWindow = 20, rtWindow = 30,
                                              withMS2 = F, listMS2 = NULL)
suspectsRefDT <- suspectsRef$suspectsDT
View(suspectsRefDT)

#Collect Referene standards info
newslStandards <- slRefStandards

caliStandards <- newslStandards[,"name",drop=F]
caliStandards$rmse <- NA
caliStandards$group <- NA
conc <- base::c(0.1,0.25,0.5,1,10,100)
caliCurves <- base::list()
caliModels <- base::list()

for (i in 1:base::nrow(newslStandards)) {
  standardInt <- suspectsRefDT[suspectsRefDT$name %in% newslStandards$name[i],]
  if (base::nrow(standardInt) != 0) {
    standardInt <- standardInt[standardInt$Ref6_100 == base::max(standardInt$Ref6_100),]
    newslStandards$rt[i] <- standardInt$ret
    newslStandards$int10[i] <- standardInt$Ref5_10
    caliStandards$group[i] <- standardInt$group
    curve <- base::data.frame(conc = conc, int = base::unlist(standardInt[1, base::unique(featDataRef$sample_group), drop =T]))
    curve$sd <- base::unlist(flRef[flRef$patFT %in% standardInt$group, base::paste(base::unique(featDataRef$sample_group),"_sd", sep = ""), drop =T])
    curve$sd <- curve$sd * curve$int
    curve <- curve[curve$int > 1000,]
    if (base::nrow(curve) > 3) {
      #model <- stats::lm(int ~ conc, data = curve) # linear regression
      model <- stats::lm(log(int) ~ I((log(conc))^2) + I(log(conc)), data = curve) # wagner regression
      caliModels[[newslStandards$name[i]]] <- model
      model <- base::summary(model)
      caliStandards$rmse[i] <- base::round(model$adj.r.squared, digits = 5)
      caliCurves[[newslStandards$name[i]]] <- curve
      #curve$int <- base::with(curve, (int - base::min(int)) / (base::max(int) - base::min(int)))
    } else {
      caliCurves[[newslStandards$name[i]]] <- NA
      caliStandards$rmse[i] <- NA
    }
  } else {
    caliCurves[[newslStandards$name[i]]] <- NA
  }
}

View(caliCurves)

#plot calibrations and regresion lines
library(magrittr)

colors <- ntsIUTA::getColors(15)

title <- base::list(text = "Calibration Pesticides", x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))

xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Concentration (ng/ml)",
                    titlefont = base::list(size = 12, color = "black"))

yaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                    titlefont = base::list(size = 12, color = "black"))

caliPlot <- plotly::plot_ly()

for (i in 1:base::length(newslStandards$name)) {
  if (!base::is.na(caliCurves[newslStandards$name[i]]) & caliStandards$rmse[i] > 0.9) {
    caliPlot <- caliPlot %>% plotly::add_trace(caliCurves[[newslStandards$name[i]]],
                                               x = caliCurves[[newslStandards$name[i]]]$conc,
                                               y = caliCurves[[newslStandards$name[i]]]$int,
                                               type = "scatter", mode = "markers",
                                               marker = base::list(size = 5, color = colors[i]),
                                               error_y = base::list(array = caliCurves[[newslStandards$name[i]]]$sd, color = '#000000'),
                                               name = newslStandards$name[i], legendgroup = newslStandards$name[i])
    
    caliPlot <- caliPlot %>% plotly::add_lines(x = caliCurves[[newslStandards$name[i]]]$conc,
                                               y = base::exp(stats::fitted(caliModels[[newslStandards$name[i]]])),
                                               type = "scatter", mode = "markers",
                                               line = base::list(color =  colors[i], width = 1, dash = 'dash'),
                                               legendgroup = newslStandards$name[i], showlegend = F)
  }
}

caliPlot <- caliPlot %>% plotly::layout(legend = base::list(title = base::list(text='<b> Compound: </b>')),
                                        xaxis = xaxis,yaxis = yaxis, title = title)
caliPlot



##### New sl for References -----

# #Update the suspects data.frame with the real groups
# suspectsRefDT <- suspectsRefDT[suspectsRefDT$group %in% caliStandards$group,]
# 
# #Collect MS2 of standards and add to the screeningList
# control_avgPListParams <- patRoon::getDefAvgPListParams(
#   clusterMzWindow = 0.005,
#   topMost = 50,
#   minIntensityPre = 10,
#   minIntensityPost = 10)
# ms2Ref <- base::suppressWarnings(patRoon::generateMSPeakLists(
#   patDataRef[base::which(featDataRef$sample_group == "Ref6_100"), suspectsRefDT$group], "mzr",
#   maxMSRtWindow = 5,
#   precursorMzWindow = 3,
#   avgFeatParams = control_avgPListParams, 
#   avgFGroupParams = control_avgPListParams
# ))
# for (i in 1:base::nrow(suspectsRefDT)) {
#   xgroup <- suspectsRefDT$group[i]
#   xfrag <- ms2Ref[[xgroup]]$MSMS
#   if (!base::is.null(xfrag)) {
#     newslStandards$hasMS2[newslStandards$name %in% suspectsRefDT$name[i]] <- TRUE
#     newslStandards$mzMS2[newslStandards$name %in% suspectsRefDT$name[i]] <- base::paste(xfrag$mz, collapse = ";")
#     newslStandards$intMS2[newslStandards$name %in% suspectsRefDT$name[i]] <- base::paste(xfrag$intensity, collapse = ";")
#     newslStandards$preMS2[newslStandards$name %in% suspectsRefDT$name[i]] <- base::paste(xfrag$precursor, collapse = ";")
#   }
# }
# 
# #Add m/z for positive mode to the screeningList
# newslStandards$mz <- newslStandards$neutralMass + 1.007276
# newslStandards$rt <- newslStandards$rt /60
# 
# #save the screeningList in projPath
# utils::write.csv(newslStandards, file = base::paste0(setup$projPath,"/screeningList_PesticideReferenceStandards_MS2_pos.csv"))






#### Suspects in Samples -----

#For reference standards
#slPesticides <- utils::read.csv(base::file.choose())

suspectsSample <- ntsIUTA::screenSuspectsFromCSV(patData = patDataSample[, flSampleFiltered$patFT],
                                                 screeningList = slPesticides,
                                                 polarity = "positive", 
                                                 ppmWindow = 5, rtWindow = 30,
                                                 withMS2 = TRUE, listMS2 = NULL, ppmMS2 = 5,
                                                 removeBlanks = TRUE,
                                                 blankGroups = c("Blank_Plate", "Blank"))
suspectsSampleDT <- suspectsSample$suspectsDT
View(suspectsSampleDT)
utils::write.csv(suspectsSampleDT, file = base::paste0(setup$projPath,"/results/suspectsSampleDT.csv"))

#validate suspect screening with reference samples
suspectsRefVal <- ntsIUTA::screenSuspectsFromCSV(patData = patDataRef[base::which(featDataRef$sample_group == "Ref6_100"), ],
                                                 screeningList = slPesticides,
                                                 polarity = "positive", 
                                                 ppmWindow = 5, rtWindow = 30,
                                                 withMS2 = TRUE, listMS2 = NULL, ppmMS2 = 5,
                                                 removeBlanks = TRUE,
                                                 blankGroups = c("Blank_Plate", "Blank"))
suspectsRefValDT <- suspectsRefVal$suspectsDT
View(suspectsRefValDT)




#### Suspects via Fragments -----

#slFragments <- utils::read.csv(utils::choose.files())

suspectsFragments <- ntsIUTA::screenSuspectsFromCSV(patData = patDataSample[, flSampleFiltered$patFT],
                                                    screeningList = slFragments,
                                                    polarity = "positive", 
                                                    ppmWindow = 5, rtWindow = 10,
                                                    withMS2 = TRUE, listMS2 = NULL, ppmMS2 = 10,
                                                    removeBlanks = TRUE,
                                                    blankGroups = c("Blank_Plate", "Blank"))

suspectsFragmentsDT <- suspectsFragments$suspectsDT

View(suspectsFragmentsDT)
utils::write.csv(suspectsFragmentsDT, file = base::paste0(setup$projPath,"/results/suspectsFragmentsDT.csv"))



control_avgPListParams <- patRoon::getDefAvgPListParams(
  clusterMzWindow = 0.005,
  topMost = 50,
  minIntensityPre = 10,
  minIntensityPost = 10
)

MS2 <- base::suppressWarnings(patRoon::generateMSPeakLists(
  patDataSample[which(featDataSample$sample_group == c("Sample_Plate")), "M435_R1142_23178"],
  "mzr",
  maxMSRtWindow = 5,
  precursorMzWindow = 3,
  avgFeatParams = control_avgPListParams, 
  avgFGroupParams = control_avgPListParams
))

MS2[["M435_R1142_23178"]]$MSMS



ntsIUTA::plotComponentSpectrum(
  xA = featCompSample,
  replicateGroups = c("Sample_Plate"),
  features = "M323_R746_13205",
  featData = patDataSample,
  mz = NULL,
  ppm = 5,
  mzWindowPlot = c(320,330),
  rt = NULL,
  rtWindow = 1,
  rtUnit = "min",
  comp = NULL,
  onlyAnnotated = FALSE,
  onlyRelated = TRUE,
  intval = "maxo",
  log = TRUE
)


ntsIUTA::plotFeatures(x = patDataSample,
                      features = "M323_R746_13205")

ntsIUTA::plotFeatures(x = patDataSample,
                      features = "M435_R1142_23178")

# 
# ntsIUTA::plotComponentSpectrum(xA = featCompSample,
#                                features = "M323_R746_10641",
#                                featData = patDataSample, log = F)
# 
# ntsIUTA::plotRawChrom(x = rawData,
#                       fileIndex = base::which(sampleInfo$group %in% c("Blank_Plate", "Sample_Plate")))




### Evaluate Features in Plate -----



#### Prioritization -----

#Similar/Different
flEval <- flSampleFiltered

nFeatures <- base::data.frame(group = flEval$patFT,
                              inSample = base::rep(FALSE, base::nrow(flEval)),
                              inSample_Plate = base::rep(FALSE, base::nrow(flEval)),
                              inBoth = base::rep(FALSE, base::nrow(flEval)))

nFeatures <- base::cbind(nFeatures, flEval[, base::unique(featDataSample$sample_group)])

intThresholdEval <- 600
relThresholdEval <- 20/100 #1%

nFeatures$inSample[nFeatures$Sample > 0 & (nFeatures$Sample_Plate < intThresholdEval | nFeatures$Sample_Plate < nFeatures$Sample*relThresholdEval)] <- TRUE
nFeatures$inSample_Plate[nFeatures$Sample_Plate > 0 & (nFeatures$Sample < intThresholdEval | nFeatures$Sample < nFeatures$Sample_Plate*relThresholdEval)] <- TRUE
#the remaining
nFeatures$inBoth[nFeatures$inSample == FALSE & nFeatures$inSample_Plate == FALSE] <- TRUE

#nFeatures$inBoth[nFeatures$Sample_Plate > 300 & nFeatures$Sample > 300] <- TRUE

#count Features
nFeaturesSum <- base::data.frame(cat = c("inSample","inSample_Plate","inBoth"), number = NA)
nFeaturesSum$number[1] <- base::length(nFeatures$inSample[nFeatures$inSample == TRUE])
nFeaturesSum$number[2] <- base::length(nFeatures$inSample_Plate[nFeatures$inSample_Plate == TRUE])
nFeaturesSum$number[3] <- base::length(nFeatures$inBoth[nFeatures$inBoth == TRUE])
View(nFeaturesSum)
base::sum(nFeaturesSum$number) == base::nrow(nFeatures)


flInSample_Plate <- flEval[nFeatures$inSample_Plate,]

flInSample <- flEval[nFeatures$inSample,]
flInSample <- utils::head(dplyr::arrange(flInSample, dplyr::desc(Sample)), n = base::nrow(flInSample)*0.2) #top 20%


#### Assessment of Oxydation -----

#tracking the oxidation of compounds back to the original sample, to explain the lower number of comounds in both samples

#database = utils::read.csv(utils::choose.files())

source("transformationsPathways.R")
dfTransformations <- searchTransformationsPathways_temp(x = flInSample_Plate,
                                                        y = flInSample,
                                                        rgX = "Sample_Plate",
                                                        rgY = "Sample",
                                                        rgBlank = c("Blank_Plate","Blank"),
                                                        database = database,
                                                        transformations = NULL,
                                                        ppmWindow = 5,
                                                        ppmMS2 = 10,
                                                        withMS2 = TRUE,
                                                        xPat = patDataSample)

#Count transformation categories for overall prespective of possible
countTransformations <- dplyr::count(dfTransformations[,"FD",drop=F], FD)
countTransformations$type <- base::ifelse(base::grepl(c("\\+O|\\+2O|\\+3O|\\+O2|\\+4O"), countTransformations$FD, fixed = F), "Oxydation", "Other")
countTransformationsSum <- dplyr::count(countTransformations[,"type",drop=F], type)

#How many are oxydations
countTransformationsSum$n[countTransformationsSum$type == "Oxydation"]/countTransformationsSum$n[countTransformationsSum$type == "Other"]*100
library(magrittr)
plotly::plot_ly() %>% plotly::add_bars(x = base::factor(countTransformations$FD),
                                       y = countTransformations$n,
                                       color = countTransformations$type,
                                       colors = c(Other = "gray", Oxydation = "forestgreen"))













### Other code ----

# ##### Validation with IS
# #Find oxydations of cyclophosphamide d6, diuron d3 and Diclophenac d4 from Sample in Sample Plate
# 
# View(ISSample$df)
# 
# valSampleIS <- ISSample$df[ISSample$df$name %in% c("Cyclophosphamide d6","Diuron d3","Diclophenac d4"),]
# valSampleIS <- valSampleIS[valSampleIS$d_rt < 15, "group", drop = T]
# valSampleIS <- flSample[flSample$patFT %in% valSampleIS,]
# 
# source("transformationsPathways.R")
# valTransformations <- searchTransformationsPathways_temp(x = flSample,
#                                                         y = valSampleIS,
#                                                         rgX = "Sample_Plate",
#                                                         rgY = "Sample",
#                                                         rgBlank = c("Blank_Plate","Blank"),
#                                                         database = database,
#                                                         transformations = c("\\+O|\\+2O|\\+3O|\\+O2|\\+4O"),
#                                                         ppmWindow = 5,
#                                                         ppmMS2 = 10,
#                                                         withMS2 = TRUE,
#                                                         xPat = patDataSample)





