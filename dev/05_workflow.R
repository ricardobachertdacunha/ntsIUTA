
### Load ntsIUTA --------------------------------------------------------------------------------------------

library(ntsIUTA)




### Project Setup -------------------------------------------------------------------------------------------

path <- system.file(package = "ntsIUTA", dir = "extdata")

dtall <- setupProject(path, polarity = "type", save = FALSE, makeNewProject = FALSE)




#### S4 Methods ------------------------------------------------------------------------------------------

##### Sample Groups -----

#Assign/Correct sample replicate groups
sampleGroups(dtall) <- c(rep("Blank", 3),
                         rep("IN", 3),
                         rep("OZ", 3),
                         rep("UV", 3),
                         rep("AC", 3),
                         rep("Centroid", 3),
                         rep("Profile", 3))

#See sample replicate groups
sampleGroups(dtall)

#see the samples (files)
samples(dtall) #no setter




##### Blanks -----

#Assign the blank replicate sample group/s
blanks(dtall) <- "Blank"

#see the blank replicate group/s
blanks(dtall)




##### QC -----

#When needed, assigning the QC sample replicate group
#QC(dtall) <- "QC"

#see the QC samples
QC(dtall)




##### show ntsData -----

#show method
dtall




##### sub-setting -----

#sub-setting by samples
dtcent <- dtall[16:18]
dtprof <- dtall[19:21]



# TODO Add accessors to slots, such as title, path, etc.

# TODO show the addFiles function

# TODO show the convert files function





### Import Data ---------------------------------------------------------------------------------------------

dtcent <- importRawData(dtcent,
                        rtFilter = c(13.8, 16.3),
                        rtUnit = "min",
                        centroidedData = NA,
                        removeEmptySpectra = TRUE,
                        save = FALSE)




#### TIC, BPC and EIC ---------------------------------------------------------------------------------------

exEIC <- extractEIC(dtcent,
                    samples = 1:2,
                    mz = 242.1434, ppm = 20,
                    rt = 14.8, rtWindow = 0.5,
                    rtUnit = "min")




#### Inspecting Raw Data ------------------------------------------------------------------------------------

plotRawChrom(dtcent,
             samples = NULL,
             type = "tic",
             mz = 242.1434, ppm = 20,
             rt = 14.8, rtWindow = 0.5,
             colorBy = "samples",
             rtUnit = "min",
             interactive = FALSE)

#iterative
plotRawChrom(dtcent,
             samples = NULL,
             type = "tic",
             mz = 242.1434, ppm = 20,
             rt = 14.8, rtWindow = 0.5,
             rtUnit = "min",
             colorBy = "samples",
             interactive = TRUE)

#plot correlation of replicate sample groups
plotCorReplicates(dtcent)

#plot centroids
plotTargetCentroids(obj = dtcent,
                    samples = NULL,
                    mz = 242.1434, ppm = 20,
                    rt = 14.8, rtWindow = 0.5,
                    rtUnit = "min", plotTargetMark = FALSE)


### Profile vs Centroid -------------------------------------------------------------------------------------

dtprof <- importRawData(dtprof,
                        rtFilter = c(13.8, 16.3),
                        rtUnit = "min",
                        centroidedData = FALSE,
                        removeEmptySpectra = TRUE,
                        save = FALSE)

# TODO Is centroided data working or not working
#check using raw data without centroids
table(MSnbase::isCentroidedFromFile(dtcent@MSnExp), MSnbase::msLevel(dtcent@MSnExp))
table(MSnbase::isCentroidedFromFile(dtprof@MSnExp), MSnbase::msLevel(dtprof@MSnExp))

# plot centroided and profile data
p1 <- plotTargetCentroids(dtcent, samples = 1,
                          mz = 242.1434, ppm = 150,
                          rt = 14.8, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

p2 <- plotTargetCentroids(dtprof, samples = 1,
                          mz = 242.1434, ppm = 150,
                          rt = 14.8, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

plotly::subplot(list(p1, p2), nrows = 1, margin = 0.04)




### Centroiding Data -------------------------------------------------------------------------------------------

dtprof <- centroidProfileData(obj = dtprof,
                              halfwindow = 3,
                              SNR = 3,
                              noiseMethod = "MAD",
                              methodRefineMz = "kNeighbors",
                              k = 1,
                              smoothing = FALSE,
                              save = FALSE)


p1 <- plotTargetCentroids(dtcent, samples = 1,
                          mz = 242.1434, ppm = 30,
                          rt = 14.75, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

p2 <- plotTargetCentroids(dtprof, samples = 1,
                          mz = 242.1434, ppm = 30,
                          rt = 14.75, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

plotly::subplot(list(p1, p2), nrows = 1, margin = 0.04)




### Peak Picking --------------------------------------------------------------------------------------------
# TODO when applying rtFilter for import raw data reuse when peak picking with patRoon

dt <- dtall[1:6]

dt <- importRawData(dt,
                    rtFilter = NULL,
                    rtUnit = "min",
                    centroidedData = NA,
                    removeEmptySpectra = TRUE,
                    save = FALSE)

plotRawChrom(dt)

plotCorReplicates(dt, binsize = 3)




#### xcms3 --------------------------------------------------------------------------------------------------

dtxcms <- dt

dtxcms@algorithms$peakPicking <- "xcms3"

paramXCMS <- xcms::CentWaveParam(
  ppm = 15, peakwidth = c(6, 60),
  snthresh = 5, prefilter = c(6, 300),
  mzCenterFun = "mean", integrate = 2,
  mzdiff = -0.0001, fitgauss = TRUE,
  noise = 0, verboseColumns = TRUE,
  firstBaselineCheck = FALSE,
  extendLengthMSW = TRUE)

dtxcms@parameters$peakPicking <- paramXCMS

dtxcms <- peakPicking(obj = dtxcms, save = FALSE)




#### openms -------------------------------------------------------------------------------------------------

dtopenms <- dt

dtopenms@algorithms$peakPicking <- "openms"

paramOpenms <- list(
  noiseThrInt = 250,
  chromSNR = 3,
  chromFWHM = 5,
  mzPPM = 15,
  reEstimateMTSD = TRUE,
  traceTermCriterion = "sample_rate",
  traceTermOutliers = 5,
  minSampleRate = 1,
  minTraceLength = 3,
  maxTraceLength = -1,
  widthFiltering = "fixed",
  minFWHM = 5,
  maxFWHM = 60,
  traceSNRFiltering = FALSE,
  localRTRange = 10,
  localMZRange = 6.5,
  isotopeFilteringModel = "metabolites (5% RMS)",
  MZScoring13C = FALSE,
  useSmoothedInts = TRUE,
  extraOpts = NULL,
  intSearchRTWindow = 3)

dtopenms@parameters$peakPicking <- paramOpenms

dtopenms <- peakPicking(obj = dtopenms, save = FALSE)




#### sirius -------------------------------------------------------------------------------------------------
#TODO Error with mz column

dtsirius <- dt

dtsirius@algorithms$peakPicking <- "sirius"

dtsirius@parameters$peakPicking <- list()

dtsirius <- peakPicking(obj = dtsirius, save = FALSE)




#### kpic2 -------------------------------------------------------------------------------------------------

dtkpic2 <- dt

dtkpic2@algorithms$peakPicking <- "kpic2"

dtkpic2@parameters$peakPicking <- list(
  kmeans = TRUE,
  level = 1000,
  mztol = 0.01,
  gap = 3,
  width = c(5, 60),
  alpha = 0.3,
  min_snr = 4,
  parallel = TRUE
)

dtkpic2 <- peakPicking(obj = dtkpic2, save = FALSE)
class(dtkpic2@patdata)




#### safd -------------------------------------------------------------------------------------------------

dtsafd <- dtprof

dtsafd@algorithms$peakPicking <- "safd"

dtsafd@parameters$peakPicking <- list(
  profPath = dirname(dtprof@samples$file),
  mzRange = c(0, 1200),
  maxNumbIter = 1000,
  maxTPeakW = 300,
  resolution = 20000,
  minMSW = 0.02,
  RThreshold = 0.75,
  minInt = 500,
  sigIncThreshold = 5,
  S2N = 3,
  minPeakWS = 3
)

dtsafd <- peakPicking(obj = dtsafd, save = FALSE)




#### Inspecting Peaks ---------------------------------------------------------------------------------------

#Option 1
dtxcms@peaks[1,]
dtopenms@peaks[1,]

#Option 2, using the S4 method
peaks(dtxcms, samples = NULL, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")

plotPeaks(dtxcms, samples = NULL, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min", colorBy = "peaks")

View(peaks(dtxcms, samples = NULL, mz = 441.1670, ppm = 10, rt = 15.27, rtUnit = "min"))

peak1 <- plotPeaks(dtxcms, samples = 4, mz = 441.1670, ppm = 10, rt = 15.27, rtUnit = "min", colorBy = "peaks")
peak2 <- plotPeaks(dtxcms, samples = 5, mz = 441.1670, ppm = 10, rt = 15.27, rtUnit = "min", colorBy = "peaks")
peak3 <- plotPeaks(dtxcms, samples = 6, mz = 441.1670, ppm = 10, rt = 15.27, rtUnit = "min", colorBy = "peaks")
plotly::subplot(list(peak1, peak2, peak3), nrows = 1, margin = 0.04)

mapPeaks(dtxcms, mz = 332.2200, ppm = 20, rtWindow = c(1000, 1200), rtUnit = "sec")

mapPeaks(dtxcms, samples = c(1, 4), mz = c(234, 270), ppm = NULL, rtWindow = c(920, 940), rtUnit = "sec")




### Alignment and Grouping ----------------------------------------------------------------------------------

#### xcms3 --------------------------------------------------------------------------------------------------

dtxcms@algorithms$makeFeatures <- "xcms3"

paramGroupingxcms <- list(
  rtalign = TRUE,
  loadRawData = TRUE,
  groupParam = xcms::PeakDensityParam(
    sampleGroups = "holder",
    bw = 3,
    minFraction = 0.5,
    minSamples = 1,
    binSize = 0.008,
    maxFeatures = 100),
  preGroupParam = xcms::PeakDensityParam(
    sampleGroups = "holder",
    bw = 5,
    minFraction = 0.5,
    minSamples = 1,
    binSize = 0.008,
    maxFeatures = 100),
  retAlignParam = xcms::PeakGroupsParam(
    minFraction = 1,
    extraPeaks = 0,
    smooth = "loess",
    span = 0.3,
    family = "gaussian")
)

dtxcms@parameters$peakGrouping <- paramGroupingxcms

paramFill <- xcms::FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)

paramFill <- xcms::ChromPeakAreaParam(
  mzmin = function(z) quantile(z, probs = 0.25),
  mzmax = function(z) quantile(z, probs = 0.75),
  rtmin = function(z) quantile(z, probs = 0.25),
  rtmax = function(z) quantile(z, probs = 0.75)
)

dtxcms@parameters$fillMissing <- paramFill

dtxcms

dtxcms2 <- makeFeatures(obj = dtxcms,
                        algorithm = NULL,
                        paramGrouping = NULL,
                        recursive = TRUE,
                        paramFill = NULL,
                        save = FALSE)




#### openms -------------------------------------------------------------------------------------------------

dtopenms@algorithms$makeFeatures <- "openms"

paramGroupingOpenms <- list(
  rtalign = TRUE,
  QT = FALSE,
  maxAlignRT = 5,
  maxAlignMZ = 0.008,
  maxGroupRT = 3,
  maxGroupMZ = 0.008,
  extraOptsRT = NULL,
  extraOptsGroup = NULL
)

dtopenms@parameters$peakGrouping <- paramGroupingOpenms

dtopenms2 <- makeFeatures(obj = dtopenms,
                         algorithm = NULL,
                         paramGrouping = NULL,
                         recursive = FALSE,
                         paramFill = NULL,
                         save = FALSE)




#### Ploting Alignment --------------------------------------------------------------------------------------

# TODO check plotAlignment, may only work when xcms3 is used as algorithm to align peaks
plotAlignment(features)




#### Inspecting Features ------------------------------------------------------------------------------------


##### Getter -----

#Direct access to slot
dtxcms2@features[1,]

#S4 accessor method
features(dtxcms2, samples = NULL, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")

##### plot Chroms -----

#Interactive feature/s plot
plotFeatures(obj = dtxcms2,
             samples = NULL,
             ID = NULL,
             mz = 213.1869,
             rt = 15.47,
             rtUnit = "min",
             ppm = NULL,
             rtWindow = NULL,
             colorBy = "features",
             interactive = TRUE)

#Static feature/s plot
plotFeatures(obj = dtxcms2,
             samples = NULL,
             ID = NULL,
             mz = 213.1869,
             rt = 15.47,
             rtUnit = "min",
             ppm = NULL,
             rtWindow = NULL,
             colorBy = "features",
             interactive = FALSE)


plotFeatures(obj = dtxcms2,
             samples = 3:4,
             ID = NULL,
             mz = 332.2200,
             rt = NULL,
             rtUnit = "sec",
             ppm = 20,
             rtWindow = c(900, 1200),
             colorBy = "samples",
             interactive = FALSE)

#patRoon option
patRoon::plotChroms(dtxcms2@patdata[, c("M263_R937_1829")])


##### plot Feature Peaks -----

plotFeaturePeaks(obj = dtxcms2,
                 samples = 1:5,
                 ID = NULL,
                 mz = 213.1869,
                 rt = 15.47,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 interactive = TRUE)


plotFeaturePeaks(obj = dtxcms2,
                 samples = NULL,
                 ID = NULL,
                 mz = 332.2200,
                 ppm = 20,
                 rtWindow = c(1000, 1200),
                 rtUnit = "sec",
                 rt = NULL,
                 interactive = TRUE)


plotFeaturePeaks(obj = dtxcms2,
                 samples = NULL,
                 ID = c("M263_R937_1829"),
                 mz = NULL,
                 rt = NULL,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 interactive = TRUE)


plotTargetCentroids(obj = dtxcms2,
                    samples = 4,
                    mz = 441.1670, ppm = 20,
                    rt = 15.27, rtWindow = 0.5,
                    rtUnit = "min", plotTargetMark = FALSE)




### Annotation ----------------------------------------------------------------------------------------------

dtxcms2@algorithms$annotation <- "alteredcamera"

dtxcms2@parameters$annotation <- AlteredCameraParam(
  sigma = 5,
  perfwhm = 0.3,
  cor_eic_th = 0.7,
  calcCaS = FALSE,
  cor_exp_th = 0.7,
  pval = 0.05,
  validateIsotopePatterns = TRUE,
  ppmIsotopes = 40,
  mzabs = 0.005, #not used, remove!
  noise = 300,
  searchAdducts = TRUE,
  ppmAdducts = 5,
  extendedList = FALSE
)

dtxcms2 <- annotateFeatures(dtxcms2, excludeBlanks = FALSE, save = FALSE)




#### Inspection ---------------------------------------------------------------------------------------------

# TODO make method for getting the full table
View(dtxcms2@annotation$comp)

#### Getter -----
#S4 method for getting components
View(components(dtxcms2,
                samples = NULL,
                ID = NULL,
                mz = 748.4842, ppm = 20,
                rt = 14.9, rtWindow = 2, rtUnit = "min",
                compNumber = NULL,
                entireComponents = TRUE,
                onlyAnnotated = TRUE,
                onlyRelated = TRUE))

diurond6 <- components(dtxcms2,
                       samples = NULL,
                       ID = NULL,
                       mz = 239.0628, ppm = 20,
                       rt = 15.62, rtWindow = 1, rtUnit = "min",
                       compNumber = NULL,
                       entireComponents = TRUE,
                       onlyAnnotated = TRUE,
                       onlyRelated = TRUE)

plotFeatures(obj = dtxcms2, ID = diurond6$ID, colorBy = "sampleGroups")


##### plot -----
plotComponents(obj = dtxcms2,
               samples = NULL,
               ID = NULL,
               mz = 748.4842, ppm = 20,
               rt = 14.9, rtWindow = 2, rtUnit = "min",
               comp = NULL,
               entireComponents = TRUE,
               onlyAnnotated = TRUE,
               onlyRelated = TRUE,
               log = FALSE,
               colorBy = "groups")

plotComponents(obj = dtxcms2,
               samples = c(1, 4),
               ID = NULL,
               mz = 239.0628, ppm = 20,
               rt = 15.62, rtWindow = 1, rtUnit = "min",
               comp = NULL,
               entireComponents = TRUE,
               onlyAnnotated = TRUE,
               onlyRelated = TRUE,
               log = FALSE,
               colorBy = "groups")




#### Consolidation of Annotation ----------------------------------------------------------------------------

View(dtxcms2@annotation$comp)

new("componentsCamera", xsa = dtxcms2@annotation$raw$Blank,
                        components = list(dtxcms2@annotation$raw$Blank),
                        componentInfo = data.table::data.table())




#### Annotation w/ patRoon -----

test <- patRoon::generateComponentsOpenMS(
  fGroups = dtxcms2@patdata,
  ionization = dtxcms2@polarity,
  chargeMin = 1,
  chargeMax = 3,
  chargeSpan = 3,
  qTry = "heuristic",
  potentialAdducts = patRoon::defaultOpenMSAdducts(dtxcms2@polarity),
  minRTOverlap = 0.7,
  retWindow = 1,
  absMzDev = 0.01,
  minSize = 2,
  relMinAdductAbundance = 1,
  adductConflictsUsePref = TRUE,
  NMConflicts = c("preferential"),
  prefAdducts = c("[M+H]+"),
  extraOpts = NULL
)

class(test)
View(patRoon::componentTable(test))
View(test)
testv1 <- test@featureComponents$`M200401068-r001`
testV1 <- rbindlist(testv1, idcol = "comp")
View(patRoon::findFGroup(test, "M239_R936_1314"))



test3 <- patRoon::generateComponentsCliqueMS(
  fGroups = dtxcms2@patdata,
  ionization = dtxcms2@polarity,
  maxCharge = 3,
  maxGrade = 5,
  ppm = 50,
  adductInfo = NULL,
  absMzDev = 0.01,
  minSize = 2,
  relMinAdductAbundance = 0.5,
  adductConflictsUsePref = FALSE,
  NMConflicts = c("mostIntense"),
  prefAdducts = c("[M+H]+"),
  extraOptsCli = NULL,
  extraOptsIso = NULL,
  extraOptsAnn = NULL,
  parallel = TRUE
)

class(test3)
View(patRoon::componentTable(test))
View(test3)

View(patRoon::componentTable(test3)[[patRoon::findFGroup(test3, "M239_R936_1314")]])

test3v1 <- test3@featureComponents$`M200401068-r001`
test3v1 <- rbindlist(test3v1, idcol = "comp")
View(test3v1)


comptest <- findFGroup(test3@featureComponents$`M200401068-r001`, "M441_R916_5186")
View(componentTable(test3)[[comptest]])
comptable <- rbindlist(test3@featureComponents$`M200401068-r001`, idcol = "CMP")
View(comptable)
comptable[comptable$group %in% "M239_R936_1314", ]
comptable[comptable$isogroup %in% 958, ]
View(comptable[comptable$name == "CMP386", ])


test3 <- cliqueMS::getCliques(mzdata = test3, filter = TRUE, mzerror = 5e-06, intdiff = 1e-04,
                              rtdiff = 1e-04, tol = 1e-05, silent = TRUE)
View(dtxcms3@workflows$SuspectScreening@results)

#### Build componentsCAMERA -----

test4 <- patRoon::generateComponentsCAMERA(
  fGroups = dtopenms2@patdata,
  ionization = dtopenms2@polarity,
  onlyIsotopes = FALSE,
  minSize = 2,
  relMinReplicates = 0.1,
  extraOpts = list(
    sigma = 5,
    perfwhm = 0.45,
    cor_eic_th = 0.85,
    graphMethod = "hcs",
    pval = 0.05,
    calcCiS = TRUE,
    calcIso = FALSE,
    calcCaS = FALSE, maxcharge = 2, maxiso = 5,
    ppm = 50, mzabs = 0.005, multiplier = 3, max_peaks = 100, intval = "into")
)

class(test4)

View(test4)


View(rbindlist(patRoon::componentTable(test4), idcol = "CMP"))




### QC check ------------------------------------------------------------------------------------------------
# Wrapping function for fast QC check using the workflow parameters.

QCListPath <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                     "/suspectList_MS2_pos.csv")

QCList <- getSuspectList(QCListPath)

dtxcms3 <- checkQC(obj = dtxcms2,
                   targets = QCList,
                   algorithmPeakPicking = NULL,
                   algorithmMakeFeatures = NULL,
                   algorithmAnnotation = NULL,
                   paramPeakPicking = NULL,
                   paramGrouping = NULL,
                   recursive = TRUE,
                   paramFill = NULL,
                   paramAnnotation = NULL,
                   rtWindow = 30,
                   ppm = 15,
                   MS2param = MS2param(),
                   exportResults = FALSE,
                   save = FALSE)

plotCheckQC(dtxcms3)

plotFeaturePeaks(obj = dtxcms2@QC$data,
                 ID = dtxcms2@QC$results$ID,
                 mz = NULL,
                 rt = NULL,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 names = dtxcms2@QC$results$name,
                 interactive = TRUE)

View(dtxcms3@QC$results)




### IS Check ------------------------------------------------------------------------------------------------

ISListPath <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                          "/suspectList_IS_pos.csv")

ISList <- getSuspectList(ISListPath)

dtxcms3 <- suspectScreening(obj = dtxcms3,
                            title = NULL,
                            samples = NULL,
                            suspectList = ISList,
                            ppm = 10,
                            rtWindow = 30,
                            adduct = NULL,
                            excludeBlanks = FALSE,
                            withMS2 = FALSE,
                            MS2param = MS2param())

View(dtxcms3@workflows$SuspectScreening@results)

## Add MS2 data to suspectList
ISList2 <- addMS2Info(wfobj = dtxcms3@workflows$SuspectScreening,
                      sampleGroup = "Blank",
                      suspectList = ISList,
                      MS2param = MS2param(),
                      updateRT = TRUE,
                      updateIntControl = TRUE,
                      save = "suspectList_IS_MS2_pos")

View(ISList2@data)

ISList2 <- getSuspectList(paste0(system.file(package = "ntsIUTA", dir = "extdata"), "/suspectList_IS_MS2_pos.csv"))

dtxcms3 <- checkIS(obj = dtxcms3,
                   targets = ISList2,
                   ppm = 10,
                   rtWindow = 30,
                   MS2param = MS2param(),
                   recoveryFrom = "Blank",
                   exportResults = FALSE,
                   save = FALSE)

View(dtxcms3@control$results$IN)

plotCheckIS(dtxcms3, rtWindow = 30, ppm = 10)




### Subset ntsData with features ----------------------------------------------------------------------------

# TODO update after final version of annotation workflow

test <- dtxcms2[1:2]

# filter without rebuilding mass and retention time of features
test2 <- filterFileFaster(dtxcms2, 1:2)



### Feature Metadata ---------------------------------------------------------------------------------------

# updates the features table with metadata, for improving filtering
meta <- patRoon::calculatePeakQualities(dtxcms3@patdata,
                                        weights = NULL,
                                        flatnessFactor = 0.05, #Passed to MetaClean as the flatness.factor argument to calculateJaggedness and calculateModality.
                                        avgFunc = max, # mean additional parameter for handling featureGroups
                                        parallel = TRUE)

#View(patRoon::groupQualities(meta))
#patRoon::groupScores(meta)



# TODO talk with Rick about applying a function for each indice based on the type max/min






### Filter features -----------------------------------------------------------------------------------------

dtxcms3 <- dtxcms2
# TODO make filter function/method and create filtering parameters slot
# TODO Adapt to ntsData, but integrate with features (main filter method)
filterPeaks(peaks, fileIndex = 4:5, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")
filterPeaks(peaksOpenms, fileIndex = 1:2, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")


### Workflow Parameters -------------------------------------------------------------------------------------

# TODO create framework with S4 classes and methods for assigning the workflow parameters to the ntsData object.

#### Algorithms -------------------

algorithmXCMS <- "xcms3"

algorithmOpenms <- "openms"

algorithmAnnotation <- "alteredcamera"




#### for Peak Picking -------------

paramXCMS <- xcms::CentWaveParam(
  ppm = 15, peakwidth = c(6, 60),
  snthresh = 5, prefilter = c(6, 300),
  mzCenterFun = "mean", integrate = 2,
  mzdiff = -0.0001, fitgauss = TRUE,
  noise = 0, verboseColumns = TRUE,
  firstBaselineCheck = FALSE,
  extendLengthMSW = TRUE
)

paramOpenms <- list(
  noiseThrInt = 250,
  chromSNR = 3,
  chromFWHM = 5,
  mzPPM = 15,
  reEstimateMTSD = TRUE,
  traceTermCriterion = "sample_rate",
  traceTermOutliers = 5,
  minSampleRate = 1,
  minTraceLength = 3,
  maxTraceLength = -1,
  widthFiltering = "fixed",
  minFWHM = 5,
  maxFWHM = 60,
  traceSNRFiltering = FALSE,
  localRTRange = 10,
  localMZRange = 6.5,
  isotopeFilteringModel = "metabolites (5% RMS)",
  MZScoring13C = FALSE,
  useSmoothedInts = TRUE,
  extraOpts = NULL,
  intSearchRTWindow = 3)




#### for Alignment and Grouping -----------------------------------------------------------------------------


##### Alignment ---------------------------------------------------------------------------------------------
#Not used if only one sample is given in peaks and when grouping is performed with openms

#Param preGrouping
#Only necessary if method for alignment is via PeaksGroups
param1 <- xcms::PeakDensityParam(
  sampleGroups = "holder",
  bw = 5,
  minFraction = 0.5,
  minSamples = 1,
  binSize = 0.008,
  maxFeatures = 100
)

#Parameters for alignment of retention time across samples
param2 <- xcms::PeakGroupsParam(
  minFraction = 1,
  extraPeaks = 0,
  smooth = "loess",
  span = 0.3,
  family = "gaussian"
)

paramAlignment <- list(param1, param2)


##### Grouping ----------------------------------------------------------------------------------------------
#Parameters for final grouping of peaks across samples

#With xcms3
paramGroupingxcms <- xcms::PeakDensityParam(
  sampleGroups = "holder",
  bw = 3,
  minFraction = 0.5,
  minSamples = 1,
  binSize = 0.008,
  maxFeatures = 100
)

#With openms
paramGroupingOpenms <- list(
  rtalign = TRUE,
  QT = FALSE,
  maxAlignRT = 5,
  maxAlignMZ = 0.008,
  maxGroupRT = 3,
  maxGroupMZ = 0.008,
  extraOptsRT = NULL,
  extraOptsGroup = NULL
)




#### Filling ------------------------------------------------------------------------------------------------
#Parameters for recursive integration

#Option 1
paramFill <- xcms::FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)
#Option 2
paramFill2 <- xcms::ChromPeakAreaParam(
  mzmin = function(z) quantile(z, probs = 0.25),
  mzmax = function(z) quantile(z, probs = 0.75),
  rtmin = function(z) quantile(z, probs = 0.25),
  rtmax = function(z) quantile(z, probs = 0.75)
)


### Suspect Screening WF ------------------------------------------------------------------------------------

getSuspectListTemplate(saveTo = system.file(package = "ntsIUTA", dir = "extdata"))

suspectListPath <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                          "/suspectList_MS2_pos.csv")

suspectList <- getSuspectList(suspectListPath)

dtxcms3 <- suspectScreening(obj = dtxcms2,
                            title = NULL,
                            samples = 4,
                            suspectList = suspectList,
                            ppm = 5,
                            rtWindow = 30,
                            adduct = NULL,
                            excludeBlanks = TRUE,
                            withMS2 = TRUE,
                            MS2param = MS2param())


View(dtxcms3@workflows$SuspectScreening@results)


# TODO prepare the history with info for rerun the code
# TODO add storage/drop method for processing steps
#dtprof@MSnExp@spectraProcessingQueue <- list()


