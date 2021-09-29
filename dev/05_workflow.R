

### Load ntsIUTA --------------------------------------------------------------------------------------------

library(ntsIUTA)

devtools::load_all()




### Project Setup -------------------------------------------------------------------------------------------

path <- system.file(package = "ntsIUTA", dir = "extdata")

dtall <- setupProject(path, polarity = "positive", save = FALSE, makeNewProject = FALSE)


#### Create New Session -----

setupProject(title = "Aopti_20200317",
             date = as.Date("2020-03-17"),
             polarity = "positive",
             save = TRUE, makeNewProject = TRUE)




### S4 Methods ----------------------------------------------------------------------------------------------


#### Sample Groups -------------------------------------------------------------

#Assign/Correct sample replicate groups
sampleGroups(dtall) <- c(rep("Blank", 3),
                         rep("IN", 3),
                         rep("OZ", 3),
                         rep("UV", 3),
                         rep("AC", 3),
                         rep("Centroid", 3),
                         rep("Profile", 3))

#getter for sample replicate group names
sampleGroups(dtall)

#getter for sample names (i.e. file names)
samples(dtall) #no setter




#### Blanks --------------------------------------------------------------------

#Assign the blank replicate sample group/s
blanks(dtall) <- "Blank"

#getter for blank replicate group names
blanks(dtall)




#### QC ------------------------------------------------------------------------

#When needed, assigning the QC sample replicate group
QC(dtall) <- "QC"

#getter for the QC sample names
QC(dtall)




#### show ntsData --------------------------------------------------------------

#show method
dtall


# TODO Add accessors to slots, such as title, path, etc.

# TODO show the addFiles function

# TODO show the convert files function




### Raw Data ------------------------------------------------------------------------------------------------

#sub-setting by samples
dtcent <- dtall[16:18]

dtcent <- importRawData(dtcent,
                        rtFilter = c(13.8, 16.3),
                        rtUnit = "min",
                        centroidedData = NA,
                        removeEmptySpectra = TRUE,
                        save = FALSE)




#### TIC, BPC and EIC ----------------------------------------------------------

exEIC <- extractEIC(dtcent,
                    samples = 1:2,
                    mz = 242.1434, ppm = 20,
                    rt = 14.8, rtWindow = 0.5,
                    rtUnit = "min")




#### Inspecting ----------------------------------------------------------------

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




#### Profile vs Centroid -------------------------------------------------------

dtprof <- dtall[19:21]

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




#### Centroiding Data ----------------------------------------------------------

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

#Data frame with carrelation summary
getCorReplicates(dt, exportResults = FALSE)


#### xcms3 ---------------------------------------------------------------------

dtxcms <- dt

dtxcms <- peakPickingParameters(dtxcms,
  algorithm = "xcms3",
  param = xcms::CentWaveParam(
    ppm = 15, peakwidth = c(6, 60),
    snthresh = 5, prefilter = c(6, 300),
    mzCenterFun = "mean", integrate = 2,
    mzdiff = -0.0001, fitgauss = TRUE,
    noise = 0, verboseColumns = TRUE,
    firstBaselineCheck = FALSE,
    extendLengthMSW = TRUE)
)

peakPickingParameters(dtxcms)

dtxcms <- peakPicking(obj = dtxcms, save = FALSE)




#### openms --------------------------------------------------------------------

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




#### sirius --------------------------------------------------------------------
#TODO Error with mz column

dtsirius <- dt

dtsirius@algorithms$peakPicking <- "sirius"

dtsirius@parameters$peakPicking <- list()

dtsirius <- peakPicking(obj = dtsirius, save = FALSE)




#### kpic2 ---------------------------------------------------------------------

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




#### safd ----------------------------------------------------------------------

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




#### Inspecting Peaks ----------------------------------------------------------

#Option 1
dtxcms@peaks[1, ]
dtopenms@peaks[1, ]

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


#### xcms3 ---------------------------------------------------------------------

dtxcms <- peakGroupingParameters(dtxcms,
  algorithm = "xcms3",
  param = list(
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
)

peakGroupingParameters(dtxcms)

dtxcms <- fillMissingParameters(dtxcms,
  algorithm = "xcms3",
  param = xcms::ChromPeakAreaParam(
    mzmin = function(z) quantile(z, probs = 0.25),
    mzmax = function(z) quantile(z, probs = 0.75),
    rtmin = function(z) quantile(z, probs = 0.25),
    rtmax = function(z) quantile(z, probs = 0.75)
  )
)

fillMissingParameters(dtxcms)

## second option for fillMissingParameters
xcms::FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)

dtxcms2 <- makeFeatures(dtxcms, save = FALSE)




#### openms --------------------------------------------------------------------

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




#### Ploting Alignment ---------------------------------------------------------
# TODO check plotAlignment, may only work when xcms3 is used as algorithm to align peaks

plotAlignment(features)




#### Inspecting Features -------------------------------------------------------


##### Getter -----

#Direct access to slot
dtxcms2@features[1, ]

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

dtxcms2 <- annotationParameters(dtxcms2,
  algorithm = "alteredcamera",
  param = AlteredCameraParam(
    sigma = 6,
    perfwhm = 0.4,
    cor_eic_th = 0.75,
    cor_exp_th = 0.75,
    pval = 0.1,
    validateIsotopePatterns = TRUE,
    ppmIsotopes = 40,
    noise = 300,
    searchAdducts = TRUE,
    ppmAdducts = 5,
    extendedList = FALSE
  )
)

annotationParameters(dtxcms2)

dtxcms3 <- annotateFeatures(dtxcms2, excludeBlanks = FALSE, save = FALSE)




#### Inspect -------------------------------------------------------------------

View(dtxcms3@annotation$comp)


##### Getter -----

#S4 method for getting components
View(components(dtxcms3,
                samples = NULL,
                ID = NULL,
                mz = 748.4842, ppm = 20,
                rt = 14.9, rtWindow = 2, rtUnit = "min",
                compNumber = NULL,
                entireComponents = TRUE,
                onlyAnnotated = TRUE,
                onlyRelated = TRUE))

diurond6 <- components(dtxcms3,
                       samples = NULL,
                       ID = NULL,
                       mz = 239.0628, ppm = 20,
                       rt = 15.62, rtWindow = 1, rtUnit = "min",
                       compNumber = NULL,
                       entireComponents = TRUE,
                       onlyAnnotated = TRUE,
                       onlyRelated = TRUE)

plotFeatures(obj = dtxcms3, ID = diurond6$ID, colorBy = "sampleGroups")




##### plot -----

plotComponents(obj = dtxcms3,
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

plotComponents(obj = dtxcms3,
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




#### Consolidation of Annotation -----------------------------------------------

View(dtxcms2@annotation$comp)

new("componentsCamera", xsa = dtxcms2@annotation$raw$Blank,
                        components = list(dtxcms2@annotation$raw$Blank),
                        componentInfo = data.table::data.table())




##### Annotation w/ patRoon -----

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




### Parameters ----------------------------------------------------------------------------------------------

dt <- dtall[1:6]

#### Peak Picking -----
dt <- peakPickingParameters(dt,
  algorithm = "xcms3",
  param = xcms::CentWaveParam(
    ppm = 15, peakwidth = c(6, 60),
    snthresh = 5, prefilter = c(6, 300),
    mzCenterFun = "mean", integrate = 2,
    mzdiff = -0.0001, fitgauss = TRUE,
    noise = 0, verboseColumns = TRUE,
    firstBaselineCheck = FALSE,
    extendLengthMSW = TRUE)
)

#### Peak Grouping -----
dt <- peakGroupingParameters(dt,
  algorithm = "xcms3",
  param = list(
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
)

#### Fill Missing -----
dt <- fillMissingParameters(dt,
  algorithm = "xcms3",
  param = xcms::ChromPeakAreaParam(
    mzmin = function(z) quantile(z, probs = 0.25),
    mzmax = function(z) quantile(z, probs = 0.75),
    rtmin = function(z) quantile(z, probs = 0.25),
    rtmax = function(z) quantile(z, probs = 0.75)
  )
)

#### Annotation -----
dt <- annotationParameters(dt,
  algorithm = "alteredcamera",
  param = AlteredCameraParam(
    sigma = 6,
    perfwhm = 0.4,
    cor_eic_th = 0.75,
    cor_exp_th = 0.75,
    pval = 0.1,
    validateIsotopePatterns = TRUE,
    ppmIsotopes = 40,
    noise = 300,
    searchAdducts = TRUE,
    ppmAdducts = 5,
    extendedList = FALSE
  )
)




### Wrapper PP, Grou & Anno ----------------------------------------------------------------------------------

dt <- createFeatures(dt, excludeBlanks = FALSE, save = FALSE)

dt




### QC check ------------------------------------------------------------------------------------------------
# Wrapping function for fast QC check using the workflow parameters.

QCListPath <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                     "/suspectList_MS2_pos.csv")

QCList <- getSuspectList(QCListPath)

dt1 <- checkQC(obj = dt,
                targets = QCList,
                rtWindow = NULL,
                ppm = NULL,
                exportResults = FALSE,
                save = FALSE)



#### Inspect QC ----------------------------------------------------------------

plotCheckQC(dt1)

plotFeaturePeaks(obj = qc2ntsData(dt1),
                 ID = dt1@QC@results$ID,
                 mz = NULL,
                 rt = NULL,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 names = dt1@QC@results$name,
                 interactive = TRUE)


# TODO make S4 method to get QC results
View(dt1@QC@results)




### IS Check ------------------------------------------------------------------------------------------------


#### Get MS2 of IS -------------------------------------------------------------

ISListPath <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                          "/suspectList_IS_pos.csv")

ISList <- getSuspectList(ISListPath)

dtxcms4 <- suspectScreening(obj = dtxcms3,
                            title = NULL,
                            samples = NULL,
                            suspectList = ISList,
                            ppm = 10,
                            rtWindow = 30,
                            adduct = NULL,
                            excludeBlanks = FALSE,
                            withMS2 = FALSE,
                            MS2param = MS2param())

View(dtxcms4@workflows$SuspectScreening@results)

## Add MS2 data to suspectList
ISList2 <- addMS2Info(wfobj = dtxcms4@workflows$SuspectScreening,
                      sampleGroup = "Blank",
                      suspectList = ISList,
                      MS2param = MS2param(),
                      updateRT = TRUE,
                      updateIntControl = TRUE,
                      save = "suspectList_IS_MS2_pos")

View(ISList2@data)



#### IS control ----------------------------------------------------------------

ISList2 <- getSuspectList(paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                                 "/suspectList_IS_MS2_pos.csv"))

dtxcms4 <- checkIS(obj = dtxcms3,
                   targets = ISList2,
                   ppm = NULL,
                   rtWindow = NULL,
                   MS2param = NULL,
                   recoveryFrom = "Blank",
                   exportResults = FALSE,
                   save = FALSE)

View(dtxcms4@IS@results$IN)

plotCheckIS(dtxcms4)

plotFeaturePeaks(obj = dtxcms4,
                 ID = dtxcms4@IS@results$Blank$ID,
                 mz = NULL,
                 rt = NULL,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 names = dtxcms4@IS@results$Blank$name,
                 interactive = TRUE)




### Subset ntsData  -----------------------------------------------------------------------------------------
# TODO update after final version of annotation workflow

# Subset by samples
test <- dtxcms2[1:2]

# filter without rebuilding mass and retention time of features
test2 <- filterFileFaster(dtxcms2, 1:2)

# subset features
test <- dtxcms3[, dtxcms3@control$results$IN$ID]




### Feature Metadata ---------------------------------------------------------------------------------------

dtxcms5 <- dtxcms4[4:6]

dtxcms5 <- calculateFeaturesMetadata(dtxcms5, ID = NULL)

ft <- dtxcms4@peaks

ft <- ft[ft$sample == "M200401068-r001", ]

ft <- ft[ft$intensity > 10000, ]

ft <- ft[c("ID", "mz", "rt", "intensity", "ZigZag", "ZigZagScore")]


IDs <- c(28378, 29841)


plotPeaks(obj = dtxcms4,
          ID = IDs[1],
          samples = NULL,
          mz = NULL,
          ppm = NULL,
          rt = NULL,
          rtUnit = "min",
          colorBy = "peaks")

plotPeaks(obj = dtxcms4,
          ID = IDs[2],
          samples = NULL,
          mz = NULL,
          ppm = NULL,
          rt = NULL,
          rtUnit = "min",
          colorBy = "peaks")

# TODO talk with Rick about applying a function for each indice based on the type max/min




### Filter features -----------------------------------------------------------------------------------------

dtxcms3 <- dtxcms2
# TODO add rt filter for the set, removing flushing and initial chromatographic time
# TODO make filter function/method and create filtering parameters slot
# TODO Adapt to ntsData, but integrate with features (main filter method)
filterPeaks(peaks, fileIndex = 4:5, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")
filterPeaks(peaksOpenms, fileIndex = 1:2, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")




### Suspect Screening WF ------------------------------------------------------------------------------------

getSuspectListTemplate(saveTo = system.file(package = "ntsIUTA", dir = "extdata"))

suspectListPath <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                          "/suspectList_MS2_pos.csv")

suspectList <- getSuspectList(suspectListPath)

dtxcms5 <- suspectScreening(obj = dtxcms4,
                            title = NULL,
                            samples = 4,
                            suspectList = suspectList,
                            ppm = 5,
                            rtWindow = 30,
                            adduct = NULL,
                            excludeBlanks = TRUE,
                            withMS2 = TRUE,
                            MS2param = MS2param())


View(dtxcms5@workflows$SuspectScreening@results)




### Find Fragments ------------------------------------------------------------------------------------------

dt2 <- qc2ntsData(dt1)

ID <- dt2@QC@results$ID

targets <- read.csv("F:/NTS_IUTA_Projects/art01_Aopti/tp_nitro.csv")


### Other ---------------------------------------------------------------------------------------------------




# TODO prepare the history with info for rerun the code
# TODO add storage/drop method for processing steps
#dtprof@MSnExp@spectraProcessingQueue <- list()
