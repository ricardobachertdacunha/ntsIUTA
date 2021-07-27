
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

exEIC <- extractEIC(dtcent, fileIndex = NULL,
                    mz = 242.1434, ppm = 20,
                    rt = 14.8, rtWindow = 0.5,
                    rtUnit = "min")




#### Inspecting Raw Data ------------------------------------------------------------------------------------

plotRawChrom(dtcent, fileIndex = NULL,
             type = "tic",
             mz = 242.1434, ppm = 20,
             rt = 14.8, rtWindow = 0.5,
             colorBy = "samples",
             rtUnit = "min",
             interactive = FALSE)

#iterative
plotRawChrom(dtcent, fileIndex = NULL,
             type = "tic",
             mz = 242.1434, ppm = 20,
             rt = 14.8, rtWindow = 0.5,
             rtUnit = "min",
             colorBy = "samples",
             interactive = TRUE)

#plot correlation of replicate sample groups
plotCorReplicates(dtcent)




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
p1 <- plotTargetCentroids(dtcent, fileIndex = 1,
                          mz = 242.1434, ppm = 150,
                          rt = 14.8, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

p2 <- plotTargetCentroids(dtprof, fileIndex = 1,
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


p1 <- plotTargetCentroids(dtcent, fileIndex = 1,
                          mz = 242.1434, ppm = 30,
                          rt = 14.75, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

p2 <- plotTargetCentroids(dtprof, fileIndex = 1,
                          mz = 242.1434, ppm = 30,
                          rt = 14.75, rtWindow = 0.8,
                          rtUnit = "min", plotTargetMark = TRUE)

plotly::subplot(list(p1, p2), nrows = 1, margin = 0.04)




### Peak Picking --------------------------------------------------------------------------------------------

dt <- dtall[1:6]

dt <- importRawData(dt,
                    rtFilter = NULL,
                    rtUnit = "min",
                    centroidedData = NA,
                    removeEmptySpectra = TRUE,
                    save = FALSE)

# TODO when applying rtFilter for import raw data reuse when peak picking with patRoon


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

# TODO add importRawData during peak picking when the OndiskMsnExp is not present


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




#### Inspecting Peaks ---------------------------------------------------------------------------------------

#Option 1
dtxcms@peaks[1,]
dtopenms@peaks[1,]

#Option 2, using the S4 method
peaks(dtxcms, fileIndex = NULL, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")





### Alignment and Grouping ----------------------------------------------------------------------------------


#### xcms3 --------------------------------------------------------------------------------------------------

dtxcms@algorithms$makeFeatures <- "xcms3"

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

param2 <- xcms::PeakGroupsParam(
  minFraction = 1,
  extraPeaks = 0,
  smooth = "loess",
  span = 0.3,
  family = "gaussian"
)

paramAlignment <- list(param1, param2)

dtxcms@parameters$peakAlignment <- paramAlignment

paramGroupingxcms <- xcms::PeakDensityParam(
  sampleGroups = "holder",
  bw = 3,
  minFraction = 0.5,
  minSamples = 1,
  binSize = 0.008,
  maxFeatures = 100
)

dtxcms@parameters$peakGrouping <- paramGroupingxcms

paramFill <- xcms::FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)

dtxcms@parameters$fillMissing <- paramFill

dtxcms2 <- makeFeatures(obj = dtxcms,
                       algorithm = NULL,
                       rtAlignment = TRUE,
                       paramAlignment = NULL,
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
features(dtxcms2, fileIndex = NULL, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")


##### plot Chroms -----

#Interactive feature/s plot
plotFeatures(obj = dtxcms2,
             fileIndex = 3:4,
             ID = c("M275_R622_2158" , "M210_R586_836"),
             mz = 213.1869,
             rt = 15.47,
             rtUnit = "min",
             ppm = NULL,
             rtWindow = NULL,
             colorBy = "features",
             interactive = TRUE)

#Static feature/s plot
plotFeatures(obj = dtxcms2,
             fileIndex = NULL,
             ID = c("M275_R622_2158" , "M210_R586_836"),
             mz = 213.1869,
             rt = 15.47,
             rtUnit = "min",
             ppm = NULL,
             rtWindow = NULL,
             colorBy = "features",
             interactive = FALSE)

#patRoon option
patRoon::plotChroms(dtxcms2@patdata[, c("M275_R622_2158" , "M210_R586_836")])


##### plot Feature Peaks -----

plotFeaturePeaks(obj = dtxcms2,
                 fileIndex = 3:4,
                 ID = c("M275_R622_2158" , "M210_R586_836"),
                 mz = 213.1869,
                 rt = 15.47,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 interactive = TRUE)




### QC check ------------------------------------------------------------------------------------------------
# Wrapping function for fast QC check using the workflow parameters.

targets <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                  "/QC_ScreeningList_ntsIUTA_MS2_pos.csv")

dtxcms2 <- checkQC(obj = dtxcms2,
  targets = targets,
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
  save = FALSE)

plotCheckQC(dtxcms2)

plotFeaturePeaks(obj = dtxcms2@QC$data,
                 ID = dtxcms2@QC$results$ID,
                 mz = NULL,
                 rt = NULL,
                 rtUnit = "min",
                 ppm = NULL,
                 rtWindow = NULL,
                 names = dtxcms2@QC$results$name,
                 interactive = TRUE)

### Annotation ----------------------------------------------------------------------------------------------

dtxcms2@algorithms$annotation <- "alteredcamera"

dtxcms2@parameters$annotation <- AlteredCameraParam()

dtxcms2 <- annotateFeatures(dtxcms2, save = FALSE)


#### Inspection ---------------------------------------------------------------------------------------------

#### Getter -----
#S4 method for getting components
components(dtxcms2,
           samples = 4:5,
           ID = NULL,
           mz = 748.4842, ppm = 20,
           rt = 14.9, rtWindow = 2, rtUnit = "min",
           compNumber = NULL,
           entireComponents = TRUE,
           onlyAnnotated = TRUE,
           onlyRelated = TRUE)

##### plot -----
plotComponents(obj = dtxcms2,
               samples = 4:5,
               ID = NULL,
               mz = 748.4842, ppm = 20,
               rt = 14.9, rtWindow = 2, rtUnit = "min",
               comp = NULL,
               entireComponents = TRUE,
               onlyAnnotated = TRUE,
               onlyRelated = TRUE,
               log = FALSE,
               colorBy = "groups")


# TODO prepare the history with info for rerun the code
# TODO add storage/drop method for processing steps
#dtprof@MSnExp@spectraProcessingQueue <- list()

# TODO Check the function getScreeningListTemplate() and alter in method that use screening list


### Filter features -----------------------------------------------------------------------------------------

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

suspects <- paste0(system.file(package = "ntsIUTA", dir = "extdata"),
                   "/QC_ScreeningList_ntsIUTA_MS2_pos.csv")

dtxcms2 <- screenSuspectsFromCSV(dtxcms2,
                                 name = NULL,
                                 suspects = suspects,
                                 ppm = 5,
                                 rtWindow = 30,
                                 adduct = NULL,
                                 excludeBlanks = TRUE,
                                 withMS2 = TRUE,
                                 MS2param = MS2param())
