
### Load ntsIUTA --------------------------------------------------------------------------------------------

library(ntsIUTA)




### Project Setup -------------------------------------------------------------------------------------------

path <- system.file(package = "ntsIUTA", dir = "extdata")

dtall <- setupProject(path, save = FALSE, makeNewProject = FALSE)




#### Class Methods ------------------------------------------------------------------------------------------

sampleGroups(dtall) <- c(rep("Blank", 3),
                         rep("IN", 3),
                         rep("OZ", 3),
                         rep("UV", 3),
                         rep("AC", 3),
                         rep("Centroid", 3),
                         rep("Profile", 3))

#see the samples (files)
samples(dtall) #no setter

#Assign the blank replicate sample group/s
blanks(dtall) <- "Blank"

#When needed, assigning the QC sample group after naming sample replicate groups
#Automatic detection of samples with QC in the name
#QC(dt) <- "QC"

#show method
dtall

#sub-setting by samples
dtcent <- dtall[16:18]
dtprof <- dtall[19:21]

# TODO show the addFiles function

# TODO show the convert files function


### Import Data ---------------------------------------------------------------------------------------------

dtcent <- importRawData(dtcent,
                        rtFilter = c(13.8, 16.3),
                        timeUnit = "min",
                        centroidedData = NA,
                        removeEmptySpectra = TRUE,
                        save = FALSE)

#get centroids
exEIC <- extractEIC(dtcent, fileIndex = NULL,
                    mz = 242.1434, ppm = 20,
                    rt = 14.8, rtWindow = 0.5,
                    rtUnit = "min")

#plot functions for TIC, BPC or EIC

#static
plotRawChrom(dtcent, fileIndex = NULL,
             type = "tic",
             mz = 242.1434, ppm = 20,
             rt = 14.8, rtWindow = 0.5,
             rtUnit = "min")

#iterative
plotlyRawChrom(dtcent, fileIndex = NULL,
               type = "tic",
               mz = 242.1434, ppm = 20,
               rt = 14.8, rtWindow = 0.5,
               rtUnit = "min")




### Profile vs Centroid -------------------------------------------------------------------------------------

dtprof <- importRawData(dtprof,
                        rtFilter = c(13.8, 16.3),
                        timeUnit = "min",
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

# TODO add storage/drop method for processing steps
#dtprof@MSnExp@spectraProcessingQueue <- list()

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
                    rtFilter = c(13, 17),
                    timeUnit = "min",
                    centroidedData = NA,
                    removeEmptySpectra = TRUE,
                    save = FALSE)




#### xcms3 --------------------------------------------------------------------------------------------------
dtxcms <- dt

param <- xcms::CentWaveParam(
  ppm = 15, peakwidth = c(6, 60),
  snthresh = 5, prefilter = c(6, 300),
  mzCenterFun = "mean", integrate = 2,
  mzdiff = -0.0001, fitgauss = TRUE,
  noise = 0, verboseColumns = TRUE,
  firstBaselineCheck = FALSE,
  extendLengthMSW = TRUE
)

dtxcms@parameters$peakPicking <- param

dtxcms <- peakPicking(obj = dtxcms, algorithm = "xcms3", save = FALSE)

dtxcms@patdata[1]




#### openms -------------------------------------------------------------------------------------------------
dtopenms <- dt

param2 <- list(
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

dtopenms@parameters$peakPicking <- param2

dtopenms <- peakPicking(obj = dtopenms, algorithm = "openms", save = FALSE)




#### Inspecting Peaks ---------------------------------------------------------------------------------------

dtxcms@peaks[[1]][1,]
dtopenms@peaks[[1]][1,]

# TODO Adapt to ntsData, but integrate with features
filterPeaks(peaks, fileIndex = 4:5, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")
filterPeaks(peaksOpenms, fileIndex = 1:2, mz = 748.4842, ppm = 10, rt = 14.9, rtUnit = "min")




### Alignment and Grouping ----------------------------------------------------------------------------------


#### Param Alignment ----------------------------------------------------------------------------------------
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




#### Param Grouping -----------------------------------------------------------------------------------------
#Parameters for final grouping of peaks across samples

#With xcms
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
  QT = FALSE,
  maxAlignRT = 5,
  maxAlignMZ = 0.008,
  maxGroupRT = 3,
  maxGroupMZ = 0.008,
  extraOptsRT = NULL,
  extraOptsGroup = NULL
)




#### Param Filling ------------------------------------------------------------------------------------------
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




#### xcms3 --------------------------------------------------------------------------------------------------
dtxcms <- makeFeatures(obj = dtxcms,
                       algorithm = "xcms3",
                       rtAlignment = TRUE,
                       paramAlignment = paramAlignment,
                       paramGrouping = paramGrouping,
                       recursive = TRUE,
                       paramFill = paramFill,
                       save = FALSE)




#### openms -------------------------------------------------------------------------------------------------
dtopenms <- makeFeatures(obj = dtopenms,
                         algorithm = "openms",
                         rtAlignment = TRUE,
                         paramGrouping = paramGroupingOpenms,
                         recursive = FALSE,
                         paramFill = NULL,
                         save = FALSE)

# TODO improve buildFeatureList function


#### Ploting Alignment --------------------------------------------------------------------------------------
# TODO check alignment, may only work when xcms3 is used as algorithm to align peaks
plotAlignment(features)




#### Inspecting Features ------------------------------------------------------------------------------------

patRoon::plotChroms(features[, 1:2], col = "none")

# TODO static and iterative plot for peaks and features
plotFeatures(features, features = "M201_R784_1")


