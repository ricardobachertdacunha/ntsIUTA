

### Load ntsIUTA --------------------------------------------------------------------------------------------

library(ntsIUTA)

devtools::load_all()




### Project Setup -------------------------------------------------------------------------------------------

path <- "C:\\Users\\Ricardo\\Documents\\R_Demo Project"

object <- setupProject(
  path = path,
  title = "Demo Project",
  description = "Demonstration of ntsIUTA.",
  date = Sys.Date(),
  convertFiles = FALSE,
  convertFrom = NULL,
  convertToCentroid = TRUE,
  replicates = NULL,
  polarity = "positive",
  method = NA_character_,
  save = TRUE,
  makeNewProject = FALSE
)


#### addFiles ------------------------------------------------------------

## method to add files after project setup. Files from other locations
#then the project path can be added via the addFiles.
object <- addFiles(
  files = utils::choose.files(),
  object = object,
  copy = FALSE,
  replicates = NULL,
  polarity = "positive",
  method = NULL
)


#### mzMLconverter -------------------------------------------------------

## function to convert vendor files to mzML with or without centroiding.
mzMLconverter(
  path = path(object),
  files = NULL,
  convertFrom = "agilent",
  centroidMethod = "vendor",
  outPath = NULL,
  overwrite = TRUE
)




### ntsData methods and functions ----------------------------------------------------------------------------


#### path ----------------------------------------------------------------

## getter for the project path
path(object)


#### projectInfo ---------------------------------------------------------

## setter for the projectInfo
object <- projectInfo(
  object,
  title = "Demo Project Changed Title",
  description = "Another description example.",
  date = NULL
)

#NOTE: when an argumment is missing and/or NULL, such as the date, is not updated/changed.

## getter for the project information
projectInfo(object)

#NOTE: When projectInfo is called without arguments
#it return a list with the title, description and date of the project.


#### samplesTable --------------------------------------------------------

## getter for the samples table
samplesTable(object)


#### filePaths -----------------------------------------------------------

## getter for the file paths
filePaths(object)


#### samples -------------------------------------------------------------

## getter for sample names (i.e. file names)
samples(object)

#NOTE: there is no setter for samples method, as samples reflects the file name


#### replicates ----------------------------------------------------------

## setter for sample replicate names
replicates(object) <- c(
  rep("Blank", 3),
  rep("IN", 3),
  rep("OZ", 3),
  rep("UV", 3),
  rep("AC", 3),
  rep("QC", 3),
  rep("Centroid", 3),
  rep("Profile", 3),
  "Centroid_mzMLconverter"
)

## getter for sample replicate names
replicates(object)


#### blanks --------------------------------------------------------------

## setter for the blank sample replicate/s
blanks(object) <- "Blank"

## getter for the blank sample replicate/s
blanks(object)


#### polarity ------------------------------------------------------------

## setter for the polarity mode of each sample/file
polarity(object) <- "positive"


## getter for the polarity mode of each sample or replicate
polarity(object)
#or
polarity(object, groupby = "replicates")


#### acquisitionMethods --------------------------------------------------

## setter for acquisition method names
acquisitionMethods(object) <- "NTS_MethodNameExample"

#NOTE: For different method in the set, the names should be given per sample.

## getter for acquisition method names
acquisitionMethods(object)


#### metadata ------------------------------------------------------------

## function to add metadata in an ntsData object
object <- addMetadata(
  object,
  var = c(rep("WW", 5), "QC", rep("Dev", 3)),
  varname = "datatype"
)

## getter for metadata
metadata(object, varname = "datatype")

## remove metadata
object <- removeMetadata(
  object,
  varname = "datatype"
)


#### QC ------------------------------------------------------------------

## setter of QC sample replicate/s
QC(object) <- "QC"

#NOTE: Sample name/s can be given instead of replicates
#by adding argument nametype = "samples".

## getter for the QC samples date table
QC(object)

## method for restauring QC samples/replicates back to the samples slot
QC(object, remove = TRUE) <- "QC"

#NOTE: The same argumment nametype = "samples" can be applied
#for moving samples from the QC slot to the samples slot.


#### show ntsData --------------------------------------------------------

## method to show details from the ntsData object
object


#### sub-setting simple [#] ----------------------------------------------

object[1:3]




### Inspect Raw Data ----------------------------------------------------------------------------------------

## sub-setting by samples
example01 <- object[16:18]


#### Raw Info ------------------------------------------------------------

## getter for the data table with raw information for all or specified files
getRawInfo(example01, samples = 1)

## adds the raw information directly into the ntsData object
example01 <- addRawInfo(example01)
samplesTable(example01)


#### TICs -----------------------------------------------------------------

## samples can be defined by index or name
tic_ex <- TICs(example01, samples = NULL)

##### visualization ---------

plotTICs(example01, interactive = TRUE)

plotTICs(tic_ex)


#### EICs -----------------------------------------------------------------

## extract EIC based on target mz and rt pair with fixed mass and time deviations
mz_01 <- c(247.1651, 239.0628)
rt_01 <- c(839, 937)

eic_ex1 <- EICs(example01, samples = NULL, mz = mz_01, rt = rt_01, ppm = 10, sec = 30)

## extract EIC based on minimum and maximum m\z and retention time
mz_02 <- data.frame(mzmin = c(247.1626, 239.0604), mzmax = c(247.1676, 239.0652))
rt_02 <- data.frame(rtmin = c(809, 907), rtmax = c(869, 967))

eic_ex2 <- EICs(example01, samples = NULL, mz = mz_02, rt = rt_02, ppm = 10, sec = 30)
#Note: ppm and sec are ignored as both mass and time ranges are given by mz and rt argumments

mz_03 <- data.frame(
  id = c("target1", "target2"),
  mz = c(247.1651, 239.0628),
  rt = c(839, 937)
)
eic_ex3 <- EICs(example01, samples = NULL, mz = mz_03, rt = NULL, ppm = 10, sec = 30)
#Note: rt taken from mz and ppm and sec are used to calculate deviations

mz_04 <- data.frame(
  id = c("target1", "target2"),
  mzmin = c(247.1626, 239.0604), mzmax = c(247.1676, 239.0652),
  rtmin = c(809, 907), rtmax = c(869, 967)
)
eic_ex4 <- EICs(example01, samples = NULL, mz = mz_04, rt = NULL, ppm = 10, sec = 30)
#Note: the ppm and sec deviations are not used, as ranges are already given


##### visualization ---------

#static
plotEICs(
  example01,
  samples = NULL,
  mz = mz_04, ppm = 10,
  rt = NULL, sec = 30,
  colorBy = "targets",
  legendNames = NULL,
  interactive = FALSE
)

#static, with NULL rt which gets the full time range
plotEICs(
  example01,
  samples = NULL,
  mz = mz_01, ppm = 10,
  rt = NULL, sec = 30,
  colorBy = "samples",
  legendNames = NULL,
  interactive = FALSE
)

#iterative
plotEICs(
  example01,
  samples = 1:2,
  mz = mz_01, ppm = 10,
  rt = rt_01, sec = 30,
  colorBy = "targets",
  legendNames = NULL,
  interactive = TRUE
)

## Note: The colorBy argument can also be set to targets to color
# against each EIC or a legendNames character vector with the same length as
# the number of EICs (i.e., mz/rt pairs) can be given to use as legend.
# TODO read about replayPlot(obj) not needed

#with a data table from the EICs method
plotEICs(
  eic_ex1,
  samples = NULL,
  targets = unique(eic_ex1$id),
  colorBy = "samples",
  legendNames = NULL,
  interactive = FALSE
)

## Note: Ploting from a existing EIC data table, the mz and rt do not work
# as when object is an ntsData object. Instead, the targets argument can be defined with the id
# of targets for plotting using the id column of the EIC data table.


#### XICs -----------------------------------------------------------------

xic_ex1 <- XICs(example01, samples = NULL, mz = mz_01, rt = rt_01, ppm = 10, sec = 30)

xic_ex2 <- XICs(example01, samples = NULL, mz = mz_02, rt = rt_02, ppm = 10, sec = 30)

##### Visualisation ---------

## with defined mz and rt pairs and standard deviations
plotXICs(
  object,
  samples = 3:4,
  mz = mz_01, ppm = 20,
  rt = rt_01, sec = 30,
  legendNames = NULL,
  plotTargetMark = TRUE,
  targetsMark = NULL,
  ppmMark = 5,
  secMark = 10,
  numberRows = 2
)

## with defined mz and rt spaces using a data.frame/data.table
# Note that targets are not plotted as mz and rt are NULL
plotXICs(
  object,
  samples = 3:4,
  mz = mz_02, ppm = NULL,
  rt = rt_02, sec = NULL,
  legendNames = c("Target01", "Target02"),
  plotTargetMark = TRUE,
  targetsMark = data.frame(mz = mz_01, rt = rt_01),
  ppmMark = 5,
  secMark = 10,
  numberRows = 2
)

## via a data.table produced using XICs()
# Note that target mark is build as mean of the given mz_id and rt_id ranges in xic_ex2
plotXICs(
  xic_ex2,
  samples = 1:2,
  targets = NULL,
  legendNames = c("Target01", "Target02"),
  plotTargetMark = TRUE,
  targetsMark = NULL,
  ppmMark = 5,
  secMark = 10,
  numberRows = 2
)


#### MS2 -----------------------------------------------------------------

# same functionality as EICs for setting mz/rt pairs
# below the example for a combined table of rt and mz predefined deviations for two targets
ms2_ex01 <- MS2s(
  object = object,
  samples = 3:4,
  mz = mz_04, ppm = 20,
  rt = NULL, sec = 60,
  ppmClustering = 50,
  minIntensityPre = 250,
  minIntensityPost = 250,
  mergeCEs = TRUE,
  mergeBy = NULL
)

##### visualization ---------

# ploting through an ntsData object
plotMS2s(
  object = object,
  samples = 3:4,
  mz = mz_04, ppm = 20,
  rt = NULL, sec = 60,
  clusteringUnit = "ppm",
  clusteringMethod = "distance",
  clusteringUnit = "ppm",
  clusteringWindow = 15,
  minIntensityPre = 250,
  minIntensityPost = 250,
  mergeCEs = TRUE,
  mergeBy = NULL,
  colorBy = "targets"
)

# plots already produce table from MS2s
plotMS2s(
  object = ms2_ex01,
  colorBy = "targets",
  interactive = TRUE
)

# when merging by replicates and ploting the same target for two replicates
plotMS2s(
  object = object,
  samples = 3:4,
  mz = mz_04[1, ], ppm = 20,
  rt = NULL, sec = 60,
  clusteringMethod = "distance",
  clusteringUnit = "ppm",
  clusteringWindow = 15,
  minIntensityPre = 100,
  minIntensityPost = 150,
  mergeCEs = TRUE,
  mergeBy = "replicates",
  colorBy = "replicates",
  interactive = TRUE
)








































### Basic Workflow -------------------------------------------------------------------------------------------

dt <- object[1:6]


#### add metadata --------------------------------------------------------

var <- data.table::data.table(matrix = c("clean", "dirty"))
dt <- addMetadata(dt, var)

metadata(dt)


#### check TICs ----------------------------------------------------------

tic <- TICs(dt, samples = NULL)

plotTICs(tic, samples = NULL, colorBy = "replicates")



#### replicates quality --------------------------------------------------

# TODO Make an alternative using mzR
plotCorReplicates(dt, binsize = 3)
getCorReplicates(dt, exportResults = FALSE)




#### Peak Picking --------------------------------------------------------------------------------------------

##### xcms3 ---------------------------------------------------------------------

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
    extendLengthMSW = TRUE
  )
)

peakPickingParameters(dtxcms)

dtxcms <- peakPicking(obj = dtxcms, save = FALSE)




##### openms --------------------------------------------------------------------

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
  intSearchRTWindow = 3
)

dtopenms@parameters$peakPicking <- paramOpenms

dtopenms <- peakPicking(obj = dtopenms, save = FALSE)




##### sirius --------------------------------------------------------------------
#TODO Error with mz column

dtsirius <- dt

dtsirius@algorithms$peakPicking <- "sirius"

dtsirius@parameters$peakPicking <- list()

dtsirius <- peakPicking(obj = dtsirius, save = FALSE)




##### kpic2 ---------------------------------------------------------------------

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




##### safd ----------------------------------------------------------------------

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




##### Inspecting Peaks ----------------------------------------------------------

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




#### Alignment and Grouping ----------------------------------------------------------------------------------


##### xcms3 ---------------------------------------------------------------------

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




##### openms --------------------------------------------------------------------

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
  paramFill = NULL,
  save = FALSE
)




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
  interactive = TRUE
)

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
  interactive = FALSE
)

plotFeatures(obj = dtxcms2,
  samples = 3:4,
  ID = NULL,
  mz = 332.2200,
  rt = NULL,
  rtUnit = "sec",
  ppm = 20,
  rtWindow = c(900, 1200),
  colorBy = "samples",
  interactive = FALSE
)

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
  interactive = TRUE
)


plotFeaturePeaks(obj = dtxcms2,
  samples = NULL,
  ID = NULL,
  mz = 332.2200,
  ppm = 20,
  rtWindow = c(1000, 1200),
  rtUnit = "sec",
  rt = NULL,
  interactive = TRUE
)


plotFeaturePeaks(obj = dtxcms2,
  samples = NULL,
  ID = c("M263_R937_1829"),
  mz = NULL,
  rt = NULL,
  rtUnit = "min",
  ppm = NULL,
  rtWindow = NULL,
  interactive = TRUE
)


plotTargetCentroids(obj = dtxcms2,
  samples = 4,
  mz = 441.1670, ppm = 20,
  rt = 15.27, rtWindow = 0.5,
  rtUnit = "min", plotTargetMark = FALSE
)




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
    onlyRelated = TRUE
  )
)

diurond6 <- components(dtxcms3,
  samples = NULL,
  ID = NULL,
  mz = 239.0628, ppm = 20,
  rt = 15.62, rtWindow = 1, rtUnit = "min",
  compNumber = NULL,
  entireComponents = TRUE,
  onlyAnnotated = TRUE,
  onlyRelated = TRUE
)

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
  colorBy = "groups"
)

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
  colorBy = "groups"
)




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
    extendLengthMSW = TRUE
  )
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
      maxFeatures = 100
    ),
   preGroupParam = xcms::PeakDensityParam(
      sampleGroups = "holder",
      bw = 5,
      minFraction = 0.5,
      minSamples = 1,
      binSize = 0.008,
      maxFeatures = 100
    ),
   retAlignParam = xcms::PeakGroupsParam(
      minFraction = 1,
      extraPeaks = 0,
      smooth = "loess",
      span = 0.3,
      family = "gaussian"
    )
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

QCListPath <- paste0(
  system.file(package = "ntsIUTA", dir = "extdata"),
  "/suspectList_MS2_pos.csv"
)

QCList <- getSuspectList(QCListPath)

dt1 <- checkQC(obj = dt,
  targets = QCList,
  rtWindow = NULL,
  ppm = NULL,
  exportResults = FALSE,
  save = FALSE
)



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
  interactive = TRUE
)


# TODO make S4 method to get QC results
View(dt1@QC@results)




### IS Check ------------------------------------------------------------------------------------------------


#### Get MS2 of IS -------------------------------------------------------------

ISListPath <- paste0(
  system.file(package = "ntsIUTA", dir = "extdata"),
  "/suspectList_IS_pos.csv"
)

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
  MS2param = MS2param()
)

View(dtxcms4@workflows$SuspectScreening@results)

## Add MS2 data to suspectList
ISList2 <- addMS2Info(wfobj = dtxcms4@workflows$SuspectScreening,
  sampleGroup = "Blank",
  suspectList = ISList,
  MS2param = MS2param(),
  updateRT = TRUE,
  updateIntControl = TRUE,
  save = "suspectList_IS_MS2_pos"
)

View(ISList2@data)



#### IS control ----------------------------------------------------------------

ISList2 <- getSuspectList(
  paste0(system.file(package = "ntsIUTA", dir = "extdata"),
    "/suspectList_IS_MS2_pos.csv"
  )
)

dtxcms4 <- checkIS(obj = dtxcms3,
  targets = ISList2,
  ppm = NULL,
  rtWindow = NULL,
  MS2param = NULL,
  recoveryFrom = "Blank",
  exportResults = FALSE,
  save = FALSE
)

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
  interactive = TRUE
)




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

dtxcms5 <- dtxcms5[, dtxcms5@IS@results$IN$ID]

dtxcms5 <- calculateSNR(dtxcms5)

dtxcms5 <- calculateFeaturesMetadata(dtxcms5)



ft <- dtxcms5@features

pk <- dtxcms5@peaks

View(ft)

View(pk)

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
  colorBy = "peaks"
)

plotPeaks(obj = dtxcms4,
  ID = IDs[2],
  samples = NULL,
  mz = NULL,
  ppm = NULL,
  rt = NULL,
  rtUnit = "min",
  colorBy = "peaks"
)

# TODO talk with Rick about applying a function for each indice based on the type max/min




### Filter features -----------------------------------------------------------------------------------------

dtxcms6 <- dtxcms4

dtxcms6 <- filterFeatures(obj = dtxcms6, filterMinInt = 5000)

dtxcms7 <- removeFilteredFeatures(dtxcms6)

dtxcms7@features$sn <- NA

dtxcms8 <- filterFeatures(obj = dtxcms7, filterBlank = 3)

dtxcms8 <- removeFilteredFeatures(dtxcms8)

dtxcms9 <- restoreFilteredFeatures(dtxcms8)

all.equal(dtxcms8@features$ID, dtxcms6@features$ID)

View(dtxcms7@removed)

View(dtxcms8@features)




### Workflows ------------------------------------------------------------------------------------------------


#### Suspect Screening WF ------------------------------------------------------------------------------------

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




#### Find Fragments ------------------------------------------------------------------------------------------

dt2 <- qc2ntsData(dt1)

targets <- read.csv("F:/NTS_IUTA_Projects/art01_Aopti/tp_nitro.csv")


dt2 <- findFragments(dt2,
                     targets = targets,
                     replicates = NULL,
                     ID = dt2@QC@results$ID,
                     ppm = 10,
                     intMin = 10)

class(dt2@workflows$MS2FragmentsScreening)

View(dt2@workflows$MS2FragmentsScreening@results)




#### Monitoring ---------------------------------------------------------------------------------------------

dt3 <- dtall[1:15]

dt3@parameters <- dt2@parameters

dt3 <- checkQC(dt3, targets = QCList,  save = FALSE)

plotCheckQC(dt3)

dt3 <- createFeatures(dt3)

dt3 <- checkIS(dt3, targets = ISList, save = FALSE)

plotCheckIS(dt3)

dt3 <- filterFeatures(obj = dt3, filterMinInt = 3000, filterBlank = 2)

dt3 <- removeFilteredFeatures(dt3)

dt3

dt3 <- processMonitoring(dt3, sequences = list(Aopti = c("IN", "OZ", "UV", "AC")),
                         title = NULL,
                         constantLevel = 5000,
                         doClustering = TRUE,
                         minClust = 10,
                         sizeClus = 2:100,
                         exportResults = FALSE)

object <- dt3@workflows$processMonitoring

catplot <- plotProcessCategories(dt3@workflows$processMonitoring, sequences = "Aopti", yUnit = "rt")

effplot <- plotProcessEfficiency(dt3@workflows$processMonitoring, sequences = "Aopti")

clusplot <- plotProcessClusters(dt3@workflows$processMonitoring, sequences = "Aopti")

exportPlots(dt3@workflows$processMonitoring, path = "C:/Users/MZmine/Desktop", interactive = TRUE)



### Other ---------------------------------------------------------------------------------------------------

#### Profile vs Centroid -------------------------------------------------------

dtprof <- dtall[19:21]

dtprof <- importRawData(dtprof,
  rtFilter = c(13.8, 16.3),
  rtUnit = "min",
  centroidedData = FALSE,
  removeEmptySpectra = TRUE,
  save = FALSE
)

# TODO Is centroided data working or not working
#check using raw data without centroids
table(MSnbase::isCentroidedFromFile(dtcent@MSnExp), MSnbase::msLevel(dtcent@MSnExp))
table(MSnbase::isCentroidedFromFile(dtprof@MSnExp), MSnbase::msLevel(dtprof@MSnExp))

# plot centroided and profile data
p1 <- plotTargetCentroids(dtcent, samples = 1,
  mz = 242.1434, ppm = 150,
  rt = 14.8, rtWindow = 0.8,
  rtUnit = "min", plotTargetMark = TRUE
)

p2 <- plotTargetCentroids(dtprof, samples = 1,
  mz = 242.1434, ppm = 150,
  rt = 14.8, rtWindow = 0.8,
  rtUnit = "min", plotTargetMark = TRUE
)

plotly::subplot(list(p1, p2), nrows = 1, margin = 0.04)




#### Centroiding Data ----------------------------------------------------------

dtprof <- centroidProfileData(obj = dtprof,
  halfwindow = 3,
  SNR = 3,
  noiseMethod = "MAD",
  methodRefineMz = "kNeighbors",
  k = 1,
  smoothing = FALSE,
  save = FALSE
)

p1 <- plotTargetCentroids(dtcent, samples = 1,
  mz = 242.1434, ppm = 30,
  rt = 14.75, rtWindow = 0.8,
  rtUnit = "min", plotTargetMark = TRUE
)

p2 <- plotTargetCentroids(dtprof, samples = 1,
  mz = 242.1434, ppm = 30,
  rt = 14.75, rtWindow = 0.8,
  rtUnit = "min", plotTargetMark = TRUE
)

plotly::subplot(list(p1, p2), nrows = 1, margin = 0.04)


#### Real FullData -----

# object <- setupProject(title = "Test Project", date = Sys.Date(),
#   polarity = "positive", save = FALSE, makeNewProject = FALSE
# )



# TODO prepare the history with info for rerun the code
# TODO add storage/drop method for processing steps
#dtprof@MSnExp@spectraProcessingQueue <- list()
