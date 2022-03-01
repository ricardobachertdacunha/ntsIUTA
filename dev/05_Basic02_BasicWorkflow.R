

### setup ---------------------------------------------------------------------------------------------------

library(ntsIUTA)
devtools::load_all()

path <- "C:\\Users\\Ricardo\\Documents\\R_Demo Project"

dt <- setupProject(
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

dt <- dt[c(1:6, 16:18)]

replicates(dt) <- c(
  rep("Blank", 3),
  rep("IN", 3),
  rep("QC", 3)
)

blanks(dt) <- rep("Blank", 9)

samplesTable(dt)




### add metadata --------------------------------------------------------------------------------------------

var <- data.table::data.table(matrix = c(rep("clean", 3), rep("dirty", 3), rep("control", 3)))
dt <- addMetadata(dt, var)

metadata(dt)




### check TICs ----------------------------------------------------------------------------------------------

tic <- TICs(dt, samples = NULL)

plotTICs(tic, samples = NULL, colorBy = "replicates")




### peak picking --------------------------------------------------------------------------------------------

#### xcms3 ---------------------------------------------------------------

dtxcms <- dt

dtxcms <- pickingParameters(
  dtxcms,
  algorithm = "xcms3",
  settings = xcms::CentWaveParam(
    ppm = 15, peakwidth = c(5, 60),
    snthresh = 10, prefilter = c(6, 5000),
    mzCenterFun = "mean", integrate = 2,
    mzdiff = -0.0001, fitgauss = TRUE,
    noise = 250, verboseColumns = TRUE,
    firstBaselineCheck = FALSE,
    extendLengthMSW = TRUE
  )
)

pickingParameters(dtxcms)

dtxcms <- peakPicking(dtxcms, save = FALSE)




#### openms --------------------------------------------------------------

dtopenms <- dt

dtopenms <- pickingParameters(
  dtopenms,
  algorithm = "openms",
  settings = list(
    noiseThrInt = 500,
    chromSNR = 10,
    chromFWHM = 10,
    mzPPM = 15,
    reEstimateMTSD = TRUE,
    traceTermCriterion = "sample_rate",
    traceTermOutliers = 5,
    minSampleRate = 0.5,
    minTraceLength = 3,
    maxTraceLength = -1,
    widthFiltering = "fixed",
    minFWHM = 1,
    maxFWHM = 30,
    traceSNRFiltering = FALSE,
    localRTRange = 10,
    localMZRange = 6.5,
    isotopeFilteringModel = "metabolites (5% RMS)",
    MZScoring13C = FALSE,
    useSmoothedInts = TRUE,
    extraOpts = NULL,
    intSearchRTWindow = 3,
    useFFMIntensities = FALSE
  )
)

pickingParameters(dtopenms)

dtopenms <- peakPicking(dtopenms, save = FALSE)


#### inspecting peaks ----------------------------------------------------
#using the S4 method and the mz/rt filtering as presented earlier
mz_01 <- c(247.1651, 239.0628)
rt_01 <- c(839, 937)

peaks(dtxcms, samples = 2, mz = mz_01, ppm = 5, rt = rt_01, sec = 10)

peaks(dtopenms, samples = 2, mz = mz_01, ppm = 5, rt = rt_01, sec = 10)

#Note that targets can be given as peaks and/or features IDs.
#when targets are defined, mz/rt arguments are ignored.

plotPeaks(dtxcms, samples = c(1, 4), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "targets")

plotPeaks(dtopenms, samples = c(1:2, 4:5), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "replicates", interactive = TRUE)

mapPeaks(dtopenms, samples = c(1:2), mz = mz_01[1], ppm = 5, rt = rt_01[1], sec = 10, colorBy = "targets", ylim = 0.02)


### peak grouping and alignment -----------------------------------------------------------------------------

#### xcms3 ---------------------------------------------------------------

dtxcms <- groupingParameters(
  dtxcms,
  algorithm = "xcms3",
  settings = list(
    rtalign = TRUE,
    loadRawData = TRUE,
    groupParam = xcms::PeakDensityParam(
      sampleGroups = "holder",
      bw = 4,
      minFraction = 0.5,
      minSamples = 1,
      binSize = 0.008,
      maxFeatures = 100),
    preGroupParam = xcms::PeakDensityParam(
      sampleGroups = "holder",
      bw = 6,
      minFraction = 0.5,
      minSamples = 1,
      binSize = 0.01,
      maxFeatures = 100),
    retAlignParam = xcms::PeakGroupsParam(
      minFraction = 1,
      extraPeaks = 0,
      smooth = "loess",
      span = 0.3,
      family = "gaussian")
  )
)

groupingParameters(dtxcms)

dtxcms <- peakGrouping(dtxcms, save = FALSE)


#### openms --------------------------------------------------------------

dtopenms <- groupingParameters(
  dtopenms,
  algorithm = "openms",
  settings = list(
    rtalign = TRUE,
    QT = FALSE,
    maxAlignRT = 6,
    maxAlignMZ = 0.01,
    maxGroupRT = 4,
    maxGroupMZ = 0.008,
    extraOptsRT = NULL,
    extraOptsGroup = NULL
  )
)

groupingParameters(dtopenms)

dtopenms <- peakGrouping(dtopenms, save = FALSE)


#### inspecting features -------------------------------------------------

features(dtxcms, samples = 2, mz = mz_01, ppm = 5, rt = rt_01, sec = 10)

features(dtopenms, samples = 2, mz = mz_01, ppm = 5, rt = rt_01, sec = 10)

plotPeaks(dtxcms, samples = c(1, 4), targets = c("M239_R936_135", "M247_R840_150"), colorBy = "targets")

plotFeatures(dtopenms, samples = c(2, 5), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "replicates", interactive = TRUE)

#patRoon option
patRoon::plotChroms(dtxcms@pat[c(2, 4), "M239_R936_135"])









































### Code Trials ---------------------------------------------------------------------------------------------


dtxcms <- fragmentsParameters(dtxcms, algorithm = "ntsiuta")

mz_02 <- data.frame(mzmin = c(207, 233, 242, 748), mzmax = c(208, 234, 243, 748.8))
rt_02 <- data.frame(rtmin = c(931, 939, 880, 880), rtmax = c(935, 944, 888, 888))
targets <- features(dtxcms, mz = mz_02, rt = rt_02)[, id]

MS2 <- generateMS2(dtxcms[which(replicates(dtxcms) == "QC"), targets$id],
  algorithm = "ntsiuta",
  settings = list(
    isolationTimeWindow = 10,
    isolationMassWindow = 1.3,
    clusteringMethod = "distance",
    clusteringUnit = "ppm",
    clusteringWindow = 25,
    minIntensityPre = 300,
    minIntensityPost = 400,
    asPatRoon = TRUE,
    mergeCEs = TRUE,
    mergeBy = NULL
  )
)

MS2[["M233_R940_217"]]$MSMS
MS2[["M207_R932_124"]]$MSMS

MS2[["M233_R940_217"]]$MS
MS2[["M207_R932_124"]]$MS

MS2 <- patRoon::filter(
  MS2,
  isolatePrec = list(
    maxIsotopes = 6,
    mzDefectRange = c(-0.005, 0.005),
    intRange = c(0.0005, 2),
    z = 1,
    maxGap = 2
  ),
  retainPrecursorMSMS = TRUE
)

MS2[["M233_R940_217"]]$MS
MS2[["M207_R932_124"]]$MS


xic_ex1 <- XICs(dtxcms, samples = 7, mz = targets[1], rt = NULL, ppm = 10, sec = 30)


object <- dtxcms[which(replicates(dtxcms) == "QC"), targets$id]

names(MS2@peakLists[[1]][[1]])

names(MS2@metadata[[1]][[1]])

names(MS2@averagedPeakLists[[1]])


plotMS2s(
  object = dtxcms,
  samples = 7,
  mz = targets[2, ], ppm = 20,
  rt = NULL, sec = 60,
  isolationTimeWindow = 10,
  isolationMassWindow = 1.3,
  clusteringMethod = "distance",
  clusteringUnit = "ppm",
  clusteringWindow = 25,
  minIntensityPre = 200,
  minIntensityPost = 200,
  mergeCEs = TRUE,
  mergeBy = NULL,
  colorBy = "targets",
  interactive = TRUE
)


test_sirius <- findFragments(
  dtxcms,
  database = data.table::fread(paste0(path(object), "//tp_nitro.csv")),
  replicates = "QC",
  targets = c("M207_R932_124", "M233_R941_217", "M748_R883_1287", "M242_R884_264"),
  title = NULL,
  ppm = 15,
  ppmLoss = 40,
  minFeatureIntensity = 5000,
  inSilico = "sirius",
  MS2settings = NULL
)

test_genform <- findFragments(
  dtxcms,
  database = data.table::fread(paste0(path(object), "//tp_nitro.csv")),
  replicates = "QC",
  targets = c("M207_R932_124", "M233_R941_217", "M748_R883_1287", "M242_R884_264"),
  title = NULL,
  ppm = 15,
  ppmLoss = 40,
  minFeatureIntensity = 5000,
  inSilico = "genform",
  MS2settings = NULL
)

test_sirius@workflows[[1]]@results

test_genform@workflows[[1]]@results

