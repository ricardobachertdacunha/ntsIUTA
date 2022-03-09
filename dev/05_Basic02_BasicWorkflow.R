

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
    snthresh = 10, prefilter = c(6, 1000),
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


#### kpic2 --------------------------------------------------------------

dtkpic2 <- dt

dtkpic2 <- pickingParameters(
  dtkpic2,
  algorithm = "kpic2",
  settings = list(
    kmeans = FALSE,
    level = 1000,
    mztol = 0.01,
    gap = 3,
    width = c(5),
    #alpha = 0.3,
    min_snr = 4,
    parallel = TRUE
  )
)

pickingParameters(dtkpic2)

dtkpic2 <- peakPicking(dtkpic2, save = FALSE)


##### safd ----------------------------------------------------------------------
# Only mzXML data in profile mode are possible

dtsafd <- dt

dtsafd <- pickingParameters(
  dtsafd,
  algorithm = "safd",
  settings = list(
    profPath = dirname(dtsafd@samples$file),
    mzRange = c(0, 1200),
    maxNumbIter = 1000,
    maxTPeakW = 300,
    resolution = 12000,
    minMSW = 0.02,
    RThreshold = 0.75,
    minInt = 500,
    sigIncThreshold = 5,
    S2N = 3,
    minPeakWS = 3
  )
)

pickingParameters(dtsafd)

dtsafd <- peakPicking(dtsafd, save = FALSE)


##### sirius --------------------------------------------------------------------
#Note, note working with trimed mzML files

dtsirius <- dt

dtsirius <- pickingParameters(
  dtsirius,
  algorithm = "sirius",
  settings = list()
)

pickingParameters(dtsirius)

dtsirius <- peakPicking(dtsirius, save = FALSE)



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

plotPeaks(dtxcms, samples = c(1, 4), targets = c("M239_R936_251", "M247_R840_281"), colorBy = "targets")

plotFeatures(dtxcms, samples = c(2, 5), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "replicates", interactive = TRUE)

plotFeatures(dtopenms, samples = c(2, 5), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "replicates", interactive = TRUE)

#patRoon option
patRoon::plotChroms(dtxcms@pat[c(2, 4), "M239_R936_135"])

#### inspecting alignment ------------------------------------------------

plotAlignment(dtxcms)


### peak filling --------------------------------------------------------------------------------------------

#### xcms ----------------------------------------------------------------

dtxcms <- fillingParameters(
  dtxcms,
  algorithm = "xcms",
  settings = xcms::ChromPeakAreaParam(
    mzmin = function(z) quantile(z, probs = 0.25),
    mzmax = function(z) quantile(z, probs = 0.75),
    rtmin = function(z) quantile(z, probs = 0.25),
    rtmax = function(z) quantile(z, probs = 0.75)
  )
)

fillingParameters(dtxcms)

dtxcms <- peakFilling(dtxcms, save = FALSE)

#### openms --------------------------------------------------------------

dtopenms <- fillingParameters(
  dtopenms,
  algorithm = "xcms",
  settings = xcms::FillChromPeaksParam(
    expandMz = 0,
    expandRt = 0,
    ppm = 0,
    fixedMz = 0,
    fixedRt = 0
  )
)

fillingParameters(dtopenms)

dtopenms <- peakFilling(dtopenms, save = FALSE)


#### inspecting filling --------------------------------------------------

showfill <- features(dtxcms)[hasFilled == TRUE, ]
setorder(showfill, -IN)
showfill <- showfill[1:5, ]

plotFeaturePeaks(
  object <- dtxcms,
  samples = c(1:6),
  targets = showfill$id,
  mz = NULL, ppm = 20,
  rt = NULL, sec = 30,
  legendNames = NULL
)

peaks(dtxcms, targets = showfill$id)


### annotation ----------------------------------------------------------------------------------------------

dt_cliquems <- annotationParameters(
  dtxcms,
  algorithm = "cliquems",
  settings = list(
    ionization = "positive",
    maxCharge = 3,
    maxGrade = 5,
    ppm = 20,
    adductInfo = NULL,
    absMzDev = 0.01,
    minSize = 2,
    relMinAdductAbundance = 0.5,
    adductConflictsUsePref = TRUE,
    NMConflicts = c("preferential"),
    prefAdducts = c("[M+H]+"),
    extraOptsCli = NULL,
    extraOptsIso = NULL,
    extraOptsAnn = NULL,
    parallel = TRUE
  )
)

dt_ramclustr <- annotationParameters(
  dtxcms,
  algorithm = "ramclustr",
  settings = list(
    ionization = "positive",
    st = NULL,
    sr = NULL,
    maxt = 12,
    hmax = 0.3,
    normalize = "TIC",
    absMzDev = 0.01,
    relMzDev = 20,
    minSize = 2,
    relMinReplicates = 0.5,
    RCExperimentVals = list(design = list(platform = "LC-MS"), instrument =
      list(ionization = "positive", MSlevs = 1)),
    extraOptsRC = NULL,
    extraOptsFM = NULL
  )
)

dt_xcms <- annotationParameters(
  dtxcms,
  algorithm = "camera",
  settings = list(
    ionization = polarity(dtxcms),
    onlyIsotopes = FALSE,
    minSize = 2,
    relMinReplicates = 0.5,
    extraOpts = list(
      sigma = 6,
      perfwhm = 0.35,
      cor_eic_th = 0.3,
      graphMethod = "hcs",
      pval = 0.05,
      calcCiS = TRUE,
      calcIso = TRUE,
      calcCaS = TRUE,
      maxcharge = 3,
      maxiso = 5,
      ppm = 20,
      mzabs = 0.01,
      rules =  data.table::fread(system.file("rules/primary_adducts_pos.csv", package = "CAMERA"), header = TRUE),
      multiplier = 3,
      max_peaks = 500,
      intval = "maxo"
    )
  )
)

dt_cliquems <- peakAnnotation(dt_cliquems, save = FALSE)
dt_ramclustr <- peakAnnotation(dt_ramclustr, save = FALSE)
dt_xcms <- peakAnnotation(dtxcms, save = FALSE)

#### inspect annotation --------------------------------------------------

tg <- features(dt_xcms, targets = "M239_R936_638")

#cliqueMS
eval_cliquems <- features(dt_cliquems)
eval_cliquems <- eval_cliquems[neutralMass %in% tg$neutralMass, ]

#ramclustr
eval_ramclustr <- features(dtramclustr)
eval_ramclustr <- eval_ramclustr[neutralMass %in% tg$neutralMass, ]

#xcms
eval_xcms <- features(dt_xcms)
eval_xcms <- eval_xcms[neutralMass %in% tg$neutralMass, ]
eval_xcms_p <- peaks(dt_xcms)
eval_xcms_p <- eval_xcms_p[feature %in% annoInsp$id, ]

mapPeaks(object, targets = eval_xcms_p$id)





































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


xcms::plotAdjustedRtime(object@pat@xdata) # col = group_colors[object@pat@xdata$sample_group]


rtAdj <- xcms::adjustedRtime(object@pat@xdata)
head(rtAdj)


pkAdj <- xcms::processHistory(object@pat@xdata)[[3]]
pkAdj <- pkAdj@param
pkAdj <- xcms::peakGroupsMatrix(pkAdj)



rtAdj_dt <- object@scans
rtAdj_dt <- lapply(rtAdj_dt, function(x) x[, .(seqNum, acquisitionNum, retentionTime)])
names(rtAdj_dt) <- samples(object)
rtAdj_dt <- lapply(seq_len(length(rtAdj_dt)), function(x, rtAdj, rtAdj_dt) {

  rts <- names(rtAdj)
  rts <- stringr::str_detect(rts, paste0("F", x))
  rts <- rtAdj[rts]

  rtAdj_dt[[x]][, adjustedRetentionTime := rts]
  rtAdj_dt[[x]][, adjustment := adjustedRetentionTime - retentionTime]

  return(rtAdj_dt[[x]])

}, rtAdj = rtAdj, rtAdj_dt = rtAdj_dt)


plot(rtAdj_dt[[4]]$retentionTime, rtAdj_dt[[4]]$adjustment, type = "l")

points(
  unique(pkAdj[, 4][pkAdj[, 4] %in% rtAdj_dt[[4]]$retentionTime]),
  unique(rtAdj_dt[[4]]$adjustment[rtAdj_dt[[4]]$retentionTime %in% pkAdj[, 4]])
)


length(pkAdj[, 4])
length(rtAdj_dt[[4]]$adjustment[rtAdj_dt[[4]]$adjustedRetentionTime %in% pkAdj[, 4]])

rtAdj_dt[[1]]

View(xcms::adjustedRtime(object@pat@xdata))

View(test)
test1 <- test@phenoData

test <- xdata@msFeatureData

names(xdata@msFeatureData$adjustedRtime)

