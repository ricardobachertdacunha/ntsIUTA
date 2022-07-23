

### setup ---------------------------------------------------------------------------------------------------

library(ntsIUTA)
devtools::load_all()

path <- "C:\\Users\\Ricardo\\Documents\\R_DemoProject"

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

dt <- dt[1:9]

replicates(dt) <- c(
  rep("QC", 3),
  rep("Blank", 3),
  rep("IN", 3)
)

blanks(dt) <- rep("Blank", 9)

samplesTable(dt)




### add metadata --------------------------------------------------------------------------------------------

var <- data.table::data.table(matrix = c(rep("control", 3), rep("clean", 3), rep("dirty", 3)))
dt <- addMetadata(dt, var)

metadata(dt)


### check TICs ----------------------------------------------------------------------------------------------

tic <- TICs(dt, samples = NULL)

plotTICs(tic, samples = NULL, colorBy = "replicates")




### peak picking --------------------------------------------------------------------------------------------

dt <- pickingParameters(
  dt,
  algorithm = "xcms3",
  settings = xcms::CentWaveParam(
    ppm = 15, peakwidth = c(5, 60),
    snthresh = 5, prefilter = c(6, 500),
    mzCenterFun = "mean", integrate = 2,
    mzdiff = -0.0001, fitgauss = TRUE,
    noise = 250, verboseColumns = TRUE,
    firstBaselineCheck = FALSE,
    extendLengthMSW = TRUE
  )
)

pickingParameters(dt)

dt <- peakPicking(dt, save = FALSE)


#### inspecting peaks ----------------------------------------------------

#using the S4 method and the mz/rt filtering as presented earlier
mz_01 <- c(247.1651, 239.0628)
rt_01 <- c(839, 937)

peaks(dt, samples = c(1:6), mz = mz_01, ppm = 5, rt = rt_01, sec = 10)[, 1:10]

plotPeaks(dt, samples = c(1:6), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "samples")
plotPeaks(dt, samples = c(1:6), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "targets", interactive = TRUE)

mapPeaks(dt, samples = c(1:6), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "targets", ylim = 0.02)

#Note that targets can be given as peaks and/or features IDs/UFIs.
#when targets are defined, mz/rt arguments are ignored.




### peak grouping and alignment -----------------------------------------------------------------------------

dt <- groupingParameters(
  dt,
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

groupingParameters(dt)

dt <- peakGrouping(dt, save = FALSE)


#### inspecting features -------------------------------------------------

fts <- features(dt, samples = (1:6), mz = mz_01, ppm = 5, rt = rt_01, sec = 10)
targets <- fts$id

plotFeatures(dt, samples = c(1:6), targets = targets, colorBy = "targets")
plotFeatures(dt, samples = c(1:6), targets = targets, colorBy = "replicates", interactive = TRUE)

plotFeaturePeaks(dt, samples = c(1:6), targets = targets, heights = c(0.7, 0.3))

#patRoon option
patRoon::plotChroms(dt@pat[c(2, 4), targets])


#### inspecting alignment ------------------------------------------------

plotAlignment(dt)




### peak filling --------------------------------------------------------------------------------------------

dt_fill <- fillingParameters(
  dt,
  algorithm = "xcms",
  settings = xcms::ChromPeakAreaParam(
    mzmin = function(z) quantile(z, probs = 0.25),
    mzmax = function(z) quantile(z, probs = 0.75),
    rtmin = function(z) quantile(z, probs = 0.25),
    rtmax = function(z) quantile(z, probs = 0.75)
  )
)

fillingParameters(dt_fill)

dt_fill <- peakFilling(dt_fill, save = FALSE)


#### inspecting filling --------------------------------------------------

showfill <- features(dt_fill)[hasFilled == TRUE, ]
setorder(showfill, -IN)
showfill <- showfill[1:5, ]

plotFeaturePeaks(dt_fill, samples = c(1:6), targets = showfill$id)




### annotation ----------------------------------------------------------------------------------------------

dt <- annotationParameters(
  dt,
  algorithm = "camera",
  settings = list(
    ionization = polarity(dt),
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

annotationParameters(dt)

dt <- peakAnnotation(dt, save = FALSE)


#### inspect annotation --------------------------------------------------

components(dt, targets = targets, all = TRUE)

plotComponents(dt, targets = targets, all = FALSE, colorBy = "isotopes")




### unique feature identifier -------------------------------------------------------------------------------

dt <- makeUFI(dt)

ufi_ex <- features(dt,
  mz = data.frame(
    mz = c(239.0628, 247.1651, 213.1869, 267.0698, 275.23346),
    rt = c(936, 840, 930, 720, 624)
  )
)

ufi_ex[, .(ufi, id, mz, rt)]
dt@unified[id %in% targets[1], ]



### wrapper workflow ----------------------------------------------------------------------------------------

#### save and load parameters --------------------------------------------

saveParameters(dt)

dt <- loadParameters(dt)


#### wrapper function ----------------------------------------------------

dt2 <- createFeatures(dt, save = FALSE)




### filter features -----------------------------------------------------------------------------------------

dt <- filter(dt,
  minIntensity = 5000,
  blankThreshold = 3,
  maxReplicateIntensityDeviation = 40,
  minReplicateAbundance = 3
)

nrow(features(dt))

dt <- removeFilteredFeatures(dt)

nrow(features(dt))

dt <- restoreFilteredFeatures(dt)

nrow(features(dt))

dt <- removeFilteredFeatures(dt)



### calculate quality of features ---------------------------------------------------------------------------

#### SNR -----------------------------------------------------------------

sn_targets <- features(dt)
sn_targets <- head(sn_targets, 700)
sn_targets <- sn_targets$id

dt <- calculateSNR(dt, targets = sn_targets, rtExpand = 200)

features(dt)[1:5]

dt <- filter(
  dt,
  snRatio = 3
)

dt <- removeFilteredFeatures(dt)

#### MetaClean -----------------------------------------------------------
dt <- calculateFeaturesMetadata(dt, targets = targets)
features(dt, targets = targets)




### subset on features --------------------------------------------------------------------------------------

dt3 <- dt[7:9, features(dt)[, id]]

samples(dt)
samples(dt3)

#Note that after subseting on samples might be good to update the feature table
#as the subset only does a lazy update (i.e., only the averaged mz, rt and intensity as updated)
dt3 <- updateFeatureTable(dt3, fast = FALSE)

### load MS2 ------------------------------------------------------------------------------------------------

dt3 <- fragmentsParameters(
  dt3,
  algorithm = "ntsiuta",
  settings = list(
    isolationTimeWindow = 5,
    isolationMassWindow = 1.3,
    clusteringMethod = "distance",
    clusteringUnit = "ppm",
    clusteringWindow = 25,
    minIntensityPre = 250,
    minIntensityPost = 300,
    asPatRoon = TRUE,
    mergeCEs = TRUE,
    mergeBy = NULL
  )
)

fragmentsParameters(dt3)

dt3 <- loadMS2(dt3)


#### inspection ----------------------------------------------------------

dt3@ms2[["M242_R884_1145"]]$MS

test <- features(dt3)


patRoon::plotSpectrum(
  dt3@ms2,
  groupName = c("M242_R884_1145", "M242_R884_1145"), #, "M326_R706_2656"
  analysis1 = 1,
  analysis2 = 2,
  MSLevel = 2,
  title = NULL,
  specSimParams = patRoon::getDefSpecSimParams(),
  xlim = NULL,
  ylim = NULL
)



showMethods("plotSpectrum")



















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


fragmentsParameters

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
