

### setup ---------------------------------------------------------------------------------------------------

library(ntsIUTA)
devtools::load_all()

path <- "C:\\Users\\Ricardo\\Documents\\R_NTS_article"

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

rm(path)

replicates(dt) <- c(
  rep("orb_mzml_cent_01", 1),
  rep("orb_mzml_cent_1", 1),
  rep("orb_mzml_cent_10", 1),
  rep("orb_mzml_prof_10", 1),
  rep("orb_mzxml_cent_10", 1),
  #rep("tof_mzml_cent_10", 3),
  rep("trim_orb_mzml_cent_10", 1),
  rep("trim_orb_mzml_prof_10", 1),
  rep("trim_orb_mzxml_cent_10", 1)
)

samplesTable(dt)



### add metadata --------------------------------------------------------------------------------------------

var <- data.table::data.table(system = c(rep("orbitrap", 5), rep("tof", 0), rep("orbitrap", 3)))
dt <- addMetadata(dt, var)
rm(var)
metadata(dt)




### check TICs ----------------------------------------------------------------------------------------------

tic <- TICs(dt[c(4:6)], samples = NULL)
plotTICs(tic, samples = NULL, colorBy = "samples", interactive = TRUE)




### database ------------------------------------------------------------------------------------------------

db <- data.table::fread("C:/Users/Ricardo/Documents/CodeProjects/ntsIUTA/inst/extdata/suspectListTemplate.csv")
targets <- data.table::data.table(id = db$name, mz = db$mz)




### EICs ----------------------------------------------------------------------------------------------------

#### orbitrap ------------------------------------------------------------

orb_targets <- copy(targets)
orb_eics <- EICs(dt, samples = 3, mz = orb_targets, ppm = 5)
plotEICs(orb_eics[intensity > 0, ], interactive = TRUE, title = "orbitrap", colorBy = "targets")

orb_targets[, rt := 0]
for (i in orb_targets$id) {
  orb_targets[id %in% i, rt := orb_eics[id %in% i & intensity == max(orb_eics[id %in% i, intensity]), rt]]
}

setorder(orb_targets, rt)
#write.csv(orb_targets, "C:\\Users\\Ricardo\\Documents\\CodeProjects\\ntsIUTA\\inst\\extdata\\mix1_RTs.csv")

#### tof -----------------------------------------------------------------

tof_targets <- copy(targets)
tof_eics <- EICs(dt, samples = 8, mz = tof_targets, ppm = 10)
plotEICs(tof_eics, interactive = TRUE, title = "tof")

tof_targets[, rt := 0]
for (i in tof_targets$id) {
  tof_targets[id %in% i, rt := tof_eics[id %in% i & intensity == max(tof_eics[id %in% i, intensity]), rt]]
}

setorder(tof_targets, rt)

#### rt comparison -------------------------------------------------------

mergedTargets <- dplyr::left_join(orb_targets[, .(id, mz, rt)], tof_targets[, .(id, rt)], by = "id")
setnames(mergedTargets, c("id", "mz", "orb_rt", "tof_rt"))
mergedTargets[, `:=`("diff" = orb_rt - tof_rt)]
mergedTargets <- mergedTargets[!id %in% "Ibuprofen", ]

plot(mergedTargets$orb_rt, type = "n",
     xlim = c(min(mergedTargets$orb_rt), max(mergedTargets$orb_rt)),
     ylim = c(min(mergedTargets$diff), max(mergedTargets$tof_rt)))
points(x = mergedTargets$orb_rt, y = mergedTargets$tof_rt, col = "black")
points(x = mergedTargets$orb_rt, y = mergedTargets$diff, col = "red")




#### mz & rt dev --------------------------------------------------------------------------------------------

orb_XIC <-  XICs(dt, samples = c(5), mz = orb_targets[5, ], ppm = 20)
tof_XIC <-  XICs(dt, samples = c(8), mz = tof_targets[5, ], ppm = 20)

min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

orb_XIC$intensity <- min_max_norm(orb_XIC$intensity)
tof_XIC$intensity <- min_max_norm(tof_XIC$intensity)

XIC <- rbind(orb_XIC, tof_XIC)
plotXICs(XIC)


#### tof high res -------------------------------------------------------------------------------------------

# Systematic eval of the different tunning modes
dt2 <- setupProject(save = FALSE, makeNewProject = FALSE)


temp_eics <- EICs(dt2, samples = 1:2, mz = targets, ppm = 10)
plotEICs(temp_eics, interactive = TRUE, title = "temp", colorBy = "targets")



high_XIC <-  XICs(dt2, samples = c(6), mz = tof_targets[16, ], ppm = 20)
norm_XIC <-  XICs(dt2, samples = c(4), mz = tof_targets[16, ], ppm = 20)
XIC2 <- rbind(high_XIC, norm_XIC)
plotXICs(XIC2)




### peak picking --------------------------------------------------------------------------------------------

#### trim files ----------------------------------------------------------

#Files were trimmed by mzR using the code in 05_Basic04_trimMSfiles.R
#Samples 3, 4 and 5 were trimmed as requested for each algorithm to be tested

t_dt <- dt[6:8]
t_orb_targets <- orb_targets[rt > 300 & rt < 700, ]

trim_orb_eics <- EICs(t_dt, samples = NULL, mz = t_orb_targets, ppm = 10)
plotEICs(trim_orb_eics[intensity > 0, ], interactive = TRUE, title = "orbitrap", colorBy = "samples")





#### optimization -------------------------------------------------------------------------------------------












#### xcms3 ---------------------------------------------------------------

dtxcms <- t_dt[1]

dtxcms <- pickingParameters(
  dtxcms,
  algorithm = "xcms3",
  settings = xcms::CentWaveParam(
    ppm = 3, peakwidth = c(5, 120),
    snthresh = 10, prefilter = c(6, 60000),
    mzCenterFun = "mean", integrate = 2,
    mzdiff = -0.001, fitgauss = TRUE,
    noise = 20000, verboseColumns = TRUE,
    firstBaselineCheck = TRUE,
    extendLengthMSW = TRUE
  )
)

dtxcms <- peakPicking(dtxcms, save = FALSE)


pks <- peaks(dtxcms, mz = orb_targets, ppm = 3, rt = 5)

plotPeaks(dtxcms, targets = pks[1:3, id], interactive = TRUE)








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

dtsirius <- dt[22]

dtsirius <- pickingParameters(
  dtsirius,
  algorithm = "sirius",
  settings = list()
)

pickingParameters(dtsirius)

dtsirius <- peakPicking(dtsirius, save = FALSE)

pks <- peaks(dtsirius, mz = orb_targets, ppm = 5, rt = 10)
plotPeaks(dtsirius, targets = pks[5:6, id], interactive = TRUE)


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

plotPeaks(dtxcms, samples = c(1, 4), targets = c("M239_R936_638", "M247_R840_725"), colorBy = "targets")

plotFeatures(dtxcms, samples = c(1, 4), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "replicates", interactive = TRUE)

plotFeatures(dtopenms, samples = c(2, 5), mz = mz_01, ppm = 5, rt = rt_01, sec = 10, colorBy = "replicates", interactive = TRUE)

plotFeaturePeaks(
  object <- dtxcms,
  samples = c(1:6),
  targets = c("M239_R936_638", "M247_R840_725"),
  mz = NULL, ppm = 20,
  rt = NULL, sec = 30,
  legendNames = NULL,
  heights = c(0.8, 0.2)
)

#patRoon option
patRoon::plotChroms(dtxcms@pat[c(2, 4), "M239_R936_135"])

#### inspecting alignment ------------------------------------------------

plotAlignment(dtxcms)


### peak filling --------------------------------------------------------------------------------------------

#### xcms ----------------------------------------------------------------

dtxcms_fill <- fillingParameters(
  dtxcms,
  algorithm = "xcms",
  settings = xcms::ChromPeakAreaParam(
    mzmin = function(z) quantile(z, probs = 0.25),
    mzmax = function(z) quantile(z, probs = 0.75),
    rtmin = function(z) quantile(z, probs = 0.25),
    rtmax = function(z) quantile(z, probs = 0.75)
  )
)

fillingParameters(dtxcms_fill)

dtxcms_fill <- peakFilling(dtxcms_fill, save = FALSE)

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

showfill <- features(dtxcms_fill)[hasFilled == TRUE, ]
setorder(showfill, -IN)
showfill <- showfill[1:5, ]

plotFeaturePeaks(
  object <- dtxcms_fill,
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
dt_xcms <- peakAnnotation(dt_xcms, save = FALSE)

#### inspect annotation --------------------------------------------------

#cliqueMS
tg <- features(dt_cliquems, targets = "M239_R936_638")
eval_cliquems <- features(dt_cliquems)
eval_cliquems <- eval_cliquems[neutralMass %in% tg$neutralMass, ]
eval_cliquems

#ramclustr
tg <- features(dt_ramclustr, targets = "M239_R936_638")
eval_ramclustr <- features(dt_ramclustr)
eval_ramclustr <- eval_ramclustr[neutralMass %in% tg$neutralMass, ]
eval_ramclustr

#xcms
tg <- features(dt_xcms, targets = c("M239_R936_638",  "M247_R840_725"))
eval_xcms <- features(dt_xcms)
eval_xcms <- eval_xcms[neutralMass %in% tg$neutralMass, ]
eval_xcms

eval_xcms_p <- peaks(dt_xcms)
eval_xcms_p <- eval_xcms_p[feature %in% eval_xcms$id, ]

mapPeaks(dt_xcms, targets = eval_xcms_p$id)

plotFeatures(dt_xcms, targets = unique(sort(eval_xcms_p$feature)), interactive = TRUE)

targets <- eval_xcms$id




### unique feature identifier -------------------------------------------------------------------------------

dt_xcms <- makeUFI(dt_xcms)

targets <- features(dt_xcms,
  mz = data.frame(
    mz = c(239.0628, 247.1651, 213.1869, 267.0698, 275.23346),
    rt = c(936, 840, 930, 720, 624)
  )
)

targets[, .(ufi, id, mz, rt)]

test <- features(dt_xcms)






### calculate quality of features ---------------------------------------------------------------------------

#### SNR -----------------------------------------------------------------

object <- dt_xcms
targets <- features(object)
targets <- head(targets, 500)
targets <- targets$id

object <- calculateSNR(object, targets = targets, rtExpand = 200)


#### MetaClean -----------------------------------------------------------
object <- calculateFeaturesMetadata(object, targets = targets)

corQuality <- features(object)
plot(log(corQuality$sn_value), corQuality$totalScore, type = "p")
abline(v = log(3))
plot(log(corQuality$sn_value), corQuality$GaussianSimilarityScore, type = "p")
abline(v = log(3))
plot(log(corQuality$sn_value), corQuality$ZigZagScore, type = "p")
abline(v = log(3))
plot(log(corQuality$sn_value), corQuality$SharpnessScore, type = "p")
abline(v = log(3))


### filter features -----------------------------------------------------------------------------------------

object2 <- filter(object,
  minIntensity = 5000,
  blankThreshold = 3,
  maxReplicateIntensityDeviation = 40,
  minReplicateAbundance = 3,
  snRatio = 3
)

test <- features(object2)

object3 <- removeFilteredFeatures(object2)

features(object3)

object3 <- restoreFilteredFeatures(object3)




























### Code Trials ---------------------------------------------------------------------------------------------

dt <- fragmentsParameters(
  dt,
  algorithm = "ntsiuta",
  settings = list(
    isolationTimeWindow = 5,
    isolationMassWindow = 0.005,
    clusteringMethod = "distance",
    clusteringUnit = "ppm",
    clusteringWindow = 2,
    minIntensityPre = 1000,
    minIntensityPost = 2000,
    asPatRoon = FALSE,
    mergeVoltages = TRUE,
    mergeBy = NULL
  )
)

plotMS2s(dt, samples = 1, mz = db[name %in% "Sotalol", mz], ppm = 5, interactive = TRUE)





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
