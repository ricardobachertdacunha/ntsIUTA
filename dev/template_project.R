
### Library -------------------------------------------------------------------------------------------------
library(ntsIUTA)


### Create project ------------------------------------------------------------------------------------------
# Only when ntsData is not yet created

#### Demo project path ---------------------------------------------------
path <- system.file(package = "ntsIUTA", dir = "extdata")

#### Choose a project path -----------------------------------------------
path <- utils::choose.dir(getwd(), "Select or create a project folder")

#### Setup project -------------------------------------------------------
dt <- setupProject(
  path = path,
  title = "Project name",
  description = "Project desription",
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


#### Convert raw files -------------------------------------------------------
# When raw MS files in project path need to be converted to mzML
mzMLconverter(
  path = path(dt),
  files = NULL,
  convertFrom = "agilent",
  centroidMethod = "vendor",
  outPath = NULL,
  overwrite = TRUE
)


#### Add (more) files ----------------------------------------------------
# Only when more files need to be added
dt <- addFiles(
  files = utils::choose.files(),
  object = dt,
  copy = FALSE,
  replicates = NULL,
  polarity = "positive",
  method = NULL
)

#### Save ntsData --------------------------------------------------------
dt <- dt[1:18] #filter files not needed for demo project
saveObject(dt = dt)



### Load ntsData --------------------------------------------------------------------------------------------
# When ntsData is already created
path <- system.file(package = "ntsIUTA", dir = "extdata") #only for demo project
dt <- readRDS(paste0(path, "rData/ntsData.rds"))



### Update/Correct ntsData -----------------------------------------------------------------------------------

#### Project info -------------------------------------------------
dt <- projectInfo(
  dt,
  title = "Demo Project",
  description = "Description example",
  date = Sys.Date()
)

#### Replicates ----------------------------------------------------------
replicates(dt) <- c(
  rep("Blank", 3),
  rep("Replicate1", 3),
  rep("Replicate2", 3),
  rep("Replicate3", 3),
  rep("Replicate4", 3),
  rep("QC", 3)
)

#### Assign blanks -------------------------------------------------------
blanks(dt) <- rep("Blank", 18)

#### Polarity ------------------------------------------------------------
polarity(dt) <- "positive"

#### Acquisition Method --------------------------------------------------
acquisitionMethods(dt) <- "NTS_MethodName"

#### Experimental metadata -----------------------------------------------
var <- data.table::data.table(
  matrix = c(
    rep("clean", 3),
    rep("dirty", 12),
    rep("control", 3)
  )
)

dt <- addMetadata(dt, var)

#### Assign QCs ----------------------------------------------------------
QC(dt) <- "QC"



### TICs ----------------------------------------------------------------------------------------------------
plotTICs(dt, samples = NULL, colorBy = "replicates")


### Basic Workflow ------------------------------------------------------------------------------------------

#### Settings ------------------------------------------------------------
dt <- pickingParameters(
  dt,
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

dt <- fillingParameters(
  dt,
  algorithm = "xcms",
  settings = xcms::ChromPeakAreaParam(
    mzmin = function(z) quantile(z, probs = 0.25),
    mzmax = function(z) quantile(z, probs = 0.75),
    rtmin = function(z) quantile(z, probs = 0.25),
    rtmax = function(z) quantile(z, probs = 0.75)
  )
)


#### Check QC ------------------------------------------------------------
# Before running the basic workflow for the project samples, the QC samples should be checked






