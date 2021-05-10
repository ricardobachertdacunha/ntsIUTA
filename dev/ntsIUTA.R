## ----settings, include = FALSE------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(magrittr)
library(plotly)
library(kableExtra)


## ----pkgLoad------------------------------------------------------------------
library(ntsIUTA)

## ----patRoon installation, eval=FALSE-----------------------------------------
#  writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
#  install.packages("remotes")
#  install.packages("BiocManager")
#  BiocManager::install(c("CAMERA","mzR","xcms","Biobase","MSnbase"))
#  remotes::install_github("rickhelmus/patRoon", upgrade = "never")
#  

## ----install_ntsIUTA, eval=FALSE----------------------------------------------
#  remotes::install_github("ricardobachertdacunha/ntsIUTA", auth_token = "<auth_token>", upgrade = "never")
#  
#  #or once the package goes public
#  remotes::install_github("ricardobachertdacunha/ntsIUTA")
#  

## ----setup--------------------------------------------------------------------
projPath <- system.file(package = "ntsIUTA",dir = "extdata")
setup <- ntsIUTA::makeSetup(projPath, save = F, makeNewProject = F)
#Filter for selecting the six first files in the sampleInfo
setup$sampleInfo <- setup$sampleInfo[1:6,]

## ----sampleInfoExample, cache=FALSE, echo=FALSE, message=FALSE, results='axis'----
sampleInfo_temp <- setup$sampleInfo
sampleInfo_temp$filePath <- "C:/FilePath/FileName.mzML"
sampleInfo_temp %>% kbl(caption = "sampleInfo Example") %>% kable_paper("hover", full_width = F) %>% kable_styling(position = "left")

## ----showPeakPickingObj, warning=FALSE----------------------------------------
peaksData[[2]]

## -----------------------------------------------------------------------------
FT_diuron <- xcms::featureDefinitions(featData,
                                      mz = 233.0243, ppm = 5,
                                      rt = c(15.7*60-10,15.7*60+10),
                                      type = "within", msLevel = 1)
FT_diuron

## ----chromPeaksInFeat---------------------------------------------------------
xcms::chromPeaks(featData, isFilledColumn = TRUE)[base::unlist(FT_diuron[ ,"peakidx"]), ]

## -----------------------------------------------------------------------------
featComp

## ----checkComponents----------------------------------------------------------
ntsIUTA::checkComponents(xA = featComp, replicateGroups = 1,
                         features = NULL, featData = NULL,
                         mz = 233.0243, ppm = 5,
                         rt = 15.7, rtWindow = 0.8,
                         rtUnit = "min",
                         onlyRelated = TRUE)

## ----covertToPat--------------------------------------------------------------
patData <- ntsIUTA::getPatData(featData,
                               setup$sampleInfo,
                               save = FALSE)

## ----seeFl--------------------------------------------------------------------
base::t(fl[fl$FT == "FT0107",])

## -----------------------------------------------------------------------------


## ----prepQC-------------------------------------------------------------------
#screeningList of the reference standards present in the sample
sl_qc <- base::paste0(system.file(package = "ntsIUTA", dir = "extdata"),"/QC_ScreeningList_ntsIUTA_MS2_pos.csv")
utils::head(utils::read.csv(sl_qc), n = 2)

## ----screeningListtemplate, eval=FALSE, cache=FALSE---------------------------
#  ntsIUTA::getScreeningListTemplate()

## ----QCdf---------------------------------------------------------------------
utils::head(QC$df, n = 2)

