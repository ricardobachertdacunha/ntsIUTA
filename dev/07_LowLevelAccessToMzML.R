
### Get Paths -----

pathLowLevel <- system.file(package = "ntsIUTA", dir = "extdata")

msFilesLowLevel <- list.files(path = pathLowLevel,
  pattern = ".mzML|.mzXML",
  recursive = TRUE,
  full.names = TRUE,
  no.. = TRUE)

setup <- ntsIUTA::setupProject(path = pathLowLevel,
                               makeNewProject = FALSE,
                               save = FALSE)
setup <- setup[1:3]
system.time(setup <- ntsIUTA::importRawData(setup))

### Goal -------

#Extract MS1 EIC and MS2 data for Diuron d3 with m/z of 239.0628 and rt of 937 sec.
#The data is extracted for the blank replicate files (x3).
#EIC is extracted with +/- 20 ppm and +/- 30 sec.

mz <- 239.0628
rt <- 937
ppm <- 20
rtw <- 30
mzr <- c(mz - ((ppm / 1E6) * mz), mz + ((ppm / 1E6) * mz))
rtr <- c(rt - rtw, rt + rtw)

### with mzR -------------------------------------------------------------------


#BiocManager::install("sneumann/mzR", ref = "feature/updatePwiz_3_0_21263")


library(mzR)
library(data.table)

ms <- as.list(msFilesLowLevel[1])
system.time(ms <- lapply(ms, function(x) openMSfile(x, backend = "pwiz")))
hd <- lapply(ms, header)
hd <- lapply(hd, as.data.table)

#### get MS1 EIC -----


ms <- mzR::openMSfile("C:/Users/Ricardo/Documents/R_Dev_mzML_example/QC-r001.mzML")
hd <- data.table::as.data.table(mzR::header(ms))

runInfo(ms)


instrumentInfo(ms)
## Individual accessors from instrumentInfo
# ionisation(ms)
# softwareInfo(ms)
# sampleInfo(ms)
# sourceInfo(ms)
# model(ms)
# analyzer(ms)
# detector(ms)

# mzidInfo(ms)
# modifications(ms)
# psms(ms)
# tolerance(ms)
# score(ms)
# para(ms)
# specParams(ms)


EIC <- chromatogramHeader(ms[[1]])

nChrom(ms[[1]])

head(tic(ms[[1]]))

class(ms[[1]])

mzR::chromatogram(ms[[1]])

mzR::peaks(ms[[1]], 1:2)

mzR::chromatogram()


system.time(
EIC <- lapply(seq_len(length(hd)), function(x, ms, hd, rtr, mzr) {
  temp <- hd[[x]][between(hd[[x]][, retentionTime], rtr[1], rtr[2]), ]
  rtv <- temp[, retentionTime]
  scansN <- temp[, seqNum]
  spt <- mzR::peaks(ms[[x]], scansN)
  names(spt) <- round(rtv, digits = 3)
  spt <- parallel::mclapply(spt, function(z, mzr) {
    dt <- data.table(mz = z[, 1], int = z[, 2])
    dt <- dt[between(dt[, mz], mzr[1], mzr[2]), ]
  }, mzr = mzr)
  spt <- rbindlist(as.list(spt), idcol = "rt")
  spt <- spt[, file := x]
  return(spt)
}, ms = ms, hd = hd, rtr = rtr, mzr = mzr)
)

EIC[[1]]

system.time(
EIC_ntsIUTA <- ntsIUTA::extractEIC(setup,
                                   mz = mz,
                                   ppm = ppm,
                                   rt = rt,
                                   rtWindow = rtw)
)

EIC_ntsIUTA[EIC_ntsIUTA$file == 1, ]


ms <- openMSfile(msFilesLowLevel[1])

hd <- header(ms)

class(ms)

class(hd)

View(hd)

#scan 1 and 2
mzR::peaks(ms, c(480,500))

instrumentInfo(ms)
runInfo(ms)
analyzer(ms)



hd2 <- hd[hd$msLevel == 2, ]
i <- which.max(hd2$basePeakIntensity)
hd2[i, ]

pi <- peaks(ms, hd2[i, "seqNum"])
plot(pi, type = "h")

mz <- hd2[i, "basePeakMZ"]
plot(pi, type = "h", xlim = c(mz - 0.5, mz + 0.5), ylim = c(0, 20000))


## Zooming into spectrum 300 (an MS1 spectrum).
j <- 17
pj <- peaks(ms, j)
plot(pj, type = "l")

mz <- hd[j, "basePeakMZ"]
plot(pj, type = "h", xlim = c(mz - 0.5, mz + 0.5))



### with xcmsRaw ------------------------------------------------------------------

library(xcms)

xRaw <- as.list(msFilesLowLevel[1:3])
system.time(xRaw <- lapply(xRaw, function(x) xcmsRaw(x, includeMSn = TRUE)))

system.time(
EIC_xcms <- lapply(xRaw, function(x) {
  eic <- getEIC(x, mzrange = mzr, rtrange = rtr)
  eic <- as.data.table(eic@eic$xcmsRaw)
  return(eic)
})
)


EIC_xcms[[1]]



### with xcmsnExp --------------------------------------------------------------------------------------------

library(xcms)
library(MSnbase)


Exp <- setup@MSnExp

bpis <- chromatogram(Exp, aggregationFun = "max")
plot(bpis)

system.time(chr_raw <- chromatogram(filterFile(Exp, 1), mz = mzr, rt = rtr))

fData(Exp)


Exp[["F1.S0003"]]

head(rtime(Exp))

EIC_msn <- filterRt(Exp, rt = rtr)
EIC_msn <- filterMz(EIC_msn, mzlim = mzr)
EIC_msn <- filterMsLevel(EIC_msn, msLevel. = 1)
system.time(EIC_msn <- chromatogram(EIC_msn))

scans <- fData(EIC_msn)


plot(EIC_msn)


### with patRoon -----

#ms <- patRoon:::loadSpectra(msFilesLowLevel[1])

ms <- as.list(msFilesLowLevel[1:3])
system.time(ms <- lapply(ms, function(x) patRoon:::loadSpectra(x)))
system.time(EIC <- lapply(ms, function(x) patRoon:::loadEICs(x, mzMins = c(mzr[1], mzr[1]), mzMaxs = c(mzr[2], mzr[2]),
                                                             rtMins = c(rtr[1], rtr[1]), rtMaxs = c(rtr[2], rtr[2]))))
EIC[[1]]

ms <- example01@samples$file[1] #path of the msFile
ms <- patRoon:::loadSpectra(ms)
ms <- patRoon:::loadEICs(ms, mzMins = c(233.023, 242.142), mzMaxs = c(233.025, 242.145), rtMins = c(15 * 60, 14.2 * 60), rtMaxs = c(16 * 60, 15.4 * 60))

ms <- patRoon:::loadEICs(ms, mzMins = c(200), mzMaxs = c(300), rtMins = c(800), rtMaxs = c(900))


### raMS -------

library(RaMS)

msdata <- grabMSdata(msFilesLowLevel[19])

saveRDS(msdata, file = paste0(choose.dir(), "/msdata.rds"))

library(DBI)
library(magrittr)

con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
DBI::dbWriteTable(con, "MS1data", msdata$MS1)


library(xml2)
require(XML)

install.packages("xml2")

test <- read_xml(msFilesLowLevel[19]) %>% as.list()

driver_tb = tibble::as_tibble(test) %>% tidyr::unnest_longer('mzML')

test2 <- xml2::xml_find_all(test, xpath = "//mzML")


test2 <-  xmlParse(msFilesLowLevel[19])
test2 <- xmlToList(test2)
test3 <- test2[["mzML"]][["run"]][[1]][[1]][["binaryDataArrayList"]]

