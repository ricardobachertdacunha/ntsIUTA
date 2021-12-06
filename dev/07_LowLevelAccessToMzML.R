
### Get Paths -----

pathLowLevel <- system.file(package = "ntsIUTA", dir = "extdata")

msFilesLowLevel <- list.files(path = path,
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


library(mzR)
library(data.table)

ms <- as.list(msFilesLowLevel[1:3])
system.time(ms <- lapply(ms, openMSfile))
hd <- lapply(ms, header)
hd <- lapply(hd, as.data.table)

#### get MS1 EIC -----

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


ms <- openMSfile(msFilesLowLevel[1:3])

hd <- header(ms)

class(hd)

View(hd)

#scan 1 and 2
mzR::peaks(ms, 1:2)
instrumentInfo(hd)

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



### with xcms ------------------------------------------------------------------

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





getEIC(setup@MSnExp, mzrange = mzr, rtrange = rtr)








