
pathLowLevel <- system.file(package = "ntsIUTA", dir = "extdata")

msFilesLowLevel <- list.files(path = path,
  pattern = ".mzML|.mzXML",
  recursive = TRUE,
  full.names = TRUE,
  no.. = TRUE)

library(mzR)

ms <- openMSfile(msFilesLowLevel[1])

View(ms)

hd <- header(ms)

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
plot(pi, type = "h", xlim = c(mz - 0.5, mz + 0.5), ylim = c(0, 50))


## Zooming into spectrum 300 (an MS1 spectrum).
j <- 17
pj <- peaks(ms, j)
plot(pj, type = "l")

mz <- hd[j, "basePeakMZ"]
plot(pj, type = "h", xlim = c(mz - 0.5, mz + 0.5))




