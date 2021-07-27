

# filterPeaks <- function(peaks = peaks,
#                         fileIndex = NULL,
#                         mz = NULL,
#                         ppm = NULL,
#                         rt = NULL,
#                         rtWindow = NULL,
#                         rtUnit = "sec",
#                         absMinIntensity = NULL,
#                         chromWidthRange = NULL,
#                         negate = FALSE,
#                         dataframe = TRUE) {
#
#   x <- peaks
#
#   x <- x[fileIndex]
#
#   mzr <- mzrBuilder(mz = mz, ppm = ppm)
#
#   rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
#
#   x <- patRoon::filter(obj = x,
#                        absMinIntensity = absMinIntensity,
#                        relMinIntensity = NULL,
#                        retentionRange = rtr,
#                        mzRange = mzr,
#                        mzDefectRange = NULL,
#                        chromWidthRange = chromWidthRange,
#                        negate = negate)
#
#   if (dataframe) x <- patRoon::as.data.frame(x)
#
#   return(x)
#
# }
