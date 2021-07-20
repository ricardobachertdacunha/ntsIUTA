

setGeneric("pullSamples", function(object) standardGeneric("pullSamples"))

setGeneric("pullSamplegroups", function(object) standardGeneric("pullSamplegroups"))
setGeneric("pushSamplegroups", function(object, value) standardGeneric("pushSamplegroups"))
setGeneric("pushSamplegroups<-", function(object, value) standardGeneric("pushSamplegroups<-"))

setGeneric("pullBlanks", function(object) standardGeneric("pullBlanks"))
setGeneric("pushBlanks", function(object, value) standardGeneric("pushBlanks"))
setGeneric("pushBlanks<-", function(object, value) standardGeneric("pushBlanks<-"))

setGeneric("pullQC", function(object) standardGeneric("pullQC"))
setGeneric("pushQC", function(object, value) standardGeneric("pushQC"))
setGeneric("pushQC<-", function(object, value) standardGeneric("pushQC<-"))

setGeneric("pullPeaks", function(obj,
                             fileIndex = NULL, ID = NULL,
                             mz = NULL, ppm = 20,
                             rt = NULL, rtWindow = NULL,
                             rtUnit = "sec") standardGeneric("pullPeaks"))

setGeneric("pullFeatures", function(obj,
                                fileIndex = NULL, ID = NULL,
                                mz = NULL, ppm = 20,
                                rt = NULL, rtWindow = NULL,
                                rtUnit = "sec") standardGeneric("pullFeatures"))

setGeneric("plotRawChrom", function(obj, fileIndex = NULL,
                                    mz = NULL, ppm = 20,
                                    rt = NULL, rtWindow = NULL,
                                    rtUnit = "sec",
                                    msLevel = 1,
                                    type = "tic",
                                    colorBy = "samplegroups",
                                    interactive = FALSE) standardGeneric("plotRawChrom"))

setGeneric("plotCorReplicates", function(obj, binsize = 2) standardGeneric("plotCorReplicates"))

setGeneric("plotFeatures", function(obj, fileIndex = NULL, ID = NULL,
                                    mz = NULL, ppm = 20,
                                    rt = NULL, rtWindow = NULL,
                                    rtUnit = "sec",
                                    msLevel = 1,
                                    colorBy = "features",
                                    interactive = FALSE) standardGeneric("plotFeatures"))
