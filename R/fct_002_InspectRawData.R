

#' @title getRawInfo
#'
#' @description Get raw information for each sample (i.e., file) in an
#' \linkS4class{ntsData} object.
#'
#' @param ojbect An \linkS4class{ntsData} object.
#' @param samples A numeric or character vector with the index or name of
#' samples to collect raw inforamtion, respectively. The default is \code{NULL}
#' which considers all samples in the \code{object}.
#'
#' @return A data table with raw information for each specified file.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom data.table data.table as.data.table copy
#' @importFrom mzR openMSfile header runInfo
#'
getRawInfo <- function(object, samples = NULL) {

  assertClass(object, "ntsData")

  if (is.character(samples)) {
    if (FALSE %in% (samples %in% object@samples$sample)) {
      warning("Given sample name/s not found in the ntsData object.")
      return(rawinfo)
    }
    samples <- which(x@samples$sample %in% sn)
  }

  spt <- samplesTable(object)
  if (!is.null(samples)) spt <- spt[samples]

  fl <- spt$file

  rawinfo <- data.table(
    sample = spt$sample,
    scans = 0,
    centroided = FALSE,
    msLevels = list(rep(0, 2)),
    rtStart = 0,
    rtEnd = 0,
    mzLow = 0,
    mzHigh = 0,
    numberPrecursos = 0,
    collisionEnergy = 0
  )

  for (f in fl) {
    msf <- mzR::openMSfile(f, backend = "pwiz")
    hd <- as.data.table(mzR::header(msf))
    rInfo <- mzR::runInfo(msf)
    spn <- which(spt$file == f)
    rawinfo[spn, scans := rInfo$scanCount]
    rawinfo[spn, msLevels := rInfo$msLevels]
    rawinfo[spn, rtStart := rInfo$dStartTime]
    rawinfo[spn, rtEnd := rInfo$dEndTime]
    rawinfo[spn, mzLow := rInfo$lowMz]
    rawinfo[spn, mzHigh := rInfo$highMz]
    rawinfo[spn, centroided := TRUE %in% unique(hd$centroided)]
    rawinfo[spn, numberPrecursos := length(unique(hd$precursorMZ))]
    ce <- unique(hd$collisionEnergy)
    ce <- ce[!ce %in% NA]
    if (length(ce) > 1) {
      if (which(f == fl) == 1) {
        rawinfo[, collisionEnergy := as.list(collisionEnergy)]
        rawinfo[, collisionEnergy := list(rep(0, length(ce)))]
      }
    }
    rawinfo[spn, collisionEnergy := ce]
  }
  return(copy(rawinfo))
}




#' @title addRawInfo
#'
#' @description Adds raw info for each sample (i.e., file)
#' in the samples slot of an \linkS4class{ntsData} object.
#'
#' @param ojbect An \linkS4class{ntsData} object.
#'
#' @return An \linkS4class{ntsData} object with raw information added
#' per sample in the samples slot.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#'
addRawInfo <- function(object) {
  assertClass(object, "ntsData")
  rinfo <- getRawInfo(object)
  if (TRUE %in% colnames(rinfo)[-1] %in% colnames(object@samples)) {
    object@samples <- object@samples[, colnames(rinfo) := rinfo][]
  } else {
  object@samples <- object@samples[rinfo, on = .(sample = sample)]
  }

  return(object)
}




#' @title openCacheDBScope
#'
#' @description Opens the cache database for storing processing data.
#' The function is loaded from \pkg{patroon}.
#'
#' @importFrom utils getFromNamespace
#'
openCacheDBScope <- utils::getFromNamespace("openCacheDBScope", "patRoon")




#' @title loadSpectra
#'
#' @description Loads spectra data for a given file.
#' The function is imported from the \pkg{patRoon}.
#'
#' @param filepath A charcater string with the complete file path to load the spectra.
#' @param rtRange A vector with length two specifying the time window
#' (i.e., mininum and maximum) to extract the spectra.
#' @param verbose Logical, set to \code{TRUE} to print processing information.
#' @param cacheDB The database for caching information.
#' Use openCacheDBScope() to set the cacheDB and enable integration within \pkg{patRoon}.
#'
#' @importFrom utils getFromNamespace
#'
loadSpectra <- utils::getFromNamespace("loadSpectra", "patRoon")




#' @title loadEICs
#'
#' @description Loads extracted ion chromatograms
#' from pre-loaded spectra data.
#' The function is imported from the \pkg{patRoon}.
#'
#' @param spectra The spectra data to build the EICs.
#' @param rtMins A numeric vector with the minimum retention times of EICs to build. 
#' @param rtMaxs A numeric vector with the maximum retention times of EICs to build.
#' @param mzMins A numeric vector with the minimum mass-to-charge ratio of EICs to build.
#' @param mzMaxs A numeric vector with the maximum mass-to-charge ratio of EICs to build.
#'
#' @importFrom utils getFromNamespace
#'
loadEICs <- utils::getFromNamespace("loadEICs", "patRoon")




#' @title makeTargets
#'
#' @description Helper function to build \emph{m/z} and retention time (in seconds)
#' target pairs for searching within other functions. Each target is composed of an
#' id and mass (\emph{m/z}) and time (seconds) ranges. When mass is defined without
#' time, the entire time range is used.
#'
#' @param mz A vector with target \emph{m/z} values or
#' a two columns data.table or data.frame with minimum and maximum \emph{m/z} values.
#' Alternatively, \emph{m/z} and retention time values can be given as one data.table/data.frame
#' and the deviations given as \code{ppm} and \code{sec} are used to calculate the ranges.
#' The same also works for min and max values of \emph{m/z} and retention time targets.
#' Note that when mass/time ranges are given \code{ppm} and \code{sec} are not used.
#' An id column with target identifiers can be added with \code{mz} as column.
#' @param rt A vector with target retention time values or
#' a two columns data.table/data.frame with minimum and maximum retention time values.
#' @param ppm A numeric vector of length one with the mass deviation, in ppm.
#' @param sec A numeric vector of length one with the time deviation, in seconds.
#'
#' @return A data.table with \emph{m/z} and retention time target pairs.
#'
#' @export
#'
#' @importFrom data.table data.table is.data.table as.data.table
#'
makeTargets <- function(mz, rt, ppm = 20, sec = 60) {

  if (is.null(mz)) {
    mzrts <- data.table(
      id = NA_character_,
      mz = 0,
      rt = 0,
      mzmin = 0,
      mzmax = 0,
      rtmin = 0,
      rtmax = 0
    )
    return(mzrts)

  } else if (length(mz) >= 1 & is.vector(mz)) {
    mzrts <- data.table(
      id = NA_character_,
      mz = mz,
      rt = 0,
      mzmin = 0,
      mzmax = 0,
      rtmin = 0,
      rtmax = 0
    )
    mzrts[, mzmin := mz - ((ppm / 1E6) * mz)]
    mzrts[, mzmax := mz + ((ppm / 1E6) * mz)]

    if (is.vector(rt) & length(rt) == length(mz)) {
      mzrts$rt <- rt
      mzrts[, rtmin := rt - sec]
      mzrts[, rtmax := rt + sec]
    }

    mzrts[, id := paste(
      round(mzmin, digits = 4),
      "-",
      round(mzmax, digits = 4),
      "/", rtmin,
      "-", rtmax,
      sep = ""
    )][]

    return(mzrts)

  } else if (is.data.frame(mz) | is.data.table(mz)) {
    mz <- as.data.table(mz)

    if ("mz" %in% colnames(mz)) {
      mzrts <- data.table(
        id = NA_character_,
        mz = mz$mz,
        rt = 0,
        mzmin = 0,
        mzmax = 0,
        rtmin = 0,
        rtmax = 0
      )
      mzrts[, mzmin := mz - ((ppm / 1E6) * mz)]
      mzrts[, mzmax := mz + ((ppm / 1E6) * mz)]
      if ("rt" %in% colnames(mz)) {
        mzrts$rt <- mz$rt
        mzrts$rtmin <- mz$rt - sec
        mzrts$rtmax <- mz$rt + sec
      }

    } else if ("mzmin" %in% colnames(mz)) {
      mzrts <- data.table(
        id = NA_character_,
        mz = paste(mz$mzmin,
          "-",
          mz$mzmax,
          sep = ""
          ),
        rt = 0,
        mzmin = mz$mzmin,
        mzmax = mz$mzmax,
        rtmin = 0,
        rtmax = 0
      )
      if ("rtmin" %in% colnames(mz)) {
        mzrts$rt <-  paste(mz$rtmin, "-", mz$rtmax, sep = "")
        mzrts$rtmin <- mz$rtmin
        mzrts$rtmax <- mz$rtmax
      }

    } else {
      mzrts <- data.table(
        id = NULL,
        mz = NULL,
        rt = NULL,
        mzmin = NULL,
        mzmax = NULL,
        rtmin = NULL,
        rtmax = NULL
      )

      return(mzrts)
    }

    if (is.data.frame(rt) | is.data.table(rt)) {
      rt <- as.data.table(rt)

      if ("rt" %in% colnames(rt) & nrow(rt) == nrow(mz)) {
        mzrts$rt <- rt$rt
        mzrts$rtmin <- rt$rt - sec
        mzrts$rtmax <- rt$rt + sec
      }

      if ("rtmin" %in% colnames(rt) & nrow(rt) == nrow(mz)) {
        mzrts$rt <-  paste(rt$rtmin, "-", rt$rtmax, sep = "")
        mzrts$rtmin <- rt$rtmin
        mzrts$rtmax <- rt$rtmax
      }
    }

    if ("id" %in% colnames(mz)) {
        mzrts$id <- mz$id
    } else {
      mzrts[, id := paste(
        round(mzmin, digits = 4),
        "-",
        round(mzmax, digits = 4),
        "/",
        rtmin,
        "-",
        rtmax,
        sep = ""
      )][]
    }
    return(mzrts)

  } else {
    mzrts <- data.table(
      id = NULL,
      mz = NULL,
      rt = NULL,
      mzmin = NULL,
      mzmax = NULL,
      rtmin = NULL,
      rtmax = NULL
    )

    return(mzrts)
  }
}




#' @title extractEICs
#' @description Extracts MS1 level ion chromatograms (EICs) from raw data of a specified \emph{m/z}
#' and retention time pair/s.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples A numeric or character vector with the index or names
#' of the sample/s to extract the data, respectively.
#' When \code{NULL}, all the samples in the \code{object} are used.
#' @param mz A numeric vector with target \emph{m/z} values to extract EICs.
#' Alternatively, a two columns data.table or data.frame with
#' minimum and maximum (in this order) \emph{m/z} values can be used instead of exact \emph{m/z}.
#' Note that for the latter, the argumment \code{ppm} is ignored.
#' @param ppm A numeric vector of length one with the mass deviation, in ppm, to extract the data for the EIC.
#' @param rt A numeric vector with retention time values, in seconds, to extract EICs.
#' When \code{NULL}, the complete retention time window in the sample (i.e., file) is used.
#' Alternatively, a two columns data.table or data.frame with
#' minimum and maximum (in this order) retention time values can be used instead of exact retention times.
#' Note that for the latter, the argumment \code{sec} is ignored.
#' Note that the legnth of the \code{rt} vector of number of rows of the \code{rt} table
#' should be equal to the length or number of rows of \code{mz}.
#' @param sec A numeric vector of length one with the time deviation, in seconds, to extract the data for the EIC.
#'
#' @return A \code{data.table} with the columns
#' \code{sample}, \code{replicate}, \code{id}, \code{rt}, and \code{intensity}
#' representing the sample index (i.e., file), the sample replicate name,
#' the EIC target id, the retention time and the sum of the intensities
#' for the collected \emph{m/z} in each spectrum (i.e., MS scan), respectively.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom data.table rbindlist setnames setcolorder
#'
extractEICs <- function(object = NULL,
                        samples = NULL,
                        mz = NULL, ppm = 20,
                        rt = NULL, sec = 60) {

  assertClass(object, "ntsData")

  if (is.character(samples)) {
    if (FALSE %in% (samples %in% object@samples$sample)) {
      warning("Given sample names not found in the ntsData object!")
      return(data.table())
    }
    samples <- which(object@samples$sample %in% samples)
  }

  fls <- filePaths(object)

  if (!is.null(samples)) fls <- fls[samples]

  targets <- makeTargets(mz, rt, ppm, sec)

  eicList <- list()

  eicList <- lapply(fls, function(x, targets, spt) {

    rtRange <- c(min(targets$rtmin), max(targets$rtmax))
    if (rtRange[1] == 0 & rtRange[2] == 0) rtRange <- NULL

    spectra <- loadSpectra(
      x,
      rtRange = NULL,
      verbose = TRUE,
      cacheDB = openCacheDBScope()
    )

    targets[mz == 0 & rt == 0, id := "tic"]
    targets[mz == 0 & rt == 0, mz := NA]
    targets[mz == 0 & rt == 0, rt := NA]
    targets[rtmin == 0, rtmin := min(spectra$header$retentionTime, na.rm = TRUE)]
    targets[rtmax == 0, rtmax := max(spectra$header$retentionTime, na.rm = TRUE)]
    targets[mzmin == 0, mzmin := min(spectra$header$lowMZ, na.rm = TRUE)]
    targets[mzmax == 0, mzmax := max(spectra$header$highMZ, na.rm = TRUE)]

    eic <- loadEICs(
      spectra = spectra,
      mzMins = targets[["mzmin"]],
      mzMaxs = targets[["mzmax"]],
      rtMins = targets[["rtmin"]],
      rtMaxs = targets[["rtmax"]]
    )

    names(eic) <- targets[["id"]]
    eic <- rbindlist(eic, idcol = "target")
    setnames(eic, c("id", "rt", "intensity"))
    eic[, sample := spt[file == x, sample]]
    eic[, replicate := spt[file == x, replicate]][]

    return(eic)
  }, targets = targets, spt = samplesTable(object)[file %in% fls, ])

  eics <- rbindlist(eicList)
  setcolorder(eics, c("sample", "replicate", "id", "rt", "intensity"))

  return(eics)
}




#' @title extractXICs
#' @description Extracts three dimensional MS1 level ion chromatograms (XICs)
#' from raw data of a specified \emph{m/z} and retention time pair/s.
#' The XIC extraction is slower than the \link{extractEICs}, therefore it should
#' be used for limited/narrow mass and time deviations.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples A numeric or character vector with the index or names
#' of the sample/s to extract the data, respectively.
#' When \code{NULL}, all the samples in the \code{object} are used.
#' @param mz A numeric vector with target \emph{m/z} values to extract XICs.
#' Alternatively, a two columns data.table or data.frame with
#' minimum and maximum (in this order) \emph{m/z} values can be used instead of exact \emph{m/z}.
#' Note that for the latter, the argumment \code{ppm} is ignored.
#' @param ppm A numeric vector of length one with the mass deviation, in ppm, to extract the data for the XIC.
#' @param rt A numeric vector with retention time values, in seconds, to extract XICs.
#' When \code{NULL}, the complete retention time window in the sample (i.e., file) is used.
#' Alternatively, a two columns data.table or data.frame with
#' minimum and maximum (in this order) retention time values can be used instead of exact retention times.
#' Note that for the latter, the argumment \code{sec} is ignored.
#' Note that the legnth of the \code{rt} vector of number of rows of the \code{rt} table
#' should be equal to the length or number of rows of \code{mz}.
#' @param sec A numeric vector of length one with the time deviation, in seconds, to extract the data for the XIC.
#'
#' @return A \code{data.table} with the columns
#' \code{sample}, \code{replicate}, \code{id}, \code{mz}, \code{rt}, and \code{intensity}
#' representing the sample index (i.e., file), the sample replicate name,
#' the XIC target id, the \emph{m/z}, the retention time and the intensity
#' of each \emph{m/z} and retention time pair (e.g., centroid), respectively.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom data.table rbindlist setnames setcolorder
#'
extractXICs <- function(object = NULL,
                        samples = NULL,
                        mz = NULL, ppm = 20,
                        rt = NULL, sec = 60) {

  assertClass(object, "ntsData")

  if (is.character(samples)) {
    if (FALSE %in% (samples %in% object@samples$sample)) {
      warning("Given sample names not found in the ntsData object!")
      return(data.table())
    }
    samples <- which(object@samples$sample %in% samples)
  }

  fls <- filePaths(object)

  if (!is.null(samples)) fls <- fls[samples]

  targets <- makeTargets(mz, rt, ppm, sec)

  xicList <- list()

  xicList <- lapply(fls, function(x, targets, spt) {

    rtRange <- c(min(targets$rtmin) * 0.7, max(targets$rtmax) * 1.3)
    if (rtRange[1] == 0 & rtRange[2] == 0) rtRange <- NULL

    spectra <- loadSpectra(
      x,
      rtRange = rtRange,
      verbose = TRUE,
      cacheDB = openCacheDBScope()
    )

    hd <- as.data.table(spectra$header)
    hd <- hd[retentionTime >= rtRange[1] & retentionTime <= rtRange[2], ]
    hd[, newSeqNum := seq_len(length(spectra$spectra))]

    if (length(spectra$spectra) != nrow(hd)) {
      warning("Size of header table is different than the number of spectra.")
    }

    xics <- list()

    for (i in seq_len(nrow(targets))) {
      xics[[targets$id[i]]] <- hd[
        retentionTime >= targets$rtmin[i] &
        retentionTime <= targets$rtmax[i],
      ][, .(newSeqNum, seqNum, msLevel, retentionTime)]

      xics[[targets$id[i]]] <- xics[[targets$id[i]]][msLevel == 1, ]

      xics[[targets$id[i]]] <- lapply(seq_len(nrow(xics[[targets$id[i]]])), function(x, scans, tg) {
        xic <- spectra$spectra[[scans$newSeqNum[x]]]
        xic <- as.data.table(xic)
        setnames(xic, c("mz", "intensity"))
        xic <- xic[mz >= tg$mzmin & mz <= tg$mzmax]
        xic[, rt := scans$retentionTime[x]]
        return(xic)
      }, scans = xics[[targets$id[i]]], tg = targets[i])

      xics[[targets$id[i]]] <- rbindlist(xics[[targets$id[i]]])
      xics[[targets$id[i]]][, id := targets$id[i]]
      xics[[targets$id[i]]][, mz_id := targets$mz[i]]
      xics[[targets$id[i]]][, rt_id := targets$rt[i]]
    }

    xics <- rbindlist(xics)
    xics[, sample := spt[file == x, sample]]
    xics[, replicate := spt[file == x, replicate]]

    return(xics)
  }, targets = targets, spt = samplesTable(object)[file %in% fls, ])

  xics <- rbindlist(xicList)
  setcolorder(xics, c("sample", "replicate", "id", "mz_id", "rt_id", "mz", "rt", "intensity"))

  return(xics)
}




#' @title clusterMsn
#'
#' @description Function to cluster MSn data.
#'
#' @return A data table with clustered MSn data for given targets.
#' 
#' @importFrom data.table copy setcolorder setorder
#' @importFrom fastcluster hclust
#'
clusterMSn <- function(
  ids,
  msnList,
  clusteringMethod,
  clusteringUnit,
  clusteringWindow,
  mergeCEs,
  mergeBy,
  targets) {

  msnList <- lapply(ids, function(
    x,
    msnList,
    clusteringMethod,
    clusteringUnit,
    clusteringWindow,
    mergeCEs,
    mergeBy,
    targets
  ) {

    t <- msnList[id == x, ]
    idf <- targets[id == x, ]

    if (nrow(t) > 2) {

      if (clusteringMethod == "distance") {
        setorder(t, mz)
        mzMat <- abs(diff(t$mz))
        if (clusteringUnit == "ppm") {
          mzMat <- (mzMat / t$mz[-1]) * 1E6
        }
        t[, cluster := 1 + c(0, cumsum(mzMat > clusteringWindow))]

      } else {
        mzMat <- dist(t$mz, method = clusteringMethod)
        if (clusteringUnit == "ppm") {
          mzMat <- as.data.table(as.matrix(mzMat))
          mzMat <- mzMat[, lapply(.SD, function(x, dt) x / t$mz * 1E6, dt = t), .SDcols = colnames(mzMat)]
          mzMat <- as.dist(mzMat)
        }
        hc <- fastcluster::hclust(mzMat, method = "complete")
        t[, cluster := cutree(hc, h = clusteringWindow)]

      }

      if ("CE" %in% colnames(t)) {
        if (any(t[, .(dup = anyDuplicated(seqNum)), key = c("cluster", "sample")][["dup"]] > 0)) {
          warning(paste0("MSMS traces from the same spectrum were merged for ", idf$id, "\n"))
        }
      }

      if (is.null(mergeBy)) mergeBy <- "id"

      if (mergeCEs) {

        if (mergeBy == "replicates" & "CE" %in% colnames(t)) {
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum)),
            CE = I(list(unique(CE))),
            preMZ = mean(preMZ)
            ), by = list(replicate, cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, replicate, mz)
          setcolorder(t, c("replicate", "id", "mz", "intensity", "CE", "preMZ", "precursor"))

        } else if (mergeBy != "id"  & "CE" %in% colnames(t)) { #mergeBy samples when mergeBy is not null
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum)),
            CE = I(list(unique(CE))),
            preMZ = mean(preMZ),
            replicate = unique(replicate)
            ), by = list(sample, cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, sample, mz)
          setcolorder(t, c("sample", "replicate", "id", "mz", "intensity", "CE", "preMZ", "precursor"))

        } else if ("CE" %in% colnames(t)) { #when NULL do not merge by samples
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum)),
            CE = I(list(unique(CE))),
            preMZ = mean(preMZ)
            ), by = list(cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, mz)
          setcolorder(t, c("id", "mz", "intensity", "CE", "preMZ", "precursor"))

        } else {
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum))
            ), by = list(cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, mz)
          setcolorder(t, c("id", "mz", "intensity", "precursor"))
        }

      } else {

        if (mergeBy == "replicates") {
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum)),
            preMZ = mean(preMZ)
            ), by = list(CE, replicate, cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, replicate, CE, mz)
          setcolorder(t, c("replicate", "id", "mz", "intensity", "CE", "preMZ", "precursor"))

        } else if (mergeBy != "id") { #mergeBy samples when mergeBy is not null
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum)),
            preMZ = mean(preMZ),
            replicate = unique(replicate)
            ), by = list(CE, sample, cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, sample, CE, mz)
          setcolorder(t, c("sample", "replicate", "id", "mz", "intensity", "CE", "preMZ", "precursor"))

        } else { #when NULL do not merge by samples
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(seqNum)),
            preMZ = mean(preMZ)
            ), by = list(CE, cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t[, id := idf$id][]
          setorder(t, CE, mz)
          setcolorder(t, c("id", "mz", "intensity", "CE", "preMZ", "precursor"))
        }
      }
    }

    return(t)
  },
  clusteringMethod = clusteringMethod,
  clusteringUnit = clusteringUnit,
  clusteringWindow = clusteringWindow,
  mergeCEs = mergeCEs,
  mergeBy = mergeBy,
  targets = targets,
  msnList = msnList)

  msnList <- rbindlist(msnList)

  return(msnList)
}




#' @title clusterMSnToPatRoon
#'
#' @description Function to cluster spectra and convert to
#' an \linkS4class{MSPeakLists} object.
#'
#' @return A \linkS4class{MSPeakLists} object with clustered data.
#'
#' @importFrom data.table copy setnames setcolorder setorder
#' @importClassesFrom patRoon MSPeakLists
#' @importFrom fastcluster hclust
#'
clusterMSnToPatRoon <- function(
  cl_plists,
  mlists,
  targets,
  clusteringMethod,
  clusteringUnit,
  clusteringWindow,
  minIntensityPost) {

  av_plists <- list()

  for (i in names(cl_plists)) {
    for (f in names(cl_plists[[i]])) {

      av_ms <- clusterMSn(
        ids = f,
        msnList = cl_plists[[i]][[f]]$MS,
        clusteringMethod,
        clusteringUnit,
        clusteringWindow,
        mergeCEs = TRUE,
        mergeBy = "sample",
        targets
      )

      av_ms <- av_ms[intensity > minIntensityPost, ]
      av_ms[, sample := i]

      av_msms <- clusterMSn(
        ids = f,
        msnList = cl_plists[[i]][[f]]$MSMS,
        clusteringMethod,
        clusteringUnit,
        clusteringWindow,
        mergeCEs = TRUE,
        mergeBy = "sample",
        targets
      )

      av_msms <- av_msms[intensity > minIntensityPost, ]

      av_plists[[f]][["MS"]] <- c(av_plists[[f]][["MS"]], list(av_ms))
      av_plists[[f]][["MSMS"]] <- c(av_plists[[f]][["MSMS"]], list(av_msms))

      cl_plists[[i]][[f]]$MS <- copy(av_ms)
      cl_plists[[i]][[f]]$MSMS <- copy(av_msms)

      setnames(cl_plists[[i]][[f]]$MS, "id", "ID")
      cl_plists[[i]][[f]]$MS[, sample := NULL]
      cl_plists[[i]][[f]]$MS[, ID := seq_len(nrow(cl_plists[[i]][[f]]$MS))]

      setnames(cl_plists[[i]][[f]]$MSMS, "id", "ID")
      cl_plists[[i]][[f]]$MSMS[, sample := NULL]
      cl_plists[[i]][[f]]$MSMS[, replicate := NULL]
      cl_plists[[i]][[f]]$MSMS[, CE := NULL]
      cl_plists[[i]][[f]]$MSMS[, preMZ := NULL]
      cl_plists[[i]][[f]]$MSMS[, ID := seq_len(nrow(cl_plists[[i]][[f]]$MSMS))]
    }
  }

  for (f in names(av_plists)) {
    for (lv in names(av_plists[[f]])) {

      t <- copy(rbindlist(av_plists[[f]][[lv]], fill = TRUE))
      idf <- targets[id == f, ]

      if (nrow(t) > 2) {

        if (clusteringMethod == "distance") {
          setorder(t, mz)
          mzMat <- abs(diff(t$mz))
          if (clusteringUnit == "ppm") {
            mzMat <- (mzMat / t$mz[-1]) * 1E6
          }
          t[, cluster := 1 + c(0, cumsum(mzMat > clusteringWindow))]

        } else {
          mzMat <- dist(t$mz, method = clusteringMethod)
          if (clusteringUnit == "ppm") {
            mzMat <- as.data.table(as.matrix(mzMat))
            mzMat <- mzMat[, lapply(.SD, function(x, dt) x / t$mz * 1E6, dt = t), .SDcols = colnames(mzMat)]
            mzMat <- as.dist(mzMat)
          }
          hc <- fastcluster::hclust(mzMat, method = "complete")
          t[, cluster := cutree(hc, h = clusteringWindow)]

        }

        if ("CE" %in% colnames(t)) {
          if (any(t[, .(dup = anyDuplicated(sample)), key = list(cluster)][["dup"]] > 0)) {
            warning(paste0("MSMS traces from the same sample were merged during averaging for ", idf$id, "\n"))
          }

          t <- t[, .(
          mz = mean(mz),
          intensity = sum(intensity) / length(unique(sample))
          ), by = list(cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t <- t[intensity > minIntensityPost, ]
          setorder(t, mz)
          t[, ID := seq_len(nrow(t))][]
          setcolorder(t, c("ID", "mz", "intensity", "precursor"))

        } else {
          t <- t[, .(
            mz = mean(mz),
            intensity = sum(intensity) / length(unique(sample))
          ), by = list(cluster)
          ][, cluster := NULL]

          t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
          t[is.na(precursor), precursor := FALSE]
          t <- t[intensity > minIntensityPost, ]
          setorder(t, mz)
          t[, ID := seq_len(nrow(t))]
          setcolorder(t, c("ID", "mz", "intensity", "precursor"))
        }
      }

      av_plists[[f]][[lv]] <- copy(t)
    }
  }

  final_plists <- new("MSPeakLists",
    peakLists = cl_plists,
    metadata = mlists,
    averagedPeakLists = av_plists,
    avgPeakListArgs = list(),
    origFGNames = names(av_plists),
    algorithm = "mzr"
  )

  final_plists@averagedPeakLists <- av_plists

  return(final_plists)
}




#' @title extractMSn
#' @description Extracts MSn spectra from defined isolated targets
#' defined by \emph{m/z} and retention time, including the respective deviations.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples A numeric or character vector with the index or names
#' of the sample/s to extract the data, respectively.
#' When \code{NULL}, all the samples in the \code{object} are used.
#' @param level A numeric vector with length 1 to defined the MS level.
#' Only level 2, corresponding to MS/MS was tested.
#' @param mz A numeric vector with isolation target \emph{m/z} values to extract MSn.
#' Alternatively, a two columns data.table or data.frame with
#' minimum and maximum (in this order) \emph{m/z} values can be used instead of exact \emph{m/z}.
#' Note that for the latter, the argumment \code{ppm} is ignored.
#' @param ppm A numeric vector of length one with the mass deviation, in ppm, to screen for isolation targets.
#' @param rt A numeric vector with retention time values, in seconds, to extract MSn.
#' When \code{NULL}, the complete retention time window in the sample (i.e., file) is used.
#' Alternatively, a two columns data.table or data.frame with
#' minimum and maximum (in this order) retention time values can be used instead of exact retention times.
#' Note that for the latter, the argumment \code{sec} is ignored.
#' Note that the legnth of the \code{rt} vector of number of rows of the \code{rt} table
#' should be equal to the length or number of rows of \code{mz}.
#' @param sec A numeric vector of length one with the time deviation, in seconds, to screen for isolation targets.
#' @param isolationTimeWindow A numeric vector of length one with the time deviation of the isolation window, in seconds.
#' @param isolationMassWindow A numeric vector of length one with the mass isolation window, in Da.
#' @param clusteringMethod A character vector with the method for clustering.
#' Possible values are \emph{euclidean} (the default) or \emph{distance}.
#' @param clusteringUnit A character vector specifying the clustering unit.
#' Possibel values are \emph{mz} (the default) or \emph{ppm}.
#' @param clusteringWindow A length one numeric vector with the mass deviation
#' for clustering \emph{m/z} values across different fragmentation spectra.
#' @param minIntensityPre A length one numeric vector with the minimum intensity of peaks
#' to be applied before clustering \emph{m/z} values across spectra.
#' @param minIntensityPost A length one numeric vector with the minimum intensity of peaks
#' to be applied after clustering and averaging of the \emph{m/z} and intensity values.
#' @param asPatRoon Logical, set to \code{TRUE} for return an \linkS4class{MSPeakLists} object.
#' @param mergeCEs Logical, set to \code{TRUE} to cluster different collision energies.
#' When FALSE, the spectra from different collision energies won't be clustered.
#' @param mergeBy A length one character vector with either "samples" (the default) or "replicates"
#' to merge the MS2 data for different samples or replicates, respectively.
#' When \code{NULL}, the spectra for each target is given per sample.
#'
#' @return A \code{data.table} with the columns
#' \code{sample}, \code{replicate}, \code{id}, \code{CE}, \code{mz}, \code{intensity} and \code{precursor}
#' representing the sample name (i.e., file), the sample replicate name,
#' the isolation target id, the collision energy applied, the \emph{m/z}, the intensity and the presence of the precursor
#' (i.e., \emph{m/z} matching the isolation target), respectively.
#'
#' @export
#'
#' @importFrom fastcluster hclust
#' @importFrom checkmate assertClass
#' @importFrom data.table rbindlist setnames setorder as.data.table setcolorder data.table copy
#'
extractMSn <- function(
  object = NULL,
  samples = NULL,
  level = 2,
  mz = NULL, ppm = 20,
  rt = NULL, sec = 60,
  isolationTimeWindow = 10,
  isolationMassWindow = 1.3,
  clusteringMethod = "distance",
  clusteringUnit = "ppm",
  clusteringWindow = 25,
  minIntensityPre = 200,
  minIntensityPost = 200,
  asPatRoon = TRUE,
  mergeCEs = TRUE,
  mergeBy = "samples") {

  checkmate::assertClass(object, "ntsData")

  if (is.character(samples)) {
    samples <- which(object@samples$sample %in% samples)
  }

  fls <- filePaths(object)

  if (!is.null(samples)) fls <- fls[samples]

  targets <- makeTargets(mz, rt, ppm, sec)

  if (length(targets$mz) == 1) {
    if (targets$mz == 0) {
      targets <- features(object)
      targets <- targets[, .(id, mz, rt, mzmin, mzmax, rtmin, rtmax)]
    }
  }

  if (clusteringUnit == "ppm") {
    targets$mzmin <- targets$mz - (clusteringWindow / 1E6 * targets$mz)
    targets$mzmax <- targets$mz + (clusteringWindow / 1E6 * targets$mz)
  } else {
    targets$mzmin <- targets$mz - clusteringWindow
    targets$mzmax <- targets$mz + clusteringWindow
  }


  mlists <- list()

  dummy <- data.table(
    mz = numeric(),
    intensity = numeric(),
    seqNum = numeric(),
    CE = numeric(),
    preMZ = numeric()
  )

  plists <- lapply(fls, function(
    x,
    targets,
    spt,
    level,
    isolationTimeWindow,
    isolationWindow,
    minIntensityPre,
    dummy) {

    rtRange <- c(min(targets$rtmin) * 0.7, max(targets$rtmax) * 1.3)
    if (rtRange[1] == 0 & rtRange[2] == 0) rtRange <- NULL

    spectra <- loadSpectra(
      x,
      rtRange = rtRange,
      verbose = TRUE,
      cacheDB = openCacheDBScope()
    )

    hd <- as.data.table(spectra$header)
    hd <- hd[retentionTime >= rtRange[1] & retentionTime <= rtRange[2], ]
    hd[, newSeqNum := seq_len(length(spectra$spectra))]

    if (length(spectra$spectra) != nrow(hd)) {
      warning("Size of header table is different than the number of spectra.")
    }

    pHolder <- list()
    mHolder <- list()

    for (i in seq_len(nrow(targets))) {

      idf <- targets$id[i]

      mHolder[[idf]] <- list()
      pHolder[[idf]] <- list()

      hd_msms <- copy(hd[msLevel == level, ])

      hd_msms <- hd_msms[
        retentionTime >= (targets$rt[i] - isolationTimeWindow) &
        retentionTime <= (targets$rt[i] + isolationTimeWindow) &
        precursorMZ >= (targets$mz[i] - isolationMassWindow) &
        precursorMZ <= (targets$mz[i] + isolationMassWindow),
      ]

      hd_ms <- copy(hd[msLevel == (level - 1), ])
      hd_ms <- hd_ms[acquisitionNum %in% unique(hd_msms$precursorScanNum), ]

      mHolder[[idf]][["MS"]] <- hd_ms
      mHolder[[idf]][["MSMS"]] <- hd_msms

      msms <- lapply(seq_len(nrow(hd_msms)), function(x_msms, spectra, scans) {
        prd <- as.data.table(spectra$spectra[[scans$newSeqNum[x_msms]]])
        setnames(prd, c("mz", "intensity"))
        prd[, seqNum := scans$seqNum[x_msms]]
        prd[, CE := scans$collisionEnergy[x_msms]]
        prd[, preMZ := scans$precursorMZ[x_msms]]
        return(prd)
      }, scans = hd_msms, spectra = spectra)

      msms <- rbindlist(c(msms, list(dummy)))
      msms <- msms[intensity >= minIntensityPre, ]

      ms <- lapply(seq_len(nrow(hd_ms)), function(x_ms, spectra, scans) {
        prd <- as.data.table(spectra$spectra[[scans$newSeqNum[x_ms]]])
        setnames(prd, c("mz", "intensity"))
        prd[, seqNum := scans$seqNum[x_ms]]
        return(prd)
      }, scans = hd_ms, spectra = spectra)

      ms <- rbindlist(c(ms, list(dummy[, .(mz, intensity, seqNum)])))
      ms <- ms[intensity >= minIntensityPre, ]

      msms[, id := idf]
      msms[, sample := spt[file == x, sample]]
      msms[, replicate := spt[file == x, replicate]]

      ms[, id := idf]
      ms[, sample := spt[file == x, sample]]
      ms[, replicate := spt[file == x, replicate]]

      pHolder[[idf]][["MS"]] <- ms
      pHolder[[idf]][["MSMS"]] <- msms
    }

    mlists[[spt[file == x, sample]]] <<- mHolder

    return(pHolder)

  },
  targets = targets,
  spt = samplesTable(object)[file %in% fls, ],
  level = level,
  isolationTimeWindow = isolationTimeWindow,
  isolationWindow = isolationWindow,
  minIntensityPre = minIntensityPre,
  dummy = dummy
  )

  names(plists) <- samplesTable(object)[file %in% fls, sample]

  cat("Clustering spectra... \n")

  if (asPatRoon) {

    cl_plists <- copy(plists)

    return(
      clusterMSnToPatRoon(
        cl_plists,
        mlists,
        targets,
        clusteringMethod,
        clusteringUnit,
        clusteringWindow,
        minIntensityPost
      )
    )

  } else {
    msnList <- lapply(plists, function(x) lapply(x, function(y) y[which(names(y) == "MSMS")]))
    msnList <- lapply(msnList, function(x) rbindlist(lapply(x, function(y) rbindlist(y, fill = TRUE)), fill = TRUE))
    msnList <- rbindlist(msnList, fill = TRUE)
    ids <- unique(msnList$id)

    if (length(ids) > 0) {
      msnList <- clusterMSn(
        ids,
        msnList,
        clusteringMethod,
        clusteringUnit,
        clusteringWindow,
        mergeCEs,
        mergeBy,
        targets
      )
    }

    msnList <- msnList[intensity >= minIntensityPost, ]

    return(msnList)
  }
}
