

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
    mzrts <- data.table(id = NA_character_, mz = 0, rt = 0, mzmin = 0, mzmax = 0, rtmin = 0, rtmax = 0)
    return(mzrts)

  } else if (length(mz) >= 1 & is.vector(mz)) {
    mzrts <- data.table(id = NA_character_, mz = mz, rt = 0, mzmin = 0, mzmax = 0, rtmin = 0, rtmax = 0)
    mzrts[, mzmin := mz - ((ppm / 1E6) * mz)]
    mzrts[, mzmax := mz + ((ppm / 1E6) * mz)]

    if (is.vector(rt) & length(rt) == length(mz)) {
      mzrts$rt <- rt
      mzrts[, rtmin := rt - sec]
      mzrts[, rtmax := rt + sec]
    }

    mzrts[, id := paste(round(mzmin, digits = 4), "-", round(mzmax, digits = 4), "/", rtmin, "-", rtmax, sep = "")][]
    return(mzrts)

  } else if (is.data.frame(mz) | is.data.table(mz)) {
    mz <- as.data.table(mz)

    if ("mz" %in% colnames(mz)) {
      mzrts <- data.table(id = NA_character_, mz = mz$mz, rt = 0, mzmin = 0, mzmax = 0, rtmin = 0, rtmax = 0)
      mzrts[, mzmin := mz - ((ppm / 1E6) * mz)]
      mzrts[, mzmax := mz + ((ppm / 1E6) * mz)]
      if ("rt" %in% colnames(mz)) {
        mzrts$rt <- mz$rt
        mzrts$rtmin <- mz$rt - sec
        mzrts$rtmax <- mz$rt + sec
      }

    } else if ("mzmin" %in% colnames(mz)) {
      mzrts <- data.table(
        id = NA_character_, mz = paste(mz$mzmin, "-", mz$mzmax, sep = ""), rt = 0,
        mzmin = mz$mzmin, mzmax = mz$mzmax, rtmin = 0, rtmax = 0
      )
      if ("rtmin" %in% colnames(mz)) {
        mzrts$rt <-  paste(mz$rtmin, "-", mz$rtmax, sep = "")
        mzrts$rtmin <- mz$rtmin
        mzrts$rtmax <- mz$rtmax
      }

    } else {
      mzrts <- data.table(id = NULL, mz = NULL, rt = NULL, mzmin = NULL, mzmax = NULL, rtmin = NULL, rtmax = NULL)
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
      mzrts[, id := paste(round(mzmin, digits = 4), "-", round(mzmax, digits = 4), "/", rtmin, "-", rtmax, sep = "")][]
    }
    return(mzrts)

  } else {
    mzrts <- data.table(id = NULL, mz = NULL, rt = NULL, mzmin = NULL, mzmax = NULL, rtmin = NULL, rtmax = NULL)
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

    rtRange <- c(min(targets$rtmin), max(targets$rtmax))
    if (rtRange[1] == 0 & rtRange[2] == 0) rtRange <- NULL

    spectra <- loadSpectra(
      x,
      rtRange = NULL,
      verbose = TRUE,
      cacheDB = openCacheDBScope()
    )

    xics <- list()

    for (i in seq_len(nrow(targets))) {
      xics[[targets$id[i]]] <- spectra$header[
        retentionTime >= targets$rtmin[i] &
        retentionTime <= targets$rtmax[i],
      ][, .(seqNum, msLevel, retentionTime)]

      xics[[targets$id[i]]] <- xics[[targets$id[i]]][msLevel == 1, ]

      xics[[targets$id[i]]] <- lapply(seq_len(nrow(xics[[targets$id[i]]])), function(x, scans, tg) {
        xic <- spectra$spectra[[scans$seqNum[x]]]
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
#' @param mergeCEs Logical, set to TRUE to cluster different collision energies.
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
#' @importFrom data.table rbindlist setnames setorder as.data.table setcolorder data.table
#'
extractMSn <- function(object = NULL,
                       samples = NULL,
                       level = 2,
                       mz = NULL, ppm = 20,
                       rt = NULL, sec = 60,
                       clusteringMethod = "euclidean",
                       clusteringUnit = "mz",
                       clusteringWindow = 0.008,
                       minIntensityPre = 250,
                       minIntensityPost = 100,
                       mergeCEs = FALSE,
                       mergeBy = "samples") {

  assertClass(object, "ntsData")

  if (is.character(samples)) {
    samples <- which(object@samples$sample %in% samples)
  }

  fls <- filePaths(object)

  if (!is.null(samples)) fls <- fls[samples]

  targets <- makeTargets(mz, rt, ppm, sec)

  msnList <- list()

  dummy <- data.table(
    mz = numeric(),
    intensity = numeric(),
    seqNum = numeric(),
    CE = numeric(),
    preMZ = numeric()
  )

  msnList <- lapply(fls, function(x, targets, spt, level, minIntensityPre, dummy) {

    rtRange <- c(min(targets$rtmin), max(targets$rtmax))
    if (rtRange[1] == 0 & rtRange[2] == 0) rtRange <- NULL

    spectra <- loadSpectra(
      x,
      rtRange = NULL,
      verbose = TRUE,
      cacheDB = openCacheDBScope()
    )

    msn <- list()

    for (i in seq_len(nrow(targets))) {

      msnSpec <- as.data.table(spectra$header[msLevel == level]) #redundant, but maybe needed for levels higher than 2

      msnSpec <- msnSpec[
        retentionTime >= targets$rtmin[i] &
        retentionTime <= targets$rtmax[i] &
        precursorMZ >= targets$mzmin[i] &
        precursorMZ <= targets$mzmax[i],
      ]

      msnSpec <- lapply(seq_len(nrow(msnSpec)), function(x, spectra, scans) {
        prd <- as.data.table(spectra$spectra[[scans$seqNum[x]]])
        setnames(prd, c("mz", "intensity"))
        prd[, seqNum := scans$seqNum[x]] #not needed?
        prd[, CE := scans$collisionEnergy[x]]
        prd[, preMZ := scans$precursorMZ[x]]
        return(prd)
      }, scans = msnSpec, spectra = spectra)

      msnSpec <- rbindlist(c(msnSpec, list(dummy)))

      msnSpec <- msnSpec[intensity >= minIntensityPre, ]

      msn[[targets$id[i]]] <- msnSpec
    }

    msn <- rbindlist(msn, idcol = "id")

    msn[, sample := spt[file == x, sample]]

    msn[, replicate := spt[file == x, replicate]]

    return(msn)

  },
  targets = targets,
  spt = samplesTable(object)[file %in% fls, ],
  level = level,
  minIntensityPre = minIntensityPre,
  dummy = dummy
  )

  msnList <- rbindlist(msnList, fill = TRUE)

  ids <- unique(msnList$id)

  if (length(ids) > 0) {

    msnList <- lapply(ids, function(x, msnList, clusteringMethod, clusteringUnit, clusteringWindow, mergeCEs, mergeBy, targets) {

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

        if (any(t[, .(dup = anyDuplicated(seqNum)), key = c("cluster", "sample")][["dup"]] > 0)) {
          warning("During spectral averaging multiple masses
          from the same spectrum were clustered, consider tweaking clusterMzWindow!\n")
        }

        if (mergeCEs) {

          if (mergeBy == "replicates") {
            t <- t[, .(
              mz = mean(mz),
              intensity = sum(intensity) / length(unique(seqNum)),
              CE = paste0(unique(CE), collapse = "/"),
              preMZ = mean(preMZ)
              ), by = list(replicate, cluster)
            ][, cluster := NULL]

            t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
            t[is.na(precursor), precursor := FALSE]
            t[, id := idf$id][]
            setorder(t, replicate, mz)
            setcolorder(t, c("replicate", "id", "mz", "intensity", "CE", "preMZ", "precursor"))

          } else if (!is.null(mergeBy)) { #mergeBy samples when mergeBy is not null
            t <- t[, .(
              mz = mean(mz),
              intensity = sum(intensity) / length(unique(seqNum)),
              CE = paste0(unique(CE), collapse = "/"),
              preMZ = mean(preMZ),
              replicate = unique(replicate)
              ), by = list(sample, cluster)
            ][, cluster := NULL]

            t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
            t[is.na(precursor), precursor := FALSE]
            t[, id := idf$id][]
            setorder(t, sample, mz)
            setcolorder(t, c("sample", "replicate", "id", "mz", "intensity", "CE", "preMZ", "precursor"))

          } else { #when NULL do not merge by samples
            t <- t[, .(
              mz = mean(mz),
              intensity = sum(intensity) / length(unique(seqNum)),
              CE = paste0(unique(CE), collapse = "/"),
              preMZ = mean(preMZ)
              ), by = list(cluster)
            ][, cluster := NULL]

            t[mz >= idf$mzmin[1] & mz <= idf$mzmax[1], precursor := TRUE]
            t[is.na(precursor), precursor := FALSE]
            t[, id := idf$id][]
            setorder(t, mz)
            setcolorder(t, c("id", "mz", "intensity", "CE", "preMZ", "precursor"))

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

          } else if (!is.null(mergeBy)) { #mergeBy samples when mergeBy is not null
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

  }

  msnList <- msnList[intensity >= minIntensityPost, ]

  return(msnList)

}






































     # mzMat2 <- dist(msnSpec2$mz, method = "euclidean")
      # mzMat2 <- as.data.table(as.matrix(mzMat2))
      # mzMat2 <- mzMat2[, lapply(.SD, function(x) x / msnSpec2$mz * 1E6), .SDcols = colnames(mzMat2)]

      

        #another method
        #mzdiff <- abs(diff(spcomb$mz))
        #spcomb[, cluster := 1 + c(0, cumsum(mzdiff > clusterMzWindow))]

    #   temp <- patRoon:::getSpectraHeader(spectra, rtRange, 2, scans$precursorMZ[x], 4)

    #   avgFeatParams <- patRoon::getDefAvgPListParams()
    #   avgFeatParams$minIntensityPost <- 10
    #   avgFeatParams$minIntensityPre <- 10

    #   avgFeatParamsMS <- avgFeatParamsMSMS <-
    #     avgFeatParams[setdiff(names(avgFeatParams), c("pruneMissingPrecursorMS", "retainPrecursorMSMS"))]
    # avgFeatParamsMS$retainPrecursor <- TRUE;
    # avgFeatParamsMS$pruneMissingPrecursor <- avgFeatParams$pruneMissingPrecursorMS
    # avgFeatParamsMSMS$pruneMissingPrecursor <- FALSE
    # avgFeatParamsMSMS$retainPrecursor <- avgFeatParams$retainPrecursorMSMS




    #   temp2 <- patRoon:::averageSpectraMZR(
    #       spectra = spectra,
    #       hd = temp,
    #       precursor = scans$precursorMZ[x],
    #       clusterMzWindow = 0.05,
    #       topMost = 50,
    #       minIntensityPre = 10,
    #       minIntensityPost = 10,
    #       avgFun = avgFeatParamsMS$avgFun,
    #       method = avgFeatParamsMS$method,
    #       pruneMissingPrecursor = FALSE,
    #       retainPrecursor = TRUE)












makeFileHashCopy <- function(...) digest::digest(sapply(list(...), digest::digest, file = TRUE, algo = "xxhash64"))

recursiveApplyDTCopy <- function(l, f, appl = lapply, ...) {
  rec <- function(x) {
    if (isS4(x)) {
      for (sn in slotNames(x)) {
        slot(x, sn) <- rec(slot(x, sn))
      }
    }
    else if (is.list(x)) {
      if (is.data.table(x)) {
        x <- f(x)
      } else {
        # retain attributes: https://stackoverflow.com/a/48905113
        x <- "attributes<-"(appl(x, rec, ...), attributes(x))
      }
    }
    return(x)
  }
  return(rec(l))
}

prepareDTForComparisonCopy <- function(dt) {
    setattr(dt, ".internal.selfref", NULL)
    setindex(dt, NULL)
}

makeHashCopy <- function(..., checkDT = TRUE) {
  args <- list(...)
  if (checkDT) {
    # strip DT self refs as they sometimes mess up hashing
    args <- recursiveApplyDTCopy(args, function(dt) prepareDTForComparisonCopy(copy(dt)), sapply, simplify = FALSE)
  }
  return(digest::digest(args, algo = "xxhash64"))
}






# makeFileHashCopy(filepath)
# filepath <- example01@samples$file[1]
# rtRange <- c(800, 900)
# args <- list(makeFileHashCopy(filepath), rtRange)


loadSpectraCopy <- function(filepath = NULL, rtRange = NULL, verbose = TRUE, cacheDB = NULL)
{
    hash <- makeHashCopy(makeFileHashCopy(filepath), rtRange)
    ret <- loadCacheData("specData", hash, cacheDB)
    if (!is.null(ret) && length(ret$spectra) > 1 && is.data.table(ret$spectra[[1]]))
        ret <- NULL # old (pre v1.1) format, ignore cache to avoid crashes with Rcpp interface
    if (is.null(ret))
    {
        if (verbose)
            printf("Loading raw spectra for '%s'...\n", path)
        msf <- mzR::openMSfile(path)
        hd <- as.data.table(mzR::header(msf))

        if (is.null(rtRange))
            ps <- mzR::peaks(msf) # load all
        else
            ps <- mzR::peaks(msf, hd[numGTE(retentionTime, rtRange[1]) & numLTE(retentionTime, rtRange[2]), seqNum])

        ret <- list(header = hd, spectra = ps)
        mzR::close(msf)
        saveCacheData("specData", ret, hash, cacheDB)
    }

    return(ret)
}










#' @title importRawData
#' @description Function to import MS data from the MS files listed in
#' the \linkS4class{ntsData} object. Files should contain centroided spectra.
#' The function \code{\link[MSnbase]{readMSData}} from the
#' \code{MSnbase} package is used to read the MS files.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param rtFilter A numeric vector with length 2 defining the minimum
#' and maximum chromatographic retention time for the listed MS files.
#' @param rtUnit The unit of the \code{rtFilter}.
#' Possible values are \code{min} (the default) and \code{sec}.
#' @param msLevel The MS dimensions for the rtFilter to be applied.
#' The default is both MS1 and MS2 using \code{c(1,2)}.
#' @param centroidedData Logical, set to \code{TRUE} for MS files
#' with centroided data or \code{FALSE} for profile data.
#' \code{NA} will collect all the data from the MS files.
#' @param removeEmptySpectra Logical, set to TRUE if empty spectra should be removed.
#' It is recommended to remove empty spectra as it may cause issues during creation of features.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return The \linkS4class{ntsData} object including a standard
#' \linkS4class{OnDiskMSnExp} object in the MSnExp slot. Note, that
#' the \linkS4class{OnDiskMSnExp} object can also be used within
#' the workflow of \pkg{Bioconductor} packages.
#'
#' @references
#' \insertRef{MSnbase1}{ntsIUTA}
#' \insertRef{MSnbase2}{ntsIUTA}
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom BiocParallel SnowParam register
#' @importFrom parallel detectCores
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importFrom MSnbase readMSData
#' @importMethodsFrom MSnbase filterRt filterEmptySpectra
#' @importFrom methods new
#'
#' @examples
#' path <- system.file(package = "ntsIUTA", dir = "extdata")
#' dt <- setupProject(path = path, save = FALSE)
#' dt <- importRawData(dt[1], save = FALSE, centroidedData = TRUE)
#'
importRawData <- function(
  obj = NULL,
  rtFilter = NULL,
  rtUnit = "min",
  msLevel = c(1, 2),
  centroidedData = TRUE,
  removeEmptySpectra = TRUE,
  save = FALSE
) {
  assertClass(obj, "ntsData")

  snow <- SnowParam(
    workers = detectCores() - 1,
    type = "SOCK",
    exportglobals = FALSE,
    progressbar = TRUE
  )

  register(snow, default = TRUE)

  msFiles <- obj@samples$file[drop = TRUE]
  sample_name <- obj@samples$sample
  sample_group <- obj@samples$group

  if (length(sample_name) == 0) {
    warning("There are not samples in the ntsData object.")
    return(obj)
  }

  raw <- suppressWarnings(
    readMSData(
      msFiles,
      pdata = new("NAnnotatedDataFrame",
      data.frame(
        sample_name = sample_name,
        sample_group = sample_group)
      ),
      msLevel. = NULL,
      mode = "onDisk",
      centroided. = centroidedData,
      smoothed. = FALSE
    )
  )

  if (!is.null(rtFilter)) {
    if (rtUnit == "min") rtFilter <- rtFilter * 60
    raw <- filterRt(raw, rt = rtFilter, msLevel. = msLevel)
  }

  if (removeEmptySpectra) raw <- filterEmptySpectra(raw)

  obj@MSnExp <- raw

  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)
}




#' @title centroidProfileData
#' @description Centroiding of profile data with additional possibility
#' for data smoothing before centroiding and \emph{m/z} refinement.
#' The \code{centroidProfileData} function combines functions \code{smooth}
#' and \code{pickPeaks} from the \code{MSnbase} package, see references.
#'
#' @param obj A \linkS4class{ntsData} object with profile data for centroiding.
#' @param halfwindow Sets the window size for centroiding as \code{2 * halfwindow + 1}.
#' The \code{halfwindow} should be slightly larger than the full width
#' at half maximum of the profile peak.
#' @param SNR The signal-to-noise ratio to consider a local maximum as peak.
#' @param noiseMethod The method to estimate the noise level.
#' Possible methods are "MAD" (the default) and "SuperSmoother".
#' See \code{?MSnbase::pickPeaks} for more information.
#' @param smoothing Logical, set to \code{TRUE} for applying smothing
#' to the profile data before centroiding. The default is FALSE.
#' @param methodSmoothing Method for data smoothing.
#' The possible methods are "SavitzkyGolay" (the default) and "MovingAverage".
#' See \code{?MSnbase::smooth} for more information and arguments,
#' which are passed by \code{...}.
#' @param ... Arguments for selected smoothing method.
#' See \code{?MSnbase::smooth} for possible arguments for each method.
#' @param methodRefineMz Method for refinement.
#' Possible methods are "none" (the default, for not applying \emph{m/z} refinement),
#' "kNeighbors" and "descendPeak". See \code{?MSnbase::pickPeaks} for more information.
#' @param k When refine method is "kNeighbors",
#' \code{k} is number of closest signals to the centroid.
#' @param signalPercentage When refine method is "descendPeak",
#' \code{signalPercentage} is the minimum signal percentage of centroids to refine \emph{m/z}.
#' @param stopAtTwo Logical, when refine method is "descendPeak",
#' set to \code{TRUE} for allowing two consecutive equal or higher signals.
#' \code{FALSE} will stop when one equal or higher centroid is found.
#' @param save Logical, set to \code{TRUE} to replace
#' the original files by the centroided files in disk.
#' The location is taken from the originbal file paths.
#'
#' @return Centroided \linkS4class{ntsData} object.
#' When \code{save} is set to TRUE, the profile data in the original
#' mzML or mzXML files is replaced by the centroided data.
#'
#' @references
#' \insertRef{MSnbase2}{ntsIUTA}
#'
#' @export
#'
#' @importClassesFrom MSnbase OnDiskMSnExp
#' @importMethodsFrom MSnbase fileNames smooth pickPeaks writeMSData
#'
centroidProfileData <- function(obj,
                                halfwindow = 2,
                                SNR = 0,
                                noiseMethod = "MAD",
                                smoothing = FALSE,
                                methodSmoothing = "SavitzkyGolay",
                                methodRefineMz = "kNeighbors",
                                k = 1,
                                signalPercentage = 10, stopAtTwo = TRUE,
                                save = FALSE, ...) {

  raw <- obj@MSnExp

  if (smoothing) {
    raw <- raw %>% MSnbase::smooth(method = methodSmoothing, ...)
  }

  if (methodRefineMz == "kNeighbors") {
    raw <- pickPeaks(raw,
                     halfWindowSize = halfwindow,
                     SNR = SNR,
                     noiseMethod = noiseMethod,
                     refineMz = methodRefineMz,
                     k = k)
  } else {
    if (methodRefineMz == "descendPeak") {
      raw <- pickPeaks(raw,
                       halfWindowSize = halfwindow,
                       SNR = SNR,
                       noiseMethod = noiseMethod,
                       refineMz = methodRefineMz,
                       signalPercentage = signalPercentage,
                       stopAtTwo = TRUE)
    } else {
      raw <- pickPeaks(raw,
                       halfWindowSize = halfwindow,
                       SNR = SNR,
                       noiseMethod = noiseMethod,
                       refineMz = "none")
    }
  }

  obj@MSnExp <- raw

  if (save) {
    fls_new <- fileNames(raw)
    writeMSData(raw, file = fls_new)
  }

  return(obj)

}
