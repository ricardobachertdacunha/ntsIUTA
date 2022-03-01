

#' @title importRawData_old
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
importRawData_old <- function(
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




#' @title centroidProfileData_old
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
centroidProfileData_old <- function(obj,
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



### components_Old -----

#' @describeIn ntsData Getter for components (i.e., annotated features).
#'
#' @param object An \linkS4class{ntsData} object.
#' @param samples The indice/s or name/s of samples to keep in the \code{object}.
#' @param ID The ID of features of interest.
#' @param mz Alternatively to \code{ID}, the \emph{m/z} of interest.
#' can be of length two, defining a mass range.
#' @param ppm The mass deviation, in ppm, of a given \code{mz}.
#' @param rt The retention time to find features.
#' @param rtWindow The time deviation. Can be of length two, defining a time range.
#' @param rtUnit The unit of the time arguments. Possible values are "sec" and "min".
#' @param compNumber Alternatively, the component number to find features.
#' @param entireComponents Logical, set to \code{TRUE} to extract all features
#' from the represented components.
#' @param onlyAnnotated Logical, set to \code{TRUE} to extract only annotated features.
#' @param onlyRelated Logical, set to \code{TRUE} to extract only features that are related
#' to the features of interest.
#'
#' @export
#'
#' @importFrom checkmate assertSubset
#' @importFrom dplyr filter between
#' @importFrom stats na.omit
#'
setMethod("components", "ntsData", function(object,
                                            samples = NULL,
                                            ID = NULL,
                                            mz = NULL, ppm = 5,
                                            rt = NULL, rtWindow = 1, rtUnit = "sec",
                                            compNumber = NULL,
                                            entireComponents = TRUE,
                                            onlyAnnotated = FALSE,
                                            onlyRelated = TRUE) {

  if (missing(samples)) samples <- NULL

  if (missing(ID)) ID <- NULL

  if (missing(mz)) mz <- NULL

  if (missing(compNumber)) compNumber <- NULL

  comp <- object@annotation$comp

  if (nrow(comp) == 0) return(comp)

  rg <- unique(object@samples$group)

  if (!is.null(samples)) {
    if (is.character(samples)) {
      rg <- unique(object@samples$group[object@samples$sample %in% samples])
    } else {
      rg <- unique(object@samples$group[samples])
    }
  }

  comp <- comp[comp$group %in% rg, ]

  if (nrow(comp) == 0) return(comp)

  if (!is.null(ID)) {
    comp <- comp[comp$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      if (missing(rt)) rt <- NULL
      if (missing(rtWindow)) rtWindow <- NULL
      if (missing(rtUnit)) rtUnit <- "sec"
      if (missing(ppm)) ppm <- 20
      assertSubset(rtUnit, c("sec", "min"))
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      comp <- dplyr::filter(comp,
                            between(mz, mzr[1], mzr[2]),
                            between(rt, rtr[1], rtr[2]))
    } else {
      if (!is.null(compNumber)) {
        comp <- comp[comp$comp %in% compNumber, ]
      }
    }
  }

  if (nrow(comp) == 0) return(comp)

  comp2 <- comp

  if (!missing(entireComponents)) {
    if (entireComponents) {
      comp2 <- object@annotation$comp[object@annotation$comp$comp %in% comp2$comp |
                                        object@annotation$comp$isogroup %in% comp2$isogroup, ]
      comp2 <- comp2[comp2$group %in% rg, ]
    }
  }

  if (!missing(onlyAnnotated)) {
    if (onlyAnnotated) comp2 <- comp2[!is.na(comp2$isoclass) | (comp2$nAdducts > 0), ]
  }

  if (!missing(onlyRelated)) {
    if (onlyRelated) {
      isos <- comp2$ID[comp2$Mion %in% stats::na.omit(comp$Mion)]
      comp2 <- comp2[comp2$ID %in% unique(isos), ]
    }
  }

  return(comp2)

})





#' @title plotPeaksFunc_Old
#' @description Method for plotting peaks from a \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param samples The index or name of the sample/s.
#' The default is \code{NULL} and all samples are used.
#' @param ID The identifier of the peaks of interest.
#' When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.
#' @param mz Optional target \emph{m/z} to find peaks using
#' the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the peaks
#' when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect peaks.
#' The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}.
#' The default is \code{min}.
#' @param colorBy Possible values are \code{"peaks"} (the default),
#' \code{"samples"} or \code{sampleGroups},
#' for colouring by peaks, samples or sample replicate groups, respectively.
#'
#' @return A peak/s map plot produced through \pkg{plotly}.
#'
#' @export
#'
#' @importFrom checkmate assertClass assertSubset
#' @importFrom dplyr between filter
#' @importFrom plotly toRGB plot_ly add_trace layout add_segments
#'
plotPeaksFunc_Old <- function(obj, samples = NULL,
                          ID = NULL,
                          mz = NULL, ppm = 20,
                          rt = NULL, rtWindow = NULL,
                          rtUnit = "sec",
                          colorBy = "samples") {
  
  assertClass(obj, "ntsData")
  
  assertSubset(rtUnit, c("sec", "min"))
  
  assertSubset(colorBy, c("peaks", "samples", "sampleGroups"))
  
  if (!is.null(samples)) obj <- filterFileFaster(obj, samples)
  
  rtr <- NULL
  
  if (!is.null(ID)) {
    pki <- obj@peaks[obj@peaks$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      if (is.null(rtr)) rtr <- c(min(obj@peaks$rtmin), max(obj@peaks$rtmax))
      pki <- dplyr::filter(obj@peaks,
                           dplyr::between(mz, mzr[1], mzr[2]),
                           dplyr::between(rt, rtr[1], rtr[2]))
    } else {
      return(cat("One of ID or mz should be given."))
    }
  }
  
  if (nrow(pki) == 0) return(cat("No features found."))
  
  if (is.null(rtr)) rtr <- c(min(pki$rtmin) - 60, max(pki$rtmax) + 60)
  
  if (is.null(ppm)) ppm <- 5
  
  sp <- obj@samples$sample
  
  rg <- obj@samples$group
  
  if (colorBy == "samples") {
    colors <- getColors(obj, "samples")
    col_val <- apply(pki, MARGIN = 1, FUN = function(x) colors[names(colors) %in% x["sample"]])
    leg <- pki$sample
    sleg <- !duplicated(leg)
  } else {
    if (colorBy == "peaks") {
      colors <- getColors(nrow(pki))
      col_val <- colors
      leg <- paste0( pki$ID, "/", round(pki$mz, digits = 4), "/", round(pki$rt, digits = 0))
      sleg <- !duplicated(leg)
    } else {
      colors <- getColors(obj, "sampleGroups")
      col_val <- apply(pki, MARGIN = 1, FUN = function(x) colors[names(colors) %in% x["sample"]])
      leg <- pki$group
      sleg <- !duplicated(leg)
    }
  }
    
  EIC <- extractEIC(obj,
                    mz = c(min(pki$mzmin), max(pki$mzmax)),
                    rtWindow = rtr, rtUnit = "sec")
  
  plot <- plot_ly()
  
  for (i in seq_len(nrow(pki))) {
    
    df <- EIC[EIC$mz >= pki$mzmin[i] &
                EIC$mz <= pki$mzmax[i] &
                EIC$file == which(sp == pki$sample[i]), ]
    
    
    
    plot <- plot %>% add_trace(df,
                               x = df$rt,
                               y = df$i,
                               type = "scatter", mode = "markers",
                               marker = list(size = 0.2,
                                             color = col_val[i]),
                               connectgaps = TRUE,
                               name = leg[i],
                               legendgroup = leg[i],
                               showlegend = sleg[i]
    )
    
    df <- df[df$rt >= pki$rtmin[i] & df$rt <= pki$rtmax[i], ]
    
    plot <- plot %>%  add_trace(df,
                                x = df$rt,
                                y = df$i,
                                type = "scatter", mode =  "markers",
                                fill = "tozeroy", connectgaps = TRUE,
                                fillcolor = paste(color = col_val[i], 50, sep = ""),
                                #line = list(width = 0.1, color = col_val[i]),
                                marker = list(size = 3, color = col_val[i]),
                                name = leg[i],
                                legendgroup = leg[i],
                                showlegend = FALSE,
                                hoverinfo = "text",
                                text = paste("</br> peak: ", pki$ID[i],
                                             "</br> sample: ", pki$sample[i],
                                             "</br> <i>m/z</i>: ", round(pki$mz[i], digits = 4),
                                             "</br> dppm: ", round(((pki$mzmax[i] - pki$mzmin[i]) / pki$mz[i]) * 1E6, digits = 0),
                                             "</br> rt: ", round(pki$rt[i], digits = 0),
                                             "</br> drt: ", round(pki$rtmax[i] - pki$rtmin[i], digits = 0),
                                             "</br> Int: ", round(pki$intensity[i], digits = 0),
                                             "</br> Filled: ",
                                             if ("is_filled" %in% colnames(pki)) {
                                               ifelse(pki$is_filled[i] == 1, TRUE, FALSE)
                                             } else {
                                               FALSE
                                             }))
    
    plot <- plot %>% add_segments(x = pki$rt[i], xend = pki$rt[i], y = 0, yend = pki$intensity[i],
                                  legendgroup = leg[i], showlegend = FALSE, line = list(color = col_val[i], size = 0.5))
    
  }
  
  dppm <- c(round(min(pki$mzmin), digits = 4), round(max(pki$mzmax), digits = 4))
  dppm_val <- round((dppm[2] - dppm[1])/dppm[2] * 1E6, digits = 0)
  drt <- c(round(min(pki$rtmin), digits = 0), round(max(pki$rtmax), digits = 0))
  
  title_text <- paste0("<i>m/z</i>",": ", dppm[1], " - ", dppm[2],
                       " (", dppm_val, " ppm)",
                       "  rt: ", drt[1], " - ", drt[2],
                       " (",  drt[2] - drt[1], " sec.)")
  
  title <- list(text = title_text, x = 0.1, y = 0.98,
                font = list(size = 9, color = "black"))
  
  xaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Retention Time (sec.)",
                titlefont = list(size = 12, color = "black"),
                range = c(rtr[1], rtr[2]),
                autotick = TRUE, ticks = "outside")
  
  yaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Intensity",
                titlefont = list(size = 12, color = "black"))
  
  plot <- plot %>% plotly::layout(xaxis = xaxis,
                                  yaxis = yaxis,
                                  title = title)
  
  return(plot)
  
  # TODO Idea for 3D plotly
  # if (td) {
  #   
  #   df <- EIC
  #   df$file <- sapply(df$file,FUN = function(x) sp[x])
  #   df$file <- factor(df$file, levels = sp, ordered = FALSE)
  #   
  #   plot_ly(df, x = ~rt, y = ~mz, z = ~i, color = ~file, colors = c('#BF382A', '#0C4B8E'),
  #           marker = list(size = 2), mode = "scatter3d")
  # }
  
}






#' @title removeDuplicateNamesNOTUSED
#' @description Removes duplicate sample names to avoid problems during further analysis.
#'
#' @param newJoint Original and new samples \code{data.frame} objects collated by \code{rbind}.
#'
#' @return Returns the same \code{data.frame} but with duplicate sample names
#' edited as duplicate sample names are not allowed.
#'
#' @importFrom dplyr semi_join
#' @importFrom stringr str_extract str_pad
#'
removeDuplicateNamesNOTUSED <- function(newJoint) {
  duplicates <- duplicated(newJoint$sample)
  if (TRUE %in% duplicates) {
    print(paste("Number of duplicate names found:", length(duplicates[duplicates == TRUE]), sep = " "))
    for (i in seq_len(length(duplicates))) {
      if (duplicates[i]) {
        endOfString <-  stringr::str_extract(newJoint$sample[i], "_[^_]+$")
        startOfString <-  stringr::str_extract(newJoint$sample[i], paste("^.*(?=(", endOfString, "))", sep = ""))
        if ((grepl("_[0-9]+", endOfString)) && (!is.na(endOfString))) {
          endOfString <- as.numeric(stringr::str_extract(endOfString, "[0-9]+")) + 1
          endOfString <- stringr::str_pad(endOfString, 2, pad = "0")
          newJoint$sample[i] <- paste(startOfString, endOfString, sep = "_")
        } else {
          newJoint$sample[i] <- paste(newJoint$sample[i], "_01", sep = "")
        }
      }
    }
    duplicates <- duplicated(newJoint$sample)
    if (TRUE %in% duplicates) {
      newJoint <- removeDuplicateNames(newJoint)
    }
  }
  return(newJoint)
}


#' @title annotateFeatures_Old
#'
#' @description Group features into components according to co-elution and
#' EIC similarities and annotate isotopes and adducts in each component
#' per sample replicate group.
#'
#' @param obj An \linkS4class{ntsData} object with features.
#' @param algorithm The algorithm for finding isotopes. Possible values are
#' "alteredcamera", "camera" or "ramclustr".
#' @param param The param used for annotation of isotopes and adducts.
#' @param excludeBlanks Logical, set to \code{TRUE} for excluding
#' blank replicate groups from annotation.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object, containing components with
#' annotated features.
#'
#' @export
#'
#' @importFrom checkmate testChoice assertClass
#' @importFrom patRoon getXCMSSet
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr findAdducts getPeaklist
#' @importFrom utils txtProgressBar setTxtProgressBar read.table
#' @importFrom dplyr rename select everything
#'
annotateFeatures_Old <- function(obj = NULL,
                             algorithm = NULL,
                             param = NULL,
                             excludeBlanks = FALSE,
                             save = FALSE
) {
  
  assertClass(obj, "ntsData")
  
  if (is.null(algorithm)) if (!is.na(obj@algorithms$annotation)) algorithm <- obj@algorithms$annotation
  
  if (!testChoice(algorithm, c("alteredcamera", "camera", "ramclustr", "nontarget", "intclust"))) {
    warning("Algorithm not recognized. See ?annotateFeatures for more information.")
    return(obj)
  }
  
  if (is.null(param)) if (length(obj@parameters$annotation) > 0) param <- obj@parameters$annotation
  
  if (length(param) == 0) {
    warning("Parameters for annotation must be given!")
    return(obj)
  }
  
  if (excludeBlanks) {
    rg <- unique(obj@samples$group[!(obj@samples$group %in% obj@samples$blank)])
  } else {
    rg <- unique(obj@samples$group)
  }
  
  if (!(algorithm == "alteredcamera")) {
    
    ag <- list(fGroups = obj@patdata, algorithm = algorithm)
    
    pat <- do.call(generateComponents, c(ag, param))
    
    # TODO convert to data frame from other functions and integrate with df scheme
    
  } else {
    
    xs <- obj@patdata
    
    xs <- getXCMSSet(xs, verbose = TRUE, loadRawData = TRUE)
    
    xA <- xsAnnotate(xs = xs,
                     sample = seq_len(nrow(obj@samples)),
                     polarity = obj@polarity)
    
    xA <- groupFWHM(xA, sigma = param@sigma,
                    perfwhm = param@perfwhm,
                    intval = "maxo")
    
    xA <- groupCorr(xA,
                    cor_eic_th = param@cor_eic_th,
                    cor_exp_th = param@cor_exp_th,
                    pval = param@pval,
                    graphMethod = "hcs",
                    calcIso = FALSE,
                    calcCiS = TRUE,
                    calcCaS = param@calcCaS,
                    psg_list = NULL,
                    xraw = NULL,
                    intval = "maxo")
    
    xAL <- list()
    
    pb <- txtProgressBar(min = 0, max = 100, initial = 0, char = "=", width = 80, style = 3)
    
    for (rgidx in seq_len(length(rg))) {
      
      setTxtProgressBar(pb, ((rgidx / length(rg)) * 100))
      
      sampleidxs <- which(obj@samples$group == rg[rgidx])
      
      xA_temp <- FindIsotopesWithValidationAltered(
        xA = xA,
        obj = obj,
        sampleidxs = sampleidxs,
        ppm = param@ppmIsotopes,
        mzabs = param@mzabs,
        noise = param@noise,
        maxcharge = 3,
        intval = "maxo",
        validateIsotopePatterns = param@validateIsotopePatterns
      )
      
      if (param@searchAdducts) {
        
        if (param@extendedList) {
          rules_pos <- system.file("rules/extended_adducts_pos.csv", package = "CAMERA")
          rules_neg <- system.file("rules/extended_adducts_neg.csv", package = "CAMERA")
        } else {
          rules_pos <- system.file("rules/primary_adducts_pos.csv", package = "CAMERA")
          rules_neg <- system.file("rules/primary_adducts_neg.csv", package = "CAMERA")
        }
        
        if (obj@polarity == "positive") {
          rules <- utils::read.table(rules_pos, header = TRUE, sep = ",")
          if (length(colnames(rules)) == 1) {
            rules <- utils::read.table(rules_pos, header = TRUE, sep = "")
          }
          rules <- rbind(rules, data.frame(name = "[M+]",
                                           nmol = 1,
                                           charge = 1,
                                           massdiff = 0,
                                           oidscore = 12,
                                           quasi = 1,
                                           ips = 1))
        }
        
        if (obj@polarity == "negative") {
          rules <- utils::read.table(rules_neg, header = TRUE, sep = ",")
          if (length(colnames(rules)) == 1) {
            rules <- utils::read.table(rules_neg, header = TRUE, sep = "")
          }
          rules <- rbind(rules, data.frame(name = "[M-]",
                                           nmol = 1,
                                           charge = 1,
                                           massdiff = 0,
                                           oidscore = 12,
                                           quasi = 1,
                                           ips = 1))
        }
        
        xA_temp <- findAdducts(xA_temp,
                               ppm = param@ppmAdducts,
                               mzabs = 0,
                               multiplier = 2,
                               polarity = obj@polarity,
                               rules = rules,
                               max_peaks = 100,
                               psg_list = NULL)
      }
      
      xAL[[rg[rgidx]]] <- xA_temp
      
    }
    
    #prepare data frame
    for (r in seq_len(length(xAL))) {
      df_temp <- CAMERA::getPeaklist(xAL[[r]], intval = "maxo")
      df_temp$group <- names(xAL)[r]
      df_temp$ID <- obj@features$ID
      if (r == 1 | length(xAL) == 1) {
        df <- df_temp
      } else {
        df <- rbind(df, df_temp)
      }
    }
    
  }
  
  isotopes <- df$isotopes
  isotopes <- strsplit(isotopes, split = "\\]\\[")
  
  isogroup <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][1])
  isogroup <- str_extract(isogroup, pattern = "([0-9]+)")
  noIsogroup <- is.na(isogroup)
  isogroup <- paste0(isogroup, "_", df$group)
  isogroup[noIsogroup] <- NA
  
  isoclass <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][2])
  isoclass <- as.vector(ifelse(!is.na(isoclass), paste0("[",isoclass), NA))
  
  df$isogroup <- isogroup
  df$isoclass <- isoclass
  
  Mion <- df$mz
  index <- seq_len(length(Mion))
  isMion <- grepl(pattern = "[M]", isoclass, fixed = TRUE)
  
  charge <- sapply(index, FUN =  function(x) {
    ifelse(isMion[x], as.numeric(str_extract(isoclass[x], pattern = "([0-9]+)")), NA)
  })
  charge[is.na(charge)] <- 1
  
  if (obj@polarity == "positive") {
    Mion <- (Mion - 1.007276) * charge
  } else {
    Mion <- (Mion + 1.007276) * charge
  }
  
  Mion <- round(Mion, digits = 3)
  
  adducts <- strsplit(df$adduct, split = " ")
  df$conflits <- unlist(lapply(X = seq_len(length(adducts)), function(x) ifelse(length(adducts[[x]]) > 2, TRUE, FALSE)))
  df$adductclass <- I(lapply(X = seq_len(length(adducts)), function(x) grep("M", adducts[[x]], value = TRUE)))
  df$adductMion <- I(lapply(X = seq_len(length(adducts)), function(x) as.numeric(grep("M", adducts[[x]], value = TRUE, invert = TRUE))))
  possibleAdducts <- unlist(lapply(df$adductMion, length))
  df$nAdducts <- possibleAdducts
  isAdduct <- possibleAdducts != 0
  
  Mion <- unlist(lapply(index, FUN = function(x) ifelse(isAdduct[x] & !df$conflits[x], as.numeric(unlist(df$adductMion[x])), Mion[x])))
  
  isopolog <- !isMion
  isopolog[is.na(isoclass)] <- FALSE
  
  names(Mion) <- isogroup
  
  Mion <- unlist(lapply(index, function(x) ifelse(isopolog[x], Mion[names(Mion) %in% isogroup[x]], Mion[x])))
  
  df$Mion <- as.numeric(Mion)
  
  df <- rename(df, comp = pcgroup)
  df <- select(df, ID, group, comp, Mion, isoclass, isogroup, adductclass, adductMion, conflits, nAdducts, everything())
  
  obj@annotation$comp <- df
  
  obj@annotation$raw <- xAL
  
  obj@algorithms$annotation <- algorithm
  
  obj@parameters$annotation <- param
  
  if (save) saveObject(obj = obj)
  
  if (is.character(save)) saveObject(obj = obj, filename = save)
  
  return(obj)
  
}


#' @title getScreeningListTemplate_Old
#'
#' @param projPath The project folder location. Default is \code{setup$projPath}.
#'
#' @return Pastes a template .csv file of the screeningList into the projPath.
#'
#' @export
#'
getScreeningListTemplate_Old <- function(projPath = setup$projPath) {
  base::file.copy(from = base::paste0(base::system.file(package = "ntsIUTA", dir = "extdata"), "/screeningList_template.csv"),
                  to = setup$projPath,
                  overwrite = FALSE)
}


### subsetting ntsData_Old ----

#' @describeIn ntsData Subset on samples, using sample index or name.
#'
#' @param x An \linkS4class{ntsData} object.
#' @param i The indice/s or name/s of the samples to keep in the \code{x} object.
#' @param j Ignored.
#' @param drop Ignored.
#' @param \dots Ignored.
#'
#' @export
#'
#' @importMethodsFrom MSnbase filterFile fileNames
#' @importMethodsFrom xcms filterFile
#'
setMethod("[", c("ntsData", "ANY", "missing", "missing"), function(x, i, ...) {
  
  if (!missing(i)) {
    
    if (!is.character(i)) {
      sn <- x@samples$sample[i]
      sidx <- which(x@samples$sample %in% sn)
    } else {
      sn <- i
      sidx <- which(x@samples$sample %in% sn)
    }
    
    x@samples <- x@samples[x@samples$sample %in% sn,, drop = FALSE]
    
    x@metadata <- x@metadata[x@metadata$sample %in% sn,, drop = FALSE]
    
    if (length(x@MSnExp) > 0) x@MSnExp <- MSnbase::filterFile(x@MSnExp, file = sidx)
    
    if (nrow(x@peaks) > 0) x@peaks <- x@peaks[x@peaks$sample %in% sn,, drop = FALSE]
    
    if (nrow(x@features) > 0) {
      
      #To be replaced when sub-setting in patRoon is fixed
      x@patdata <- filterFeatureGroups(x@patdata, i = sidx)
      
      # TODO subsetting rebuilds feature list, interferse with filtering of feature list
      x <- buildFeatureList(x)
      
    } else {
      
      #subset patRoon object without featureGroups
      if (length(x@patdata) > 0) {
        x@patdata@features <- x@patdata@features[sidx]
        x@patdata@analysisInfo <- x@patdata@analysisInfo[x@patdata@analysisInfo$analysis %in% sn, ]
        if (x@patdata@algorithm == "xcms3") {
          files <- basename(fileNames(x@patdata@xdata)[sidx])
          x@patdata@xdata <- filterFile(x@patdata@xdata, files)
        }
      }
    }
    
    #annotation, remove replicate groups without samples
    if (nrow(x@annotation$comp) > 0) {
      rg <- unique(x@samples$group)
      x@annotation$comp <- x@annotation$comp[x@annotation$comp$group %in% rg, ]
      x@annotation$raw <- x@annotation$raw[names(x@annotation$raw) %in% rg]
      if (nrow(x@features) > 0) x@annotation$comp <- x@annotation$comp[x@annotation$comp$ID %in% x@features$ID, ]
    }
  }
  
  return(x)
  
})


#' @title buildFeatureList_Old
#' @description Function to create a \code{data.frame} with detailed information for each feature in the given a \linkS4class{XCMSnExp} object object.
#' Optionally, a list of \linkS4class{xsAnnotate} objects per replicate group as obtained by \code{makeFeatureComponents} can be given.
#' The information from annotated isotopes and adducts for each component will be added to the \code{data.frame}.
#' 
#' @param x A \linkS4class{XCMSnExp} object.
#' @param xA A list of \linkS4class{xsAnnotate} objects per replicate group.
#' @param xPat A \linkS4class{featureGroups} object converted by \code{\link{getPatData}}.
#' @param snWindow Time in seconds to expand the peak width for calculation of the signal-to-noise ratio.
#' @param save Logical, set to \code{TRUE} to save the generated \code{fl} object in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#'
#' @return A \code{data.frame} with detailed information for each feature in the given objects.
#' 
#' @export
#' 
#' @import magrittr
#' @importFrom BiocParallel registered SnowParam register bpparam bplapply
#' @importFrom parallel detectCores
#' @importFrom xcms filterMsLevel chromPeaks featureDefinitions
#' @importFrom methods as
#' @importMethodsFrom MSnbase fileNames
#' @importFrom patRoon importFeatureGroupsXCMS3 as.data.table
#' @importFrom dplyr select mutate group_by count filter between all_of everything
#' @importFrom stats quantile sd na.omit
#' @importFrom CAMERA getPeaklist
#' @importFrom stringr str_extract
#' 
#' 
#'
#' @examples
#' 
#' 
#' 
buildFeatureList_Old <-  function(x = featData, xA = featComp, xPat = patData,
                                  snWindow = 240,
                                  save = TRUE, projPath = setup$projPath){
  
  maxMultiProcess = TRUE
  if (maxMultiProcess) {
    snow <- registered("SnowParam")
    if (snow$workers < detectCores()) {
      snow <- SnowParam(workers = detectCores() - 1,
                        type = "SOCK",
                        exportglobals = FALSE,
                        progressbar = TRUE)
      register(snow, default = TRUE)
    }
  }
  
  # x = featData
  # xA = featComp
  # xPat = patData
  # snWindow = 240
  
  #collect centroids
  base::cat("Extracting centroids...")
  base::cat("\n")
  cent <- xcms::filterMsLevel(x, msLevel. = 1)
  cent <- base::as.data.frame(methods::as(cent, "data.frame"))
  cent$rt <- base::as.numeric(cent$rt)
  cent$mz <- base::as.numeric(cent$mz)
  cent$i <- base::as.numeric(cent$i)
  base::colnames(cent) <- c("file", "rt", "mz", "into")
  
  #chromPeaks
  # TODO Peaks IDs when filled adds a zero which is not there for raw peaks IDs
  chromPeaks <- xcms::chromPeaks(x, isFilledColumn = TRUE, msLevel = 1)
  
  
  #collect features
  fl <- xcms::featureDefinitions(x)
  fl <- base::as.data.frame(fl)
  fl$FT <- base::row.names(fl)
  fl <- fl[,!(base::colnames(fl) %in% c("ms_level",x$sample_group))]
  
  if (base::is.null(xPat)) {
    patSampleInfo <- base::data.frame(path = base::dirname(MSnbase::fileNames(x)),
                                      analysis = x$sample_name,
                                      group = x$sample_group,
                                      blank = "")
    xPat <- patRoon::importFeatureGroupsXCMS3(x, patSampleInfo)
  }
  
  patfl <- base::as.data.frame(patRoon::as.data.table(xPat, average = TRUE, areas = FALSE))
  rGroups <- base::unique(x$sample_group)
  for (i in 1:base::length(rGroups)) {
    patfl[,base::paste0(rGroups[i],"_sd")] <- base::apply(
      X = patRoon::as.data.table(xPat, average = FALSE)[, .SD, .SDcols = x$sample_name[x$sample_group == rGroups[i]]],
      MARGIN = 1, function(x) {
        base::round(base::ifelse(sd(x) != 0, sd(x)/mean(x), NA), digits = 2)})
  }
  fl <- base::cbind(fl, dplyr::select(patfl, -ret))
  fl$rt <- fl$rtmed
  fl$patFT <- fl$group
  fl <- base::cbind(dplyr::select(fl, FT, patFT, mz, rt, dplyr::everything(), -mzmin, -mzmax, -rtmin, -rtmax, -group, -mzmed, -rtmed),
                    dplyr::select(fl, mzmin, mzmax, rtmin, rtmax))
  
  fl <- dplyr::mutate(fl, hasFilled = 0,
                      sn = 0, sn_max = 0, sn_sd = 0, noise = 0, noise_sd = NA,
                      egauss = 0, egauss_sd = NA, egauss_min = 0, 
                      dppm = 0, dppm_sd = 0, dppm_max = 0,
                      others_N = 0, others = base::I(base::list("")), others_R = base::I(base::list("")),
                      bg25 = NA, bg50 = NA, bg75 = NA, bg100 = NA,
                      nCent = base::I(base::list(0)))
  
  fl2 <- fl
  
  base::row.names(fl2) <- 1:base::nrow(fl2)
  
  fl2$hasFilled <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                             function(x) {1 %in% chromPeaks[base::unlist(x), "is_filled"] }))
  
  #update both mz and rt min and max values 
  fl2$mzmin <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                         function(x) {base::min(chromPeaks[base::unlist(x), "mzmin", drop = T])}))              
  fl2$mzmax <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                         function(x) {base::max(chromPeaks[base::unlist(x), "mzmax", drop = T])})) 
  fl2$rtmin <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                         function(x) {base::min(chromPeaks[base::unlist(x), "rtmin", drop = T])})) 
  fl2$rtmax <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                         function(x) {base::max(chromPeaks[base::unlist(x), "rtmax", drop = T])})) 
  
  #add gaussian fitting values
  fl2$egauss <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                          function(x) {base::round(base::mean(chromPeaks[base::unlist(x), "egauss", drop = T], na.rm = T), digits = 2)})) 
  fl2$egauss_sd <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                             function(x) {base::round(stats::sd(chromPeaks[base::unlist(x), "egauss", drop = T], na.rm = T), digits = 2)}))
  fl2$egauss_min <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                              function(x) {
                                                x <- chromPeaks[base::unlist(x), "egauss", drop = T]
                                                x <- base::ifelse(base::all(base::is.na(x)), NA, base::round(base::min(x, na.rm = T), digits = 2))}))
  fl2$egauss[base::is.nan(fl2$egauss)] <- NA
  
  
  #add dppm
  fl2$dppm <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                        function(x) {base::round(base::mean(chromPeaks[base::unlist(x), "dppm", drop = T], na.rm = T), digits = 0)})) 
  fl2$dppm_sd <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                           function(x) {base::round(stats::sd(chromPeaks[base::unlist(x), "dppm", drop = T], na.rm = T), digits = 0)}))
  fl2$dppm_max <- base::unlist(base::lapply(X = base::as.list(fl$peakidx[drop = F]),
                                            function(x) {base::round(base::max(chromPeaks[base::unlist(x), "dppm", drop = T], na.rm = T), digits = 0)}))
  
  
  #collect bg information for the 15% quantile (probability) of the intensity centroids
  if (base::max(fl2$rt)-base::min(fl2$rt) < 500) {
    rtr <- stats::quantile(base::min(fl2$rtmin):base::max(fl2$rtmax), probs = seq(0,1,1))
  } else {
    if (base::max(fl2$rt)-base::min(fl2$rt) < 1000) {
      rtr <- stats::quantile(base::min(fl2$rtmin):base::max(fl2$rtmax), probs = seq(0,1,0.5))
    } else {
      rtr <- stats::quantile(base::min(fl2$rtmin):base::max(fl2$rtmax), probs = seq(0,1,0.25))
    }
  }
  
  rGroups <- x$sample_group
  
  base::gc(verbose = FALSE, full = TRUE)
  
  #system.time({
  
  base::cat("Gathering feature details...")
  base::cat("\n")
  
  extraInfo <- BiocParallel::bplapply(X = 1:base::nrow(fl2), cent = cent, fl2 = fl2, rtr = rtr, rGroups = rGroups, snWindow = snWindow,
                                      function(i, cent, fl2, rtr, rGroups, snWindow) {
                                        temp <- fl2[,c("sn", "sn_max", "sn_sd", "noise", "noise_sd", "others_N", "others", "others_R", "bg25", "bg50", "bg75", "bg100", "nCent")][1,]
                                        temp[,c("bg25","bg50","bg75","bg100")] <- NA
                                        
                                        temp_cent <- cent[cent$mz >= fl2$mzmin[i] & cent$mz <= fl2$mzmax[i],]
                                        
                                        #base::cat(base::paste0(i,", "))
                                        
                                        nCent <- dplyr::group_by(temp_cent, file) 
                                        nCent <- dplyr::count(nCent, file)
                                        #nCent <- cbind(data.frame(rg = rGroups),nCent)
                                        temp$nCent[1] <- I(base::list(c(base::min(nCent$n),base::max(nCent$n))))
                                        
                                        temp$bg25[1] <- base::round(stats::quantile(temp_cent$into[temp_cent$rt <= rtr[2]],
                                                                                    probs = base::seq(0,1,0.15))[2], digits = 0)
                                        
                                        if(base::max(fl2$rt)-base::min(fl2$rt) >= 500) {
                                          temp$bg50[1] <- base::round(stats::quantile(temp_cent$into[temp_cent$rt <= rtr[3] & temp_cent$rt > rtr[2]],
                                                                                      probs = base::seq(0,1,0.15))[2], digits = 0)
                                        }
                                        
                                        if(base::max(fl2$rt)-base::min(fl2$rt) >= 1000) {
                                          temp$bg75[1] <- base::round(stats::quantile(temp_cent$into[temp_cent$rt <= rtr[4]  & temp_cent$rt > rtr[3]],
                                                                                      probs = base::seq(0,1,0.15))[2], digits = 0)
                                          
                                          temp$bg100[1] <- base::round(stats::quantile(temp_cent$into[temp_cent$rt <= rtr[5]  & temp_cent$rt > rtr[4]],
                                                                                       probs = base::seq(0,1,0.15))[2], digits = 0)
                                        }
                                        
                                        #Find other peaks within the same mz space as feature
                                        otherpeaks <- fl2[fl2$mz >= fl2$mzmin[i] & fl2$mz <= fl2$mzmax[i],]
                                        otherpeaks <- otherpeaks[otherpeaks$FT != fl2$FT[i],]
                                        
                                        if (base::nrow(otherpeaks) > 0) {
                                          
                                          for (j in 1:base::nrow(otherpeaks)) { #nrow(otherpeaks)
                                            
                                            #test the resolution between peaks to verify separation using R=(rt2-rt1)/((w1+w2)/2)
                                            R <- base::abs(fl2$rt[i] - otherpeaks$rt[j]) / (((fl2$rtmax[i] - fl2$rtmin[i]) + (otherpeaks$rtmax[j] - otherpeaks$rtmin[j])) / 2 )
                                            
                                            checkOtherPeaksNoise <- FALSE
                                            if (checkOtherPeaksNoise) {
                                              # test the s/n of the otherpeak to remove it or not from temp_cent, using the replicate group with the highest intensity
                                              otherpeak_noise <- cent[cent$mz >= otherpeaks$mzmin[j] & cent$mz <= otherpeaks$mzmax[j],]
                                              otherpeak_noise <- dplyr::filter(temp_cent, dplyr::between(rt, otherpeaks$rtmin[j] - snWindow, otherpeaks$rtmax[j] + snWindow))
                                              otherpeak_int <- dplyr::select(otherpeaks[j,], dplyr::all_of(base::unique(rGroups)))
                                              otherpeak_rg <- base::which(otherpeak_int == base::max(otherpeak_int))
                                              otherpeak_int <- base::max(otherpeak_int)
                                              if (base::nrow(otherpeak_noise) > 0) {
                                                otherpeak_noise <-  base::round(stats::quantile(otherpeak_noise$into[otherpeak_noise$file %in%
                                                                                                                       base::which(rGroups == base::unique(rGroups)[otherpeak_rg])],
                                                                                                probs = base::seq(0,1,0.25))[2], digits = 0)
                                                otherpeak_noise <- base::unname(otherpeak_noise)
                                              } else {otherpeak_noise <-  0}
                                              otherpeak_check <- base::ifelse(otherpeak_noise == 0, 3, otherpeak_int/otherpeak_noise)
                                              otherpeak_check <- base::unname(otherpeak_check)
                                              
                                              if (otherpeak_check > 3 & !base::is.nan(otherpeak_check)) {
                                                #if other peak has s/n higher than 3 then is rt space is subtracted from the centroids table
                                                temp_cent <- dplyr::filter(temp_cent, !dplyr::between(rt, otherpeaks$rtmin[j], otherpeaks$rtmax[j]))
                                                if (temp$others_R[1] == "") {
                                                  temp$others_R[1] <- base::I(base::list(base::round(R, digits = 1)))
                                                } else {
                                                  temp$others_R[1] <- base::I(base::list(c(base::unlist(temp$others_R[1]), base::round(R, digits = 1))))
                                                }
                                                temp$others_N[1] <- temp$others_N[1]+1
                                                if (temp$others[1] == "") {
                                                  temp$others[1] <- base::I(base::list(otherpeaks$FT[j,drop=T]))
                                                } else {
                                                  temp$others[1] <- base::I(base::list(c(base::unlist(temp$others[1]), otherpeaks$FT[j,drop=T])))
                                                }
                                              }
                                              base::rm(otherpeak_check, otherpeak_int, otherpeak_noise, otherpeak_rg)
                                              
                                            } else {
                                              
                                              temp_cent <- dplyr::filter(temp_cent, !dplyr::between(rt, otherpeaks$rtmin[j], otherpeaks$rtmax[j]))
                                              if (temp$others_R[1] == "") {
                                                temp$others_R[1] <- I(base::list(base::round(R, digits = 1)))
                                              } else {
                                                temp$others_R[1] <- I(base::list(c(base::unlist(temp$others_R[1]), base::round(R, digits = 1))))
                                              }
                                              temp$others_N[1] <- temp$others_N[1]+1
                                              if (temp$others[1] == "") {
                                                temp$others[1] <- I(base::list(otherpeaks$FT[j,drop=T]))
                                              } else {
                                                temp$others[1] <- I(base::list(c(base::unlist(temp$others[1]), otherpeaks$FT[j,drop=T])))
                                              }
                                              
                                            }
                                          }
                                          base::rm(j)
                                        }
                                        
                                        
                                        #calculate s/n for the feature in each replicate
                                        temp_cent <- dplyr::filter(temp_cent, dplyr::between(rt, fl2$rtmin[i] - snWindow, fl2$rtmax[i] + snWindow))
                                        temp_cent <- dplyr::filter(temp_cent, !dplyr::between(rt, fl2$rtmin[i], fl2$rtmax[i]))
                                        temp_int <- dplyr::select(fl2[i,], dplyr::all_of(base::unique(rGroups)))
                                        noiselevel <- base::rep(0,base::length(temp_int))
                                        for (rg in 1:base::length(temp_int)) {
                                          if (base::nrow(temp_cent) > 0) {
                                            noiselevel[rg] <-  base::round(stats::quantile(temp_cent$into[temp_cent$file %in%
                                                                                                            base::which(rGroups == base::unique(rGroups)[rg])],
                                                                                           probs = base::seq(0,1,0.25))[3], digits = 0)
                                          } else {base::message(base::paste0("sn could not be calculated for: ",fl2$FT[i]))}
                                        }
                                        
                                        sn <- base::unlist(temp_int/noiselevel)
                                        sn[base::is.infinite(sn)] <- NA
                                        sn[base::is.nan(sn)] <- NA
                                        
                                        sn_max <- base::ifelse(TRUE %in% !base::is.na(sn), base::max(sn, na.rm = T), NA)
                                        sn_sd <- stats::sd(sn, na.rm = T)
                                        noise <- base::mean(noiselevel, na.rm = T)
                                        noise_sd <- stats::sd(noiselevel, na.rm = T)
                                        
                                        temp$sn <- base::round(base::mean(sn, na.rm = T), digits = 0)
                                        temp$sn_max <- base::round(sn_max, digits = 0)
                                        temp$sn_sd <- base::round(sn_sd, digits = 0)
                                        temp$noise <- base::round(noise, digits = 0)
                                        temp$noise_sd <- base::round(noise_sd, digits = 0)
                                        
                                        i <- temp
                                        
                                      }, BPPARAM = BiocParallel::bpparam("SnowParam"))
  
  base::gc(verbose = FALSE, full = TRUE)
  
  extraInfo <- base::do.call("rbind", extraInfo)
  
  #})
  
  fl2$sn <- extraInfo$sn
  fl2$sn_max <- extraInfo$sn_max
  fl2$sn_sd <- extraInfo$sn_sd
  fl2$noise <- extraInfo$noise
  fl2$noise_sd <- extraInfo$noise_sd
  fl2$others_N <- extraInfo$others_N
  fl2$others <- extraInfo$others
  fl2$others_R <- extraInfo$others_R
  fl2$bg25 <- extraInfo$bg25
  fl2$bg50 <- extraInfo$bg50
  fl2$bg75 <- extraInfo$bg75
  fl2$bg100 <- extraInfo$bg100
  fl2$nCent <- extraInfo$nCent
  
  #add annotation to feature list
  
  fl3 <- fl2
  
  if (!base::is.null(xA)) {
    
    fl3 <- dplyr::mutate(fl3, comp = 0, Mion = 0, iso = "", isoGroup = 0, adduct = "", adductMions = "")
    
    base::cat("Adding annotation details...")
    base::cat("\n")
    
    #produce table with annotation and featureID
    rIndex <- 1:base::length(xA)
    for (r in rIndex) {
      ano_temp <- CAMERA::getPeaklist(xA[[r]], intval = "maxo")
      ano_temp$rGroup <- base::names(xA)[r]
      ano_temp$FT <- fl3$FT
      ano_temp <- dplyr::select(ano_temp, rGroup, FT, dplyr::everything())
      if (r == 1 | base::length(rIndex) == 1) {
        ano <- ano_temp
      } else {
        ano <- base::rbind(ano,ano_temp)
      }
    }
    
    rGroups <- rGroups[rGroups %in% base::names(xA)]
    
    ano2 <- BiocParallel::bplapply(X = 1:base::nrow(fl3), ano = ano, fl3 = fl3, rGroups = rGroups, function(i, ano, fl3, rGroups) {
      
      ano_temp <- fl3[c("Mion", "comp", "isoGroup", "iso", "adduct", "adductMions")][1,]
      
      FT <- fl3$FT[i]
      
      FT_ano <- ano[ano$FT == FT,]
      
      Mion <- base::round(base::unique(FT_ano$mz) - 1.007276, digits = 4)
      
      #Save comp number
      comp <- base::as.numeric(base::unique(FT_ano$pcgroup))
      ano_temp$comp <- comp
      
      #extract iso group numbers and type for each rGroup
      isotopes <- FT_ano$isotopes
      isotopes <- base::strsplit(isotopes, split = "\\]\\[")
      
      isoGroups <- base::lapply(X = 1:base::length(isotopes), isotopes = isotopes, function(x, isotopes) isotopes[[x]][1])
      isoGroups <- stringr::str_extract(isoGroups, pattern = "([0-9]+)")
      ano_temp$isoGroup <- base::I(base::list(base::as.numeric(isoGroups)))
      
      iso <- base::lapply(X = 1:base::length(isotopes), isotopes = isotopes, function(x, isotopes) isotopes[[x]][2])
      iso <- base::ifelse(!base::is.na(iso), base::paste0("[",iso), NA)
      
      #if all the same iso type, excluding NA values, recalculates Mion
      uniqueIso <- stats::na.omit(base::unique(iso))
      ano_temp$iso <-  base::I(base::list(base::unique(iso[!base::is.na(iso)])))
      
      if (base::length(uniqueIso) == 1) {
        Mion <- base::data.frame(rGroup = base::unique(rGroups)[base::which(!base::is.na(isoGroups))])
        Mion$isoGroups <- base::as.numeric(isoGroups[base::which(!base::is.na(isoGroups))])
        Mion$Mion <- base::unlist(base::lapply(X=1:base::nrow(Mion), ano = ano, Mion = Mion, function(m, ano, Mion){
          temp_ano <- ano[ano$rGroup == Mion[m,1, drop=T],]
          temp_iso <- base::strsplit(temp_ano$isotopes, split = "\\]\\[")
          temp_iso <- base::lapply(X = 1:base::length(temp_iso), temp_iso = temp_iso, function(x, temp_iso) {temp_iso[[x]][1]})
          temp_iso <- stringr::str_extract(temp_iso, pattern = "([0-9]+)")
          temp_ano <- temp_ano[base::as.numeric(temp_iso) %in% Mion$isoGroups[m],]
          temp_ano <- temp_ano[base::grepl(pattern = "[M]", temp_ano$isotopes, fixed = TRUE),]
          temp_iso <- base::unlist(base::strsplit(temp_ano$isotopes, split = "\\]\\["))[2]
          temp_iso <- stringr::str_extract(temp_iso, pattern = "([0-9]+)")
          temp_iso <- base::ifelse(base::is.na(temp_iso), 1, temp_iso)
          temp_ano <- temp_ano$mz - 1.007276 / base::as.numeric(temp_iso)
          temp_ano <- base::round(temp_ano, digits = 5)
        }))
        Mion_temp <- base::unique(Mion$Mion)
        Mion <- base::ifelse(base::length(Mion_temp) == 1, Mion_temp, Mion)
      }
      
      ano_temp$Mion <- Mion
      
      
      #extract adducts
      adducts <- FT_ano$adduct
      adducts <- base::unlist(base::strsplit(adducts, split = " "))
      
      if (base::length(adducts) > 0) {
        
        adductClass <- adducts[base::seq(1, base::length(adducts),2)]
        adductClass <- base::unique(adductClass)
        ano_temp$adduct <- base::paste(adductClass, collapse = " ")
        AdductMions <- adducts[base::seq(2,base::length(adducts),2)]
        AdductMions <- base::as.numeric(base::unique(AdductMions))
        ano_temp$adductMions <- base::I(base::list(AdductMions))
        
      }
      
      i <- ano_temp
      
    }, BPPARAM = BiocParallel::bpparam("SnowParam"))
    
    base::gc(verbose = FALSE, full = TRUE)
    
    ano2 <- base::do.call("rbind", ano2)
    
    fl3$comp <- ano2$comp
    fl3$Mion <- ano2$Mion
    fl3$iso <- ano2$iso
    fl3$isoGroup <- ano2$isoGroup
    fl3$adduct <- ano2$adduct
    fl3$adductMions <- ano2$adductMions
    
  }
  
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(fl3, file = base::paste0(rData,"\\fl.rds"))
  }
  
  return(fl3)
  
}



#' @title screenSuspectsFromCSV_Old
#' 
#' @description Method to perform suspect screening from a given
#' list of candidates in a csv file. The method uses the S4 method
#' \code{screenSuspects} from \pkg{patRoon}.
#'
#' @param obj An \linkS4class{ntsData} object with features
#' to preform suspect screening.
#' @param suspects A \code{data.frame} or a location for a csv
#' with details for each suspect compound.
#' See details for more information about the required data.frame structure.
#' @param ppm The mass deviation, in ppm, to screen for the suspects.
#' The default is 5 ppm.
#' @param rtWindow The retention time deviation, in seconds,
#' to screen for the QC target standards. The default is 30 seconds.
#' @param adduct The adduct class for screening suspects. The default is \code{NULL},
#' which uses the polarity of the \code{obj} or the \code{adduct} class defined
#' in the suspects data frame. See details for more information.
#' @param excludeBlanks Logical, set to \code{TRUE} to ignore replicate groups 
#' assigned as blanks in the \code{samples} slot of the \code{obj}.
#' @param withMS2 Logical, set to \code{TRUE} for using confirmation via MS2.
#' @param listMS2 An MS2 list, if not given it will be calculated.
#' @param ppmMS2 Optional, sets a different mass deviation (in ppm)
#' for correlation of MS2 data.
#' If \code{NULL} (the default) the \code{ppm} argument is used.
#'
#' @return A data.frame with the suspect screening results.
#' 
#' @details The \code{suspects} data.frame should follow the template requires as 
#' obtained via \code{getScreeningListTemplate()}.
#' Add other details of the template.
#' 
#' @references \insertRef{Helmus2021}{ntsIUTA}
#' 
#' @export
#' 
#' @importFrom checkmate assertClass
#' @importFrom dplyr filter arrange left_join select everything mutate distinct desc
#' @importFrom utils read.csv head
#' @importFrom patRoon screenSuspects as.data.table screenInfo getDefAvgPListParams replicateGroups generateMSPeakLists generateFormulasSIRIUS
#' @importFrom fuzzyjoin difference_inner_join
#'
#' @examples
#' 
screenSuspectsFromCSV_Old <- function(obj,
                                  suspects = utils::read.csv(base::file.choose()),
                                  ppm = 5,
                                  rtWindow = 30,
                                  adduct = NULL,
                                  excludeBlanks = TRUE,
                                  withMS2 = TRUE, listMS2 = NULL, ppmMS2 = NULL) {
  
  # TODO Adapt to ntsData frame work
  
  assertClass(obj, "ntsData")
  
  if (!is.data.frame(suspects)) suspects <- utils::read.csv(suspects)
  
  if (max(suspects$rt) < 120) suspects$rt <- suspects$rt * 60
  
  #selects top 5 or 10 fragment if MS2 data is present for suspects
  if ("hasFragments" %in% colnames(suspects)) {
    for (n in seq_len(nrow(suspects))) {
      if (suspects$hasFragments[n]) {
        top <- unlist(suspects$fragments_mz[n])
        top <- as.data.frame(unlist(strsplit(top, split = ";")))
        if (nrow(top) > 5) {
          colnames(top) <- c("mz")
          top$mz <- as.numeric(top$mz)
          top$int <- as.numeric(unlist(strsplit(suspects$fragments_int[n], split = ";")))
          #remove precursor ion from fragments list
          top <- dplyr::filter(top, mz < suspects$mz[n] - (5 / 1E6 * suspects$mz[n]))
          if (nrow(top) < 15) ntop = 5 else ntop = 10
          top <- utils::head(dplyr::arrange(top, dplyr::desc(int)), n = ntop)
          suspects$fragments_mz[n] <- paste(top$mz, collapse = ";")
          suspects$fragments_int[n] <- paste(top$int, collapse = ";")
          # TODO add top for formulas but might not be necessary        
        }
      }
    }
  }
  
  if (is.null(adduct)) {
    if (!("adduct" %in% colnames(suspects))) {
      adduct <- ifelse(obj@polarity == "positive", "[M+H]+",
                       ifelse(obj@polarity == "negative", "[M-H]-",
                              stop("polarity argument must be 'positive' or 'negative'")))
    }
  }
  
  rg <- unique(sampleGroups(obj))
  
  if (excludeBlanks) rg <- rg[!(rg %in% blanks(obj))]
  
  screen <- screenSuspects(obj@patdata, select(suspects, -mz),
                           rtWindow = rtWindow, mzWindow = 0.03,
                           adduct = adduct, onlyHits = TRUE)
  
  df <- arrange(patRoon::as.data.frame(screen, average = FALSE), group)
  df <- left_join(df, suspects[, colnames(suspects) %in% c("name", "formula", "adduct", "hasFragments", "intControl")], by = "name")
  df <- left_join(df, select(arrange(screenInfo(screen), group), group, d_mz, d_rt), by = "group")
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- dplyr::rename(df, rt = ret, ID = group)
  df <- select(df, name, formula, adduct, ID, mz, rt, everything(), -d_mz, d_ppm, d_rt)
  df <- dplyr::filter(df, d_ppm <= ppm)
  df$d_ppm <- round(df$d_ppm, digits = 1)
  df$d_rt <- round(df$d_rt, digits = 1)
  
  screen <- screen[, df$ID]
  
  elements <- gsub("[^a-zA-Z]", "", df$formula)
  elements <- paste0(elements, collapse = "")
  elements <- gsub('([[:upper:]])', ' \\1', elements)
  elements <- unique(strsplit(elements, " ")[[1]]) 
  elements <- elements[!elements %in% ""]
  elements <- paste0(elements, collapse = "")
  
  data <- list()
  results <- list()
  
  for (g in seq_len(length(rg))) {
    
    temp <- filterFeatureGroups(screen, which(obj@samples$group == rg[g]))
    
    # TODO Adduct is [M+H]+ by default but should take the value from screening
    
    control_avgPListParams <- getDefAvgPListParams(
      clusterMzWindow = 0.003,
      topMost = 50,
      minIntensityPre = 10,
      minIntensityPost = 10
    )
    
    MS2 <- suppressWarnings(generateMSPeakLists(
      temp, "mzr",
      maxMSRtWindow = 10,
      precursorMzWindow = 1.3,
      avgFeatParams = control_avgPListParams, 
      avgFGroupParams = control_avgPListParams
    ))
    
    formulas <- patRoon::generateFormulasGenForm(temp, MS2,
                                                 relMzDev = ppm,
                                                 isolatePrec = TRUE,
                                                 adduct = "[M+H]+",
                                                 elements = elements,
                                                 topMost = 25, extraOpts = NULL,
                                                 calculateFeatures = TRUE,
                                                 featThreshold = 1,
                                                 timeout = 240, hetero = TRUE, oc = TRUE)
    
    
    temp <- patRoon::annotateSuspects(temp, MSPeakLists = MS2,
                                      formulas = formulas,
                                      compounds = NULL)
    
    
    tempScreen <- temp@screenInfo[match(tempdf$ID, temp@screenInfo$group)]
    
    tempScreen$d_ppm <- (abs(tempScreen$d_mz) / tempScreen$mz) * 1E6
    
    isoScore <- tempScreen$group
    for (iso in seq_len(nrow(tempScreen))) {
      tempForm <- formulas@formulas[[isoScore[iso]]]
      tempForm <- tempForm[tempForm$neutral_formula %in% tempScreen$formula, ]
      isoScore[iso] <- round(unique(tempForm$isoScore)[1], digits = 2)
    }
    
    tempScreen$isoScore <- isoScore
    
    tempScreen$hasExpFrag <- FALSE
    for (i in seq_len(nrow(tempScreen))) {
      fragments <- MS2[[tempScreen$group[i]]]$MSMS
      if (!is.null(fragments)) tempScreen$hasExpFrag[i] <- TRUE
    }
    
    tempScreen$FragMatch <- paste0(tempScreen$maxFragMatches, "(", tempScreen$maxFrags, ")")
    
    tempScreen <- dplyr::rename(tempScreen, ID = group, IdLevel = estIDLevel)
    
    tempScreen$hasFrag <- df$hasFragments
    
    tempScreen <- select(tempScreen, name, formula, adduct, ID, mz, rt, IdLevel, d_rt, d_ppm, isoScore, FragMatch, hasExpFrag, hasFrag, everything())
    
    tempScreen <- cbind(tempScreen, df[, colnames(df) %in% obj@samples$sample[obj@samples$group == rg[g]]])
    
    data[[rg[g]]] <- temp
    results[[rg[g]]] <- tempScreen
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  if (base::is.null(ppmMS2)) ppmMS2 = ppmWindow
  
  
  
  if (withMS2) {
    
    if(base::is.null(listMS2)) listMS2 <- base::list()
    MS2_avgPListParams <- patRoon::getDefAvgPListParams(clusterMzWindow = 0.005,
                                                        topMost = 50,
                                                        minIntensityPre = 10,
                                                        minIntensityPost = 10)
    
    # Evaluation and categorization of hits per replicate group
    for (i in 1:base::length(groupNames)) {
      
      #Collecting name and prepare table
      colMS2 <- base::paste0("ms2_",groupNames[i])
      colCat <- base::paste0("cat_",groupNames[i])
      suspectsDT[, colMS2] <- 0
      suspectsDT[, colCat] <- 0
      
      #Check for lowest categories, 4 only mass and 3 mass and rt
      temp <- suspectsDT[,c("name", "group", "mz", "ret", "d_rt", groupNames[i]), drop = F]
      suspectsDT[, colCat] <- base::ifelse(base::as.numeric(temp[,groupNames[i]]) > 0, 4, 0)
      suspectsDT[, colCat] <- base::ifelse(base::as.numeric(temp[,groupNames[i]]) > 0 & !base::is.na(temp$d_rt), 3, suspectsDT[, colCat])
      
      #Extracts the MS2 data for each non 0 feature
      if(!base::is.null(listMS2[[groupNames[i]]])){
        tempMS2 <- listMS2[[groupNames[i]]]
      } else {
        tempMS2 <- temp[temp[,groupNames[i]] > 0,]
        tempMS2 <-base::suppressWarnings(
          patRoon::generateMSPeakLists(patData[base::which(rGroups == groupNames[i]), tempMS2$group[drop = TRUE]],
                                       "mzr", maxMSRtWindow = 5,
                                       precursorMzWindow = 2,
                                       avgFeatParams = MS2_avgPListParams,
                                       avgFGroupParams = MS2_avgPListParams)
        )
        listMS2[[groupNames[i]]] <- tempMS2
      }
      
      #Loop for all the rows in suspectsDT
      for (j in 1:base::nrow(temp)) {
        
        if (base::as.numeric(temp[j, groupNames[i], drop = T]) > 0) {
          
          temp1 <- tempMS2[[temp$group[j]]]$MSMS
          
          #tentative to extract MS2 again if temp1 is NULL
          if(base::is.null(temp1)) {
            lastTemptMS2 <- base::suppressWarnings(
              patRoon::generateMSPeakLists(patData[base::which(rGroups == groupNames[i]), temp$group[j, drop = TRUE]],
                                           "mzr", maxMSRtWindow = 5,
                                           precursorMzWindow = 2,
                                           avgFeatParams = MS2_avgPListParams,
                                           avgFGroupParams = MS2_avgPListParams))
            temp1 <- lastTemptMS2[[temp$group[j]]]$MSMS
          }
          
          
          if (!base::is.null(temp1)) {
            
            temp1 <- dplyr::filter(temp1, mz < temp$mz[j]+(5/1E6*temp$mz[j])) #remove mz higher than precursor, probably contamination
            if (base::nrow(temp1) < 15){
              top5_temp1 <- utils::head(dplyr::arrange(temp1, dplyr::desc(intensity)), n = 5)
            } else { top5_temp1 <- utils::head(dplyr::arrange(temp1, dplyr::desc(intensity)), n = 10) }
            
            
            # load MS2 from DB
            temp2 <- dplyr::filter(sDB, name == temp$name[j])
            #if (base::is.na(temp2$hasMS2[drop = T])) {temp2$hasMS2 <- FALSE}
            
            if (temp2$hasMS2[drop = T]) {
              
              temp3 <- temp2$mzMS2[drop = T]
              temp3 <- base::as.data.frame(base::unlist(base::strsplit(temp3, split=";")))
              base::colnames(temp3) <- c("mz")
              temp3$mz <- base::as.numeric(temp3$mz)
              temp3$intensity <- base::as.numeric(base::unlist(base::strsplit(temp2$intMS2[drop = T], split=";")))
              #remove precurssor ion from fragments list
              temp3 <- dplyr::filter(temp3, mz < temp2$mz[drop=T] - (5/1E6*temp2$mz[drop=T]))
              # select top 5 in fragments from db, or top 10 if number of fragments is above 15
              if (base::nrow(temp3) < 15){
                top5_temp3 <- utils::head(dplyr::arrange(temp3, dplyr::desc(intensity)), n = 5)
              } else { top5_temp3 <- utils::head(dplyr::arrange(temp3, dplyr::desc(intensity)), n = 10) }
              
              
              # test match for MS2 correlation and excludes in diff > +/-5ppm
              temp6 <- fuzzyjoin::difference_inner_join(temp1, temp3, by = c("mz"), max_dist = 0.1, distance_col = "d_ppm")
              temp6$d_ppm <- base::abs(temp6$d_ppm)/temp6$mz.x*1E6
              temp6 <- dplyr::filter(temp6, d_ppm <= ppmMS2)
              
              # test match only in top 5
              top5_temp6 <- fuzzyjoin::difference_inner_join(top5_temp1, top5_temp3, by = c("mz"), max_dist = 0.1, distance_col = "d_ppm")
              top5_temp6$d_ppm <- base::abs(top5_temp6$d_ppm)/top5_temp6$mz.x*1E6
              top5_temp6 <- dplyr::filter(top5_temp6, d_ppm <= ppmMS2)
              
              temp6 <- dplyr::distinct(temp6, mz.x, .keep_all= TRUE) # remove double entries for mz.x
              top5_temp6 <- dplyr::distinct(top5_temp6, mz.x, .keep_all= TRUE) # remove double entries for mz.x
              
              suspectsDT[j, colMS2] <- base::paste0(base::nrow(top5_temp6),"(",base::nrow(top5_temp3),")")
              
              if (base::nrow(top5_temp6) >= 2) {
                
                suspectsDT[j, colCat] <- 1
                
              } else {
                
                groupAndIsos <- dplyr::filter(df_patData, mz >= temp$mz[j, drop=T] & mz < 6+temp$mz[j, drop=T])
                groupAndIsos <- dplyr::filter(groupAndIsos, ret >= temp$ret[j, drop=T]-30 & ret <= temp$ret[j, drop=T]+30)
                
                tempMS2withIsotopes <- base::suppressWarnings(
                  patRoon::generateMSPeakLists(patData[base::which(rGroups == groupNames[i]),
                                                       groupAndIsos$group[drop = TRUE]],
                                               "mzr", maxMSRtWindow = 5, precursorMzWindow = 2,
                                               avgFeatParams = MS2_avgPListParams,
                                               avgFGroupParams = MS2_avgPListParams)
                )
                
                #tentative to identify by insillico fragmentation
                tempFormulas <- patRoon::generateFormulasSIRIUS(patData[base::which(rGroups == groupNames[i]), groupAndIsos$group[drop = TRUE]],
                                                                tempMS2withIsotopes, relMzDev = ppmMS2,
                                                                adduct = "[M+H]+",
                                                                elements = base::gsub("[^a-zA-Z]", "", temp2$formula[drop=T]),
                                                                profile = "qtof",
                                                                database = NULL, noise = NULL,
                                                                topMost = 20, extraOptsGeneral = NULL,
                                                                verbose = FALSE,
                                                                calculateFeatures = TRUE, featThreshold = 1)
                
                formulaResult <- tempFormulas[[temp$group[j, drop=T]]]
                if (!base::is.null(formulaResult)) {
                  formulaResult <- dplyr::filter(formulaResult, neutral_formula == temp2$formula[drop=T])
                  suspectsDT[j, colMS2] <- base::paste0(base::nrow(formulaResult),"(",base::nrow(temp1),")")
                  if (base::nrow(formulaResult) >= 2) {
                    if (base::length(base::unique(formulaResult$frag_mz[drop = T])) >= 2) {suspectsDT[j, colCat] <- 2}
                  }
                }
              }
            } else {
              
              groupAndIsos <- dplyr::filter(df_patData, mz >= temp$mz[j, drop=T] & mz < 6+temp$mz[j, drop=T])
              groupAndIsos <- dplyr::filter(groupAndIsos, ret >= temp$ret[j, drop=T]-30 & ret <= temp$ret[j, drop=T]+30)
              
              tempMS2withIsotopes <-base::suppressWarnings(
                patRoon::generateMSPeakLists(patData[base::which(rGroups == groupNames[i]),
                                                     groupAndIsos$group[drop = TRUE]],
                                             "mzr", maxMSRtWindow = 5, precursorMzWindow = 2,
                                             avgFeatParams = MS2_avgPListParams,
                                             avgFGroupParams = MS2_avgPListParams)
              )
              
              temp2 <- dplyr::filter(sDB, name == temp$name[j])
              #tentative to identify by insillico fragmentation
              tempFormulas <- patRoon::generateFormulasSIRIUS(patData[base::which(rGroups == groupNames[i]), groupAndIsos$group[drop = TRUE]],
                                                              tempMS2withIsotopes, relMzDev = ppmMS2,
                                                              adduct = "[M+H]+",
                                                              elements = base::gsub("[^a-zA-Z]", "", temp2$formula[drop=T]),
                                                              profile = "qtof",
                                                              database = NULL, noise = NULL,
                                                              topMost = 20, extraOptsGeneral = NULL,
                                                              calculateFeatures = TRUE, featThreshold = 1)
              
              # tempFormulas <- patRoon::generateFormulasGenForm(patData[which(rGroups == groupNames[i]), groupAndIsos$group[drop = TRUE]],
              #                                                 tempMS2withIsotopes, relMzDev = 10, isolatePrec = TRUE,
              #                                                 adduct = "[M+H]+", elements = gsub("[^a-zA-Z]", "", temp2$formula[drop=T]),
              #                                                 topMost = 20, extraOpts = NULL,
              #                                                 calculateFeatures = TRUE, featThreshold = 1, timeout = 240, hetero = TRUE, oc = TRUE)
              
              formulaResult <- tempFormulas[[temp$group[j, drop=T]]]
              if (!base::is.null(formulaResult)) {
                formulaResult <- dplyr::filter(formulaResult, neutral_formula == temp2$formula[drop=T])
                suspectsDT[j, colMS2] <- base::paste0(base::nrow(formulaResult),"(",base::nrow(temp1),")")
                if (base::nrow(formulaResult) >= 2) {
                  if (base::length(base::unique(formulaResult$frag_mz[drop = T])) >= 2) {suspectsDT[j, colCat] <- 2}
                }
              }
            }
          }
        }
      }
      if (exists("temp1")) base::rm(temp1)
      if (exists("temp2")) base::rm(temp2)
      if (exists("temp3")) base::rm(temp3)
      if (exists("temp6")) base::rm(temp6)
      if (exists("top5_temp1")) base::rm(top5_temp1)
      if (exists("top5_temp3")) base::rm(top5_temp3)
      if (exists("top5_temp6")) base::rm(top5_temp6)
    }
  }
  
  suspects <- base::list(patSuspects = patSuspects, suspectsDT = suspectsDT)
  
  if (!base::is.null(listMS2)) suspects[["MS2"]] <- listMS2
  
  return(suspects)
  
  # rGroups <- patRoon::analysisInfo(patData)$group
  # #if (removeBlanks) rGroups <- rGroups[rGroups != blankGroups]
  # 
  # groupNames <- patRoon::replicateGroups(patData)
  # if (removeBlanks) groupNames <- groupNames[!(groupNames %in% blankGroups)]
  # 
  # df_patData <- base::as.data.frame(patRoon::as.data.table(patData, average = T))
  # 
  # patSuspects <- patRoon::screenSuspects(patData, sDB, rtWindow = rtWindow, mzWindow = 0.03, adduct = adduct, onlyHits = TRUE)
  # 
  # suspectsDT  <- dplyr::arrange(patRoon::screenInfo(patSuspects), group)
  # suspectsDT  <- dplyr::select(suspectsDT, group, name, d_mz, d_rt)
  # suspectsDT  <- dplyr::left_join(suspectsDT, df_patData, by = "group")
  # suspectsDT <- dplyr::left_join(suspectsDT, sDB[,c("name", "formula", "comment", "int10")], by = "name")
  # suspectsDT$d_ppm <- (base::abs(suspectsDT$d_mz)/suspectsDT$mz)*1E6
  # suspectsDT <- dplyr::select(suspectsDT, group, name, formula, d_ppm, d_rt, mz, ret, dplyr::everything(), -d_mz)
  # suspectsDT <- dplyr::filter(suspectsDT, d_ppm <= ppmWindow)
  # suspectsDT <- base::as.data.frame(suspectsDT)
  
}



#' @title checkAnnotation_Old
#' @description Extract annotation details for features selected by \code{ID},
#' \emph{m/z} and retention time or component number (\code{comp}) of all
#' or selected sample replicate groups in an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object containing annotated features.
#' @param samples The samples (name or index) to extract the features.
#' The default is \code{NULL} which considers all existing sample replicate groups,
#' excluding any assigned blank replicate groups as listed in the \code{obj}.
#' Note that the features are taken from sample replicate groups, meaning that
#' features will be extracted from sample replicate groups
#' that contain the specified samples. 
#' @param ID The identifier/s of selected features. When specified,
#' it overwrites any given \code{mz} or \code{comp} values.
#' @param mz The \emph{m/z} to find features. Not used if \code{ID} is specified.
#' Can be a vector of length 2, defining the mass range to find features.
#' @param ppm The expected mass deviation (in ppm) to search
#' for features of a given \code{mz}.
#' @param rt The expected retention time of the \emph{m/z} of interest,
#' only used if \code{ID} is not specified.
#' @param rtWindow The expected retention time deviation for searching.
#' Can be a vector of length 2, giving the time range to find features.
#' @param rtUnit The time unit used.
#' Possible values are \code{sec} and \code{min}. Default is \code{sec}.
#' @param comp The component number/s to extract features.
#' Only used if both \code{ID} and \code{mz} are \code{NULL}.
#' @param entireComponents Logical, set to \code{TRUE} (The default) to give the
#' all the features in the components represented by the selected features.
#' @param onlyAnnotated Logical, set to \code{TRUE} to return only annotated features.
#' @param onlyRelated Logical, set to \code{TRUE} to return only features that are related,
#' meaning features annotated with the same molecular ion.
#'
#' @return A \code{data.frame} with annotation details for
#' the selected/found features.
#' 
#' @note If all (\code{ID}, \code{mz} and \code{comp}) are \code{NULL},
#' the returned \code{data.frame} will contain all the features
#' for each sample replicate group. Additionally, the three logical arguments
#' are applied in the follwoing order: (1) \code{entireComponents},
#' (2) \code{onlyAnnotated} and (3) \code{onlyRelated}.
#' 
#' @export
#'
#' @importFrom checkmate assertClass assertSubset
#' @importFrom dplyr filter between
#'
checkAnnotation_Old <- function(obj = NULL,
                            samples = NULL,
                            ID = NULL,
                            mz = NULL, ppm = 5,
                            rt = NULL, rtWindow = 1, rtUnit = "sec",
                            comp = NULL,
                            entireComponents = TRUE,
                            onlyAnnotated = FALSE,
                            onlyRelated = TRUE) {
  
  assertClass(obj, "ntsData")
  
  assertSubset(rtUnit, c("sec", "min"))
  
  ft <- obj@annotation$df
  
  if (nrow(ft) == 0) return(cat("Annotation not found in the ntsData object!"))
  
  #filter for a replicate group or for all
  if (!is.null(samples)) {
    if (is.character(samples)) {
      rg <- unique(obj@samples$group[obj@samples$sample %in% samples])
    } else {
      rg <- unique(obj@samples$group[samples])
    }
    ft <- ft[ft$group %in% rg, ]
  }
  
  if (!is.null(ID)) {
    ft <- ft[ft$ID %in% ID, ]
  } else {
    #When ID is NULL but specified by mz +/- ppm  
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      ft <- dplyr::filter(ft,
                          dplyr::between(mz, mzr[1], mzr[2]),
                          dplyr::between(rt, rtr[1], rtr[2]))
    } else {
      #When only the comp number is given
      if (!is.null(comp)) ft <- ft[ft$comp %in% comp, ]
    }
  }
  
  if (nrow(ft) == 0) return(cat("Features not found with the given selection parameters!"))
  
  ft2 <- ft
  
  if (entireComponents) {
    ft2 <- obj@annotation$df[obj@annotation$df$comp %in% ft$comp, ]
  }
  
  #filter ano by selecting only annotated Features
  if (onlyAnnotated) {
    ft2 <- ft2[!is.na(ft2$Mion) | !is.na(ft2$adductMion), ]
  }
  
  if (onlyRelated) {
    ft2 <- ft2[ft2$Mion %in% na.omit(ft$Mion) | ft2$adductMion %in% na.omit(ft$adductMion), ]
  }
  
  return(ft2)
  
}


#' @title plotFeaturePeaks_Old
#' @description Plots peaks for each feature in an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object.
#' @param fileIndex The index or name of the sample/s.
#' The default is \code{NULL} and all samples are used.
#' @param ID The identifier of the features of interest.
#' When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.
#' @param mz Optional target \emph{m/z} to find features using
#' the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the features
#' when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect features.
#' The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}.
#' The default is \code{min}.
#' @param msLevel The MS level to extract the data.
#' For the moment, only 1 is possible.
#' @param names A character string with names for each feature given in \code{features}.
#' Note that length should match between \code{names} and \code{features}.
#'
#' @return A double plot with peak chromatograms on the top part
#' and feature peak groups below.
#' 
#' @export
#' 
#' @importFrom checkmate assertClass assertSubset
#' @importFrom xcms filterFile hasFeatures featureDefinitions featureChromatograms chromPeaks
#' @importMethodsFrom MSnbase rtime
#' @importFrom BiocGenerics as.data.frame
#' @importFrom patRoon as.data.table
#' @importFrom plotly plot_ly add_trace layout hide_colorbar subplot toRGB
#' @importFrom stats setNames
#'
#' @examples
#' 
setMethod("plotFeaturePeaks", "ntsData", function(obj, fileIndex = NULL,
                                                  ID = NULL,
                                                  mz = NULL, ppm = 20,
                                                  rt = NULL, rtWindow = NULL,
                                                  rtUnit = "sec",
                                                  msLevel = 1,
                                                  interactive = TRUE) {
  
  assertClass(obj, "ntsData")
  
  assertSubset(rtUnit, c("sec", "min"))
  
  assertSubset(colorBy, c("features", "samples", "samplegroups"))
  
  if (!is.null(fileIndex)) obj <- filterFileFaster(obj, fileIndex)
  
  rtr <- NULL
  
  if (!is.null(ID)) {
    ft <- obj@features[obj@features$ID %in% ID, ]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      ft <- dplyr::filter(obj@features,
                          dplyr::between(mz, mzr[1], mzr[2]),
                          dplyr::between(rt, rtr[1], rtr[2]))
    } else {
      return(cat("One of ID or mz should be given."))
    }
  }
  
  if (nrow(ft) == 0) return(cat("No features found."))
  
  pk <- list()
  for (i in seq_len(nrow(ft))) {
    pk[[ft$ID[i]]] <- obj@peaks[obj@peaks$ID %in% unlist(ft$pIdx[i]), ]
  }
  
  pk <- pk[lapply(pk, nrow) > 0]
  
  if (length(pk) == 0) return(cat("No features found."))
  
  if (is.null(rtr)) rtr <- c(min(ft$rtmin) - 60, max(ft$rtmax) + 60)
  
  if (is.null(ppm)) ppm = 5
  
  EICs <- list()
  EICs <- lapply(pk, function(x) {
    extractEIC(obj,
               mz = c(min(x$mzmin) - ((ppm / 1E6) * min(x$mzmin)),
                      max(x$mzmax) + ((ppm / 1E6) * max(x$mzmax))),
               rtWindow = rtr,
               rtUnit = "sec")
  })
  
  sp <- obj@samples$sample
  rg <- obj@samples$group
  
  
  
  
  
  colors <- getColors(nrow(ft))
  
  #first plot for samples
  spleg <- lapply(pk, function(x) x$sample)
  spleg2 <- spleg[[1]]
  if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, spleg[[s]])
  legG <- spleg2
  sleg <- !duplicated(spleg2)
  
  #Second plot for features
  FlegG <- ft$ID
  Fspleg <- lapply(pk, function(x) x$sample)
  Fspleg2 <- rep(FlegG[1], length(Fspleg[[1]]))
  if (length(Fspleg) > 1) for (s in 2:length(Fspleg)) Fspleg2 <- c(Fspleg2, rep(FlegG[s], length(Fspleg[[s]])))
  Fsleg <- !duplicated(Fspleg2)
  
  
  
  # if (colorBy == "samples") {
  #   colors <- getColors(obj, "samples")
  #   spleg <- lapply(pk, function(x) x$sample)
  #   spleg2 <- spleg[[1]]
  #   if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, spleg[[s]])
  #   legG <- spleg2
  #   sleg <- !duplicated(spleg2)
  # } else {
  #   if (colorBy == "features") {
  #     colors <- getColors(nrow(ft))
  #     legG <- ft$ID
  #     spleg <- lapply(pk, function(x) x$sample)
  #     spleg2 <- rep(legG[1], length(spleg[[1]]))
  #     if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, rep(legG[s], length(spleg[[s]])))
  #     sleg <- !duplicated(spleg2)
  #   } else {
  #     colors <- getColors(obj, "samplegroups")
  #     spleg <- lapply(pk, function(x) x$group)
  #     spleg2 <- spleg[[1]]
  #     if (length(spleg) > 1) for (s in 2:length(spleg)) spleg2 <- c(spleg2, spleg[[s]])
  #     legG <- spleg2
  #     sleg <- !duplicated(spleg2)
  #   }
  # }
  
  if (interactive) {
    
    plot <- plot_ly()
    
    counter <- 1
    
    for (i in seq_len(nrow(ft))) {
      
      for (z in seq_len(nrow(pk[[i]]))) {
        
        df <- EICs[[i]][EICs[[i]]$file == which(sp == pk[[i]]$sample[z]), ]
        
        plot <- plot %>% add_trace(df,
                                   x = df$rt,
                                   y = df$i,
                                   type = "scatter", mode = "lines",
                                   line = list(width = 0.5,
                                               color = colors[i]),
                                   connectgaps = TRUE,
                                   name = legG[z],
                                   legendgroup = legG[z],
                                   showlegend = sleg[counter]
        )
        
        df <- df[df$rt >= pk[[i]]$rtmin[z] & df$rt <= pk[[i]]$rtmax[z], ]
        
        plot <- plot %>%  add_trace(df,
                                    x = df$rt,
                                    y = df$i,
                                    type = "scatter", mode =  "lines+markers",
                                    fill = 'tozeroy', connectgaps = TRUE,
                                    fillcolor = paste(color = colors[i], 50, sep = ""),
                                    line = list(width = 0.1, color = colors[i]),
                                    marker = list(size = 3, color = colors[i]),
                                    name = legG[z],
                                    legendgroup = legG[z],
                                    showlegend = FALSE,
                                    hoverinfo = 'text',
                                    text = paste('</br> feature: ', ft$ID[i],
                                                 '</br> peak: ', pk[[i]]$ID[z],
                                                 '</br> sample: ', pk[[i]]$sample[z],
                                                 '</br> <i>m/z</i>: ', round(df$mz, digits = 4),
                                                 '</br> rt: ', round(df$rt, digits = 0),
                                                 '</br> Int: ', round(df$i, digits = 0))
        )
        
        counter <- counter + 1
        
      }
    }
    
    plot2 <- plot_ly()
    
    rect <- list()
    
    dotsColor <- c("#000000","#FDFEFE00")
    dotsColor <- stats::setNames(dotsColor, c("0","1"))
    
    for (f in seq_len(nrow(ft))) {
      
      df2 <- pk[[f]]
      
      # rect[[f]] <- list(type = "rect", fillcolor = paste(colors[f], "15", sep = ""),
      #                   line = list(color = colors[f], width = 1, dash = 'dash'),
      #                   x0 = min(df2$rtmin, na.rm = T),
      #                   x1 = max(df2$rtmax, na.rm = T),
      #                   xref = "x", y0 = 0, y1 = length(sp) - 1, yref = "y")
      
      
      
      plot2 <- plot2 %>% plotly::add_trace(df2,
                                           x = df2$rt,
                                           y = df2$sample, type = "scatter", mode = "markers",
                                           error_x = list(type = "data", symmetric = FALSE,
                                                          array = df2$rtmax - df2$rt,
                                                          arrayminus = df2$rt - df2$rtmin,
                                                          color = colors[f],
                                                          width = 5),
                                           color = as.character(df2$is_filled), colors = dotsColor, size = 6,
                                           marker = list(line = list(color = colors[f], width = 3)),
                                           name = FlegG[f],
                                           legendgroup = FlegG[f],
                                           showlegend = TRUE,
                                           hoverinfo = 'text', text = paste('</br> feature: ', names(pk[f]),
                                                                            '</br> sample: ', df2$sample,
                                                                            '</br> height: ', round(df2$intensity, digits = 0),
                                                                            '</br> width: ', round(df2$rtmax - df2$rtmin, digits = 0),
                                                                            '</br> dppm: ', round(((df2$mzmax - df2$mzmin)/df2$mz)*1E6, digits = 1),
                                                                            '</br> filled: ', base::ifelse(df2$is_filled == 1, "TRUE", "FALSE")))
    }
    #plot2 <- plot2 %>% plotly::layout(shapes = rect)
    plot2 <- plotly::hide_colorbar(plot2)
    
    plot2
    
    
    plotList <- list()
    plotList[["plot"]] <- plot
    plotList[["plot2"]] <- plot2
    
    
    xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                        titlefont = list(size = 12, color = "black"),
                        range = rtr, autotick = T, ticks = "outside")
    
    yaxis1 = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                  titlefont = list(size = 12, color = "black"))
    
    yaxis2 = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Sample",
                  titlefont = list(size = 12, color = "black"), tick0 = 0, dtick = 1)
    
    plotf <- plotly::subplot(plotList, nrows = 2, margin = 0.04, shareX = TRUE, which_layout = "merge")
    
    plotf <- plotf %>% plotly::layout(xaxis = xaxis, yaxis = yaxis1, yaxis2 = yaxis2)
    
    plotf
    
    return(plot)
    
  } else {
    
    # non-iteractive
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  #check if is x is XCMSnExp (xcms) or featureGroups (patRoon)
  if (base::class(x) == "featureGroupsXCMS3") {
    if (!base::is.null(fileIndex)) {x <- x[fileIndex,]}
    y <- x@xdata
  } else {
    if (base::class(x) == "XCMSnExp") {
      if (!base::is.null(fileIndex)) {x <- xcms::filterFile(x, fileIndex, keepFeatures = TRUE,  keepAdjustedRtime = TRUE)}
      y <- x
    } else {
      stop("x must be an XCMSnExp (xcms) or featureGroups (patRoon) object.")
    }
  }
  
  if (!xcms::hasFeatures(y)) stop ("Features not found in the given XCMSnExp object.")
  
  rtr <- c(base::min(MSnbase::rtime(y)), base::max(MSnbase::rtime(y)))
  if (!base::is.null(rt)) if (rtUnit == "min") rt <- rt*60
  if (!base::is.null(rtWindow)) if (rtUnit == "min") rtWindow <- rtWindow*60
  if (!base::is.null(rt) & !base::is.null(rtWindow)) { rtr <- c((rt) - rtWindow, (rt) + rtWindow) }
  if (base::is.null(rt)) if (!base::is.null(rtWindow)) if (base::length(rtWindow) == 2) { rtr <- rtWindow }
  
  #When priotity is for feature and than mz
  if (!base::is.null(features)) {
    if (base::unique(base::grepl("M*_R", features, fixed = FALSE))) {
      gKey <- base::cbind(patRoon::as.data.table(x, average = TRUE)[,.SD, .SDcols = "group"],
                          base::data.frame(FT = base::row.names(xcms::featureDefinitions(y))))
      gKey <- base::as.data.frame(gKey)
      FT <- gKey[gKey$group %in% features, "FT", drop = T]
    } else {
      FT <- features
    }
    #When features are no given but a specific mz +/- ppm  
  } else {
    if (base::is.null(mz)) stop("If features are not given, at least mz should be given.")
    FT <- base::row.names(xcms::featureDefinitions(y, mz = mz, ppm = ppm, rt = rtr, type = "within", msLevel = 1))
  }
  
  if (base::length(FT) == 0 | base::is.null(FT)) stop ("Features were not found with the given paramters.") 
  
  defFT <- xcms::featureDefinitions(y)
  defFT <- defFT[base::row.names(defFT) %in% FT,]
  
  #Make features legend
  features <- base::character()
  for (f in 1:base::length(FT)) {
    features <- c(features, base::paste("M", base::round(defFT$mzmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                        "_R", base::round(defFT$rtmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                        "_", FT[f],sep=""))
  }
  
  featuresID <- features
  
  if (!base::is.null(names) & base::length(names) == base::length(features)) {
    features <- names
  }
  
  chrom <- xcms::featureChromatograms(y, features = FT, aggregationFun = "sum",
                                      expandRt = rtr[2]-rtr[1], include = "any",
                                      filled = TRUE, missing = NA)
  
  df <- base::data.frame(rtime = base::as.numeric(),intensity = base::as.numeric(),sample = base::as.character(),FT = base::as.character())
  for (f in 1:base::nrow(defFT)) {
    for (s in 1:base::ncol(chrom)) {
      temp <- chrom[f,s]
      temp <- BiocGenerics::as.data.frame(temp)
      temp$sample <- y$sample_name[s]
      temp$FT <- FT[f]
      df <- base::rbind(df,temp)
    }
  }
  
  if (plotBy == "samples") {
    colors <- ntsIUTA::getColors(y, "samples")
  } else {
    if (plotBy == "features") {
      colors <- ntsIUTA::getColors(base::length(features))
    } else {
      colors = ntsIUTA::getColors(y, "groups")
    }
  }
  
  p1 <- plotly::plot_ly(df)
  showlegend = base::rep(0,base::length(features))
  for (s in 1:base::length(y$sample_name)) { 
    for (f in 1:base::length(FT)) {
      rtFT <- base::as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
      if(base::unique(TRUE %in% (rtFT$sample == s))) {
        showlegend[f] = 1 + showlegend[f]
        p1 <- p1 %>% plotly::add_trace(df,
                                       x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "rtime"],
                                       y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "intensity"],
                                       type = "scatter", mode = "lines",
                                       line = list(width = 0.5,
                                                   color = ifelse(plotBy == "samples",base::unname(colors[s]),
                                                                  ifelse(plotBy == "features", colors[f],base::unname(colors[s])))),
                                       connectgaps = TRUE,
                                       name = ifelse(plotBy == "samples",y$sample_name[s],
                                                     ifelse(plotBy == "features", features[f], y$sample_group[s])),
                                       legendgroup = ifelse(plotBy == "samples",y$sample_name[s],
                                                            ifelse(plotBy == "features", features[f], y$sample_group[s])),
                                       showlegend = ifelse(plotBy == "samples",base::ifelse(f == 1, T, F),
                                                           ifelse(plotBy == "features", base::ifelse(showlegend[f] == 1, T, F),
                                                                  base::ifelse(f == 1, T, F))))
        
        
        
        rtFT <- rtFT[rtFT$sample == s, c("rtmin","rtmax")]
        p1 <- p1 %>%  plotly::add_trace(df,
                                        x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2], "rtime"],
                                        y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2],"intensity"],
                                        type = "scatter", mode = "lines", fill = 'tozeroy', connectgaps = TRUE,
                                        fillcolor = ifelse(plotBy == "samples",paste(color = base::unname(colors[s]),50, sep = ""),
                                                           ifelse(plotBy == "features", paste(color = base::unname(colors[f]),50, sep = ""),
                                                                  paste(color = base::unname(colors[s]),50, sep = ""))),
                                        line = list(width = 0.1,
                                                    color = ifelse(plotBy == "samples",base::unname(colors[s]),
                                                                   ifelse(plotBy == "features", colors[f],base::unname(colors[s])))),
                                        name = ifelse(plotBy == "samples",y$sample_name[s],
                                                      ifelse(plotBy == "features", features[f],y$sample_group[s])),
                                        legendgroup = ifelse(plotBy == "samples",y$sample_name[s],
                                                             ifelse(plotBy == "features", features[f],y$sample_group[s])),
                                        showlegend = F) #base::ifelse(showlegend == 1, T, F)
      }
    }
  }
  
  p2 <- plotly::plot_ly()
  rect <- list()
  colorsRect <- ntsIUTA::getColors(x = base::length(features))
  dotsColor <- c("#000000","#FDFEFE00")
  dotsColor <- stats::setNames(dotsColor, c("0","1"))
  for (f in 1:base::length(FT)) {
    rtFT <- base::as.data.frame(xcms::chromPeaks(y, isFilledColumn = TRUE)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
    
    rect[[f]] <- list(type = "rect", fillcolor = base::paste(colorsRect[f],"15", sep = ""),
                      line = list(color = colorsRect[f], width = 1, dash = 'dash'), #'rgba(62, 186, 32,0.15)'
                      x0 = base::min(rtFT$rtmin, na.rm = T),
                      x1 = base::max(rtFT$rtmax, na.rm = T),
                      xref = "x", y0 = 1, y1 = base::max(base::length(y$sample_name)), yref = "y")
    
    p2 <- p2 %>% plotly::add_trace(rtFT,
                                   x = rtFT$rt,
                                   y = rtFT$sample, type = "scatter", mode = "markers",
                                   color = as.character(rtFT$is_filled), colors = dotsColor, size = 6,
                                   marker = list(line = list(color = colorsRect[f], width = 3)),
                                   showlegend = FALSE,
                                   hoverinfo = 'text', text = paste('</br> feature: ', featuresID[f],
                                                                    '</br> sample: ', y$sample_name[rtFT$sample],
                                                                    '</br> height: ', base::round(rtFT$maxo, digits = 0),
                                                                    '</br> width: ', base::round(rtFT$rtmax-rtFT$rtmin, digits = 0),
                                                                    '</br> dppm: ', rtFT$dppm,
                                                                    '</br> sn: ', rtFT$sn,
                                                                    '</br> egauss: ', base::round(rtFT$egauss, digits = 3),
                                                                    '</br> filled: ', base::ifelse(rtFT$is_filled == 1, "TRUE", "FALSE")))
  }
  p2 <- p2 %>% plotly::layout(shapes = rect)
  p2 <- plotly::hide_colorbar(p2)
  
  plotList <- list()
  plotList[[paste0("p1",s)]]<- p1
  plotList[[paste0("p2",s)]]<- p2
  
  title <- base::list(text = "Coisas", x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))
  
  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                      titlefont = base::list(size = 12, color = "black"),
                      range = c(rtr), autotick = T, ticks = "outside")
  
  yaxis1 = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                titlefont = list(size = 12, color = "black"))
  
  yaxis2 = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Sample",
                titlefont = list(size = 12, color = "black"), tick0 = 0, dtick = 1)
  
  plot <- plotly::subplot(plotList, nrows = 2, margin = 0.04, shareX = TRUE, which_layout = "merge")
  plot <- plot %>% plotly::layout(xaxis = xaxis, yaxis = yaxis1, yaxis2 = yaxis2)
  
  return(plot)
  
})




#' @title checkQC_old
#' @description Function to check QC samples in \code{rawData}.
#' 
#' @param rawQC An \linkS4class{OnDiskMSnExp} object corresponding to the QC samples.
#' The function \code{filterFile} from \code{MSnbase} can be used to subset the \code{rawData}, such as \code{MSnbase::filterFile(rawData, file = 4:6)}.
#' @param sampleInfo The \code{sampleInfo} data frame generated by \code{\link{setupProject}} function.
#' @param screeningList A \code{data.frame} with details from each standard in the QC samples.
#' See details for more information about the required data.frame structure.
#' @param paramPeaks The parameters for the choosen peak picking method.
#' See documentation of \code{\link[xcms]{chromatographic-peak-detection}} for more information. 
#' @param paramPreGrouping If \code{\link[xcms]{adjustRtime}} is preformed with the method \code{PeakGroups},
#' a pre-grouping of peaks is required for alignment. The \code{paramPreGrouping} is the parameters obtained by the selected grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} for more information.
#' @param paramAlignment The parameters for the choosen alignment method. See documentation of \code{\link[xcms]{adjustRtime}} for more information.
#' @param paramGrouping The parameters for the choosen grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} for more information.
#' @param ppmForFillingGroups The mass (in ppm) to expand the \emph{mz} for filling missing peaks in imcomplete features.
#' See \code{\link[xcms]{fillChromPeaks}} for more information.
#' @param rtWindow The retention time deviation, in seconds, allowed to screen for QC reference standards.
#' @param ppmWindow The mass deviation, in ppm, allowed to screen for QC reference standards.
#' @param polarity The acquisition polarity of the QC samples. Possible values are \code{positive} or \code{negative}.
#' @param plot Logical, set to \code{TRUE} for plotting the results.
#' @param save Logical, set to \code{TRUE} for storing  the generated QC object in the rData folder. 
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#'
#' @details The \code{screeningList} template can be obtained as csv via \code{ntsIUTA::getScreeningListTemplate()}.
#' Add other details of the template.
#'
#' @return A \code{list} containing the plots and a data.frame with the summary of the
#' retention time, mass and intensity deviations as well as the quality of the MS2 data.
#' 
#' @export
#'
#' @import magrittr
#' @importFrom utils read.csv write.csv
#' @importFrom patRoon screenSuspects screenInfo getDefAvgPListParams generateMSPeakLists as.data.table
#' @importFrom dplyr select arrange left_join everything filter mutate top_n
#' @importFrom fuzzyjoin difference_inner_join
#' @importFrom stats cor
#' @importFrom gridExtra arrangeGrob
#' @import ggplot2
#' @importFrom htmlwidgets saveWidget
#' @importFrom plotly partial_bundle
#' @importFrom xcms applyAdjustedRtime
#'
#'
#' @examples
#' 
#' 
#' 
checkQC_old <- function(rawQC = rawData,
                        sampleInfo = setup$sampleInfo,
                        screeningList = utils::choose.files(base::getwd(), "Select the QC screening list"),
                        paramPeaks = NULL,
                        paramPreGrouping = NULL,
                        paramAlignment = NULL,
                        paramGrouping = NULL,
                        ppmForFillingGroups = 5,
                        rtWindow = 30,
                        ppmWindow = 15,
                        polarity = "positive",
                        plot = TRUE,
                        save = FALSE, projPath = setup$projPath) {
  
  # rawQC = MSnbase::filterFile(rawData, file = 4:6)
  # sampleInfo = setup$sampleInfo
  # screeningList = base::paste0(system.file(package = "ntsIUTA", dir = "extdata"),"/QC_ScreeningList_ntsIUTA_MS2_pos.csv")
  # paramPeaks = param
  # paramPreGrouping = param1
  # paramAlignment = param2
  # paramGrouping = param3
  # ppmForFillingGroups = 5
  # rtWindow = 30
  # ppmWindow = 15
  # polarity = "positive"
  # plot = TRUE
  # save = TRUE
  
  if (base::is.null(rawQC)) stop("rawData of QC samples should be provided using the argument rawQC. See ?checkQC for more information.")
  
  
  qcPeaks <- ntsIUTA::peakPicking(rawQC,
                                  param = paramPeaks,
                                  removeQC = FALSE,
                                  refinePeaks = FALSE,
                                  save = FALSE)
  
  qcFeat <- ntsIUTA::makeFeatures(peaksData = qcPeaks,
                                  paramPreGrouping = paramPreGrouping,
                                  paramAlignment = paramAlignment,
                                  paramGrouping = paramGrouping,
                                  ppmForFillingGroups = ppmForFillingGroups,
                                  save = FALSE)
  
  if (plot) {
    alignPlot <- ntsIUTA::plotAlignment(qcFeat)
  }
  
  qcFeat <- xcms::applyAdjustedRtime(qcFeat)
  
  qcPat <- ntsIUTA::getPatData(qcFeat, sampleInfo = sampleInfo[sampleInfo$sample %in% qcFeat$sample_name,])
  
  #load or read screeningList
  if (!base::is.data.frame(screeningList)) {
    sl <- utils::read.csv(screeningList)
  } else {
    sl <- screeningList
  }
  if (base::max(sl$rt) < 120) sl$rt <- sl$rt * 60
  
  #check polarity
  if (base::length(polarity) == 1 & !base::is.null(polarity)) {
    adduct <- base::ifelse(polarity == "positive", "[M+H]+", base::ifelse(polarity == "negative", "[M-H]-", stop("polarity argument must be 'positive' or 'negative'")))
  } else {
    stop("polarity argument must be 'positive' or 'negative'")
  }
  
  qcID <- patRoon::screenSuspects(qcPat, dplyr::select(sl, -mz), rtWindow = rtWindow, mzWindow = 0.01, adduct = adduct, onlyHits = TRUE)
  qcdf <- dplyr::arrange(patRoon::as.data.table(qcID, average = FALSE), group)
  qcdf <- dplyr::left_join(qcdf, sl[,base::c("name", "formula", "int10", "hasMS2", "mzMS2", "intMS2", "preMS2")], by = "name")
  qcdf  <- dplyr::left_join(qcdf, dplyr::select(dplyr::arrange(patRoon::screenInfo(qcID), group), group, d_mz, d_rt), by = "group")
  qcdf$av_into <- base::rowMeans(dplyr::select(qcdf, rawQC$sample_name))
  qcdf <- qcdf %>% dplyr::mutate(sd_into = base::apply(dplyr::select(., rawQC$sample_name), 1, sd))
  qcdf <- dplyr::mutate(qcdf, sd_intop = sd_into/av_into*100)
  qcdf$d_ppm <- (base::abs(qcdf$d_mz)/qcdf$mz)*1E6
  qcdf <- dplyr::select(qcdf, name, formula, d_ppm, d_rt, mz, ret, group, dplyr::everything(), -d_mz)
  qcdf<- dplyr::filter(qcdf, d_ppm <= ppmWindow)
  qcdf <- base::as.data.frame(qcdf)
  
  
  control_avgPListParams <- patRoon::getDefAvgPListParams(
    clusterMzWindow = 0.005,
    topMost = 50,
    minIntensityPre = 10,
    minIntensityPost = 10
  )
  
  MS2 <- base::suppressWarnings(patRoon::generateMSPeakLists(
    qcID, "mzr",
    maxMSRtWindow = 5,
    precursorMzWindow = 3,
    avgFeatParams = control_avgPListParams, 
    avgFGroupParams = control_avgPListParams
  ))
  
  
  #add MS2 to screening list, adds intensity at 10ng/ml to respective column and adds mz corresponding to the polarity for MS2 matching
  # for (i in 1:base::nrow(qcdf)) {
  #   xgroup <- qcdf$group[i]
  #   xfrag <- MS2[[xgroup]]$MSMS
  # 
  #   if (!base::is.null(xfrag)) {
  #     sl$hasMS2[sl$name %in% qcdf$name[i]] <- TRUE
  #     sl$mzMS2[sl$name %in% qcdf$name[i]] <- base::paste(xfrag$mz, collapse = ";")
  #     sl$intMS2[sl$name %in% qcdf$name[i]] <- base::paste(xfrag$intensity, collapse = ";")
  #     sl$preMS2[sl$name %in% qcdf$name[i]] <- base::paste(xfrag$precursor, collapse = ";")
  #   }
  # 
  #   sl$int10[sl$name %in% qcdf$name[i]] <- qcdf$av_into[i]
  #   sl$rt[sl$name %in% qcdf$name[i]] <- qcdf$ret[i]/60
  #   sl$mz[sl$name %in% qcdf$name[i]] <- sl$neutralMass[sl$name %in% qcdf$name[i]] +  base::ifelse(polarity == "positive", 1.007276, -1.007276)
  # }
  # 
  # utils::write.csv(sl, file = base::paste0(projPath,"/ScreeningList_QC_ntsIUTA_MS2_pos.csv"))
  
  qcdf$nfrag <- 0
  qcdf$pfrag <- 0
  qcdf$intoscore <- 0
  
  for (i in 1:base::nrow(qcdf)) {
    xgroup <- qcdf$group[i]
    xname <- qcdf$name[i]
    if (!base::is.na(xgroup) & qcdf$hasMS2[i]) {
      xMS2 <- MS2[[xgroup]]$MSMS
      if (!base::is.null(xMS2)) {
        dbMS2 <- base::data.frame(mz = base::as.numeric(base::unlist(base::strsplit(qcdf$mzMS2[i], split=";"))),
                                  intensity = base::as.numeric(base::unlist(base::strsplit(qcdf$intMS2[i], split=";"))),
                                  precursor = base::as.logical(base::unlist(base::strsplit(qcdf$preMS2[i], split=";"))))
        
        xMS2 <- dplyr::top_n(xMS2, 10, intensity)
        xMS2 <- dplyr::mutate(xMS2, into_ind = intensity/base::max(xMS2$intensity))
        
        dbMS2 <- dplyr::top_n(dbMS2, 10, intensity)
        dbMS2 <- dplyr::mutate(dbMS2, into_ind = intensity/base::max(dbMS2$intensity))
        
        combi <- fuzzyjoin::difference_inner_join(xMS2, dbMS2, by = c("mz"), max_dist = 0.005, distance_col = "diff")
        
        qcdf$nfrag[i] <- base::nrow(combi)
        qcdf$pfrag[i] <- base::nrow(combi)/base::nrow(dbMS2)
        
        if (qcdf$nfrag[i] == 1) { qcdf$intoscore[i] <- 1 }
        else {
          qcdf$intoscore[i] <- stats::cor(combi$into_ind.x, combi$into_ind.y, use = "everything", method = "pearson")
        }
      } else {
        qcdf$nfrag[i] <- 0
        qcdf$pfrag[i] <- 0
        qcdf$intoscore[i] <- 0
      }
    }
  }
  
  qcdf$pfrag <- base::round(qcdf$pfrag, digits = 2)
  qcdf$intoscore <- base::round(base::as.numeric(qcdf$intoscore), digits = 4)
  
  #TODO Implement isotope check
  
  QC <- list()
  QC[["df"]] <- dplyr::select(qcdf, -hasMS2, -mzMS2, -intMS2, -preMS2)
  
  if (plot) {
    evalPlot <- gridExtra::arrangeGrob(
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = -10, ymax = 10, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -15, ymax = -10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = 15, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = -15, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 15, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = d_rt)) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_text(size = 7),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("RT diff (sec)") +
        ggplot2::ylim(-rtWindow,rtWindow) +
        ggplot2::coord_flip(),
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 5, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 5, ymax = 10, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 10, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = d_ppm)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("m/z diff (ppm)") +
        ggplot2::ylim(0,ppmWindow) +
        ggplot2::coord_flip(), 
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = 2, ymax = 5, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 5, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 2, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = nfrag)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("Nr Shared Top10 MS2") +
        #ylim(0,10) +
        ggplot2::scale_y_continuous(limits = c(0, 10), breaks = base::seq(0, 10, by = 2)) +
        ggplot2::coord_flip(),
      ggplot2::ggplot(qcdf) +
        ggplot2::theme_bw() +
        ggplot2::geom_rect(ggplot2::aes(ymin = 0.95, ymax = Inf, xmin = -Inf, xmax = Inf), fill = "ForestGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = 0.9, ymax = 0.95, xmin = -Inf, xmax = Inf), fill = "PaleGreen", alpha = 0.05) +
        ggplot2::geom_rect(ggplot2::aes(ymin = -Inf, ymax = 0.9, xmin = -Inf, xmax = Inf), fill = "white", alpha = 0.01) +
        ggplot2::geom_point(stat = "identity",
                            ggplot2::aes(x = name,
                                         y = intoscore)) +
        ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                       axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_text(size = 7),
                       axis.title.x = ggplot2::element_text(size = 7)) +
        ggplot2::ylab("MS2 Intensity Corr.") +
        ggplot2::ylim(0.7,1) +
        ggplot2::coord_flip(),
      ncol = 4, widths = base::c(8,4,4,4))
    
    
    featPlot <- ntsIUTA::plotFeaturePeaks(qcPat, fileIndex = NULL,
                                          features = qcdf$group,
                                          mz = NULL, ppm = NULL,
                                          rt = NULL, rtWindow = NULL,
                                          rtUnit = "min",
                                          plotBy = "features",
                                          names = qcdf$name)
    
    if (save) {
      results <- base::paste0(projPath,"\\results")
      if (!base::dir.exists(results)) base::dir.create(results)
      ggplot2::ggsave(base::paste0(projPath,"/results/QC_Deviations.tiff"),
                      plot = evalPlot, device = "tiff", path = NULL, scale = 1,
                      width = 17, height = 10, units = "cm", dpi = 300, limitsize = TRUE)
      htmlwidgets::saveWidget(plotly::partial_bundle(alignPlot), file = base::paste0(projPath,"/results/QC_Alignment.html"))
      htmlwidgets::saveWidget(plotly::partial_bundle(featPlot), file = base::paste0(projPath,"/results/QC_Features.html"))
    }
    
    QC[["alignPlot"]] <- alignPlot
    QC[["featPlot"]] <- featPlot
    QC[["evalPlot"]] <- evalPlot
    
  }
  
  if (save) {
    results <- base::paste0(projPath,"\\results")
    if (!base::dir.exists(results)) base::dir.create(results)
    utils::write.csv(dplyr::select(qcdf, -hasMS2, -mzMS2, -intMS2, -preMS2), file = base::paste0(projPath,"/results/QC_CheckResults.csv"))
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(QC, file = base::paste0(rData,"\\QC.rds"))
  }
  
  return(QC)
  
}




#' @title makeFeatureComponents_Old
#' @description Uses the \pkg{CAMERA} package to find isotopes by grouping features over the retention time
#' and extracted ion chromotogram (EIC) for the given deviations. The isotopes are found for each replicate group.
#'
#' @param featData An \linkS4class{XCMSnExp} object containing features.
#' @param polarity The polarity of the replicate groups. Possible values are \code{positive} or \code{negative}.
#' @param sigma The multiplier of the standard deviation for grouping features by retention time.
#' @param perfwhm Percentage of the width of the FWHM.
#' @param cor_eic_th Minimum correlation index (from 0 to 1) for the EICs of features within the same sample.  
#' @param cor_exp_th Minimum correlation index (from 0 to 1) for the EICs of features across samples within the replicate group.
#' @param pval p-value threshold for testing correlation of significance.
#' @param calcCaS Logical, set to \code{TRUE} to calculate correlation accross samples. Deafault is \code{TRUE}.
#' @param calcIso Logical, set to \code{TRUE} to include isotope detection informationen for graph clustering. Deafault is \code{FALSE}.
#' @param ppmIsotopes The expected mass deviation to find isotopes.
#' @param mzabs The expected deviation of the \emph{m/z} to find isotopes.
#' @param noise The extimated intensity threshold for the noise level used for find isotopes.
#' @param validateIsotopePatterns Logical, set to \code{TRUE} for validating the annoatated isotopes with the kegg database. 
#' @param searchAdducts Logical, set to \code{TRUE} to screen for adducts after finding isotopes.
#' @param ppmAdducts The expected mass deviation to find adducts.
#' @param extendedList Logical, set to \code{TRUE} to use the extended list of adducts. The default (\code{FALSE} uses a shorter list.)
#' @param excludeBlanks Set to \code{TRUE} to not screen blank samples for isotopes and adducts.
#' This option is intended for saving processing time as the replicate blank groupes are expected to be removed before final feature list creation.
#' @param blankGroups A character vector with the name of the blank replicate groups in the given \code{featData} object.
#' For simplification use \code{base::unique(sampleInfo$blank)}.
#' @param save Logical, set to \code{TRUE} to save the generated \code{list} of \linkS4class{xsAnnotate} objects in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#'
#' @return A \code{list} with an \linkS4class{xsAnnotate} object per replicate group as defined in the given \linkS4class{XCMSnExp} object.
#' 
#' @references
#' \insertRef{CAMERA}{ntsIUTA}
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#' 
#' @importFrom methods as
#' @importFrom xcms filterMsLevel sampnames
#' @importFrom Biobase pData
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr findAdducts
#' @importFrom utils txtProgressBar setTxtProgressBar read.table
#'
#' @examples
#' 
#' 
#' 
makeFeatureComponents_Old <- function(featData = featData, polarity = "positive",
                                  sigma = 5, perfwhm = 0.5, cor_eic_th = 0.85,
                                  cor_exp_th = 0.75, pval = 0.05,
                                  calcCaS = TRUE, calcIso = TRUE,
                                  ppmIsotopes = 50, mzabs = 0.01, noise = 350,
                                  validateIsotopePatterns = TRUE,
                                  searchAdducts = TRUE, ppmAdducts = 5, extendedList = FALSE,
                                  excludeBlanks = TRUE, blankGroups = "Blank",
                                  save = TRUE, projPath = setup$projPath) {
  
  #Examples
  # setup <- ntsIUTA::makeSetup(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)
  # setup$sampleInfo[1:3,"group"] <- "Sample"
  # rawDataExample <- ntsIUTA::importRawData(setup$sampleInfo[1:3,], save = FALSE, centroidedData = TRUE)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = param, save = FALSE)
  # paramList <- ntsIUTA::paramListExample
  # instParam <- ntsIUTA::getInstParam(paramList)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = instParam$PP, save = FALSE)
  # featDataExample <- ntsIUTA::makeFeatures(peaksDataExample, paramPreGrouping = instParam$preGrouping, paramAlignment = instParam$alignment, paramGrouping = instParam$grouping, save = FALSE)
  # featCompExample <- ntsIUTA::makeFeatureComponents(featData = featDataExample)
  
  
  # featData = featData
  # polarity = "positive"
  # sigma = 5
  # perfwhm = 0.45
  # cor_eic_th = 0.85
  # cor_exp_th = 0.85
  # pval = 0.05
  # ppmIsotopes = 50
  # mzabs = 0.01
  # noise = 350
  # searchAdducts = TRUE
  # ppmAdducts = 5
  # extendedList = FALSE
  # excludeBlanks = FALSE
  # calcCaS = TRUE
  # calcIso = TRUE
  
  #convert featData to xcmsSet class
  xSet <- methods::as(xcms::filterMsLevel(featData,  msLevel. = 1), "xcmsSet")
  xcms::sampnames(xSet) <- Biobase::pData(featData)$sample_name
  xcms::sampclass(xSet) <- Biobase::pData(featData)$sample_group
  #xSet <- xcms::fillPeaks(xSet) #not used as filling was already performed
  
  #Replicate sample group divider
  groups <- base::unique(featData$sample_group)
  
  
  if (excludeBlanks) groups <- groups[groups != blankGroups]
  
  #Holder for isotopes from each replicate group as list() and indices for which replicate group
  featComp <- base::list()
  
  xA <- CAMERA::xsAnnotate(xs = xSet, sample = c(1:base::length(featData$sample_group)),
                           polarity = polarity)
  
  xA <- CAMERA::groupFWHM(xA, sigma = sigma, #the multiplier of the standard deviation
                          perfwhm = perfwhm, #percentage of the width of the FWHM
                          intval = "maxo")
  
  xA <- CAMERA::groupCorr(xA, cor_eic_th = cor_eic_th, # Correlation threshold for EIC correlation
                          pval = pval, # p-value threshold for testing correlation of signiﬁcance
                          graphMethod = "hcs",
                          calcIso = calcIso, calcCiS = TRUE, calcCaS = calcCaS, psg_list = NULL, xraw = NULL,
                          cor_exp_th = cor_exp_th, # Threshold for intensity correlations across samples
                          intval = "maxo") #suppressWarnings()
  
  # test <- CAMERA::getPeaklist(xA)
  # test <- dplyr::filter(test, pcgroup == 8)
  # #testx <- as.data.frame(xcms::featureDefinitions(featData))
  # #testX <-base::row.names(xcms::featureDefinitions(featData))
  
  
  #Make peaksData by replicate groups as defined in the setup experiment
  pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, char = "=", width = 80, style = 3)
  
  for (rgidx in 1:base::length(groups)) {
    
    utils::setTxtProgressBar(pb, ((rgidx/base::length(groups))*100))
    
    sampleidxs <- base::which(featData$sample_group == groups[rgidx])
    
    xA_temp <- ntsIUTA::FindIsotopesWithValidationAltered(object = xA,
                                                          featData = featData,
                                                          sampleidxs = sampleidxs,
                                                          ppm = ppmIsotopes,
                                                          mzabs = mzabs,
                                                          noise = noise,
                                                          maxcharge = 3, #maxcharge set to match small molecules and lipids
                                                          intval = "maxo",
                                                          validateIsotopePatterns = validateIsotopePatterns)
    
    if (searchAdducts)
    {
      if (extendedList) {
        rules_pos <- base::system.file('rules/extended_adducts_pos.csv', package = "CAMERA")
        rules_neg <- base::system.file('rules/extended_adducts_neg.csv', package = "CAMERA")
      } else {
        rules_pos <- base::system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
        rules_neg <- base::system.file('rules/primary_adducts_neg.csv', package = "CAMERA")
      }
      
      
      if (polarity == "positive")
      {rules <- utils::read.table(rules_pos, header = TRUE, sep = ",")}
      
      if (polarity == "negative")
      {rules <- utils::read.table(rules_neg, header = TRUE, sep = ",")}
      
      
      xA_temp <- CAMERA::findAdducts(xA_temp,
                                     ppm = ppmAdducts,
                                     mzabs = 0,
                                     multiplier = 2, # highest number(n) of allowed clusterion [nM+ion]
                                     polarity = polarity,
                                     rules = rules,
                                     max_peaks = 100, # If run in parralel mode, this number deﬁnes how much peaks will be calculated in every thread
                                     psg_list = NULL) # Vector of pseudospectra indices. The correlation analysis will be only done for those groups
    }
    
    featComp[[groups[rgidx]]] <- xA_temp
    
  }
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(featComp, file = base::paste0(rData,"\\featComp.rds"))
  }
  
  return(featComp)
  
}




#' @title plotComponentSpectrum_Old
#' @description Plots the spectra for given features, \emph{m/z} and rt or components from a list of \linkS4class{xsAnnotate} objects.
#' 
#' 
#' @param xA The list of \linkS4class{xsAnnotate} objects obtained via \code{\link{makeFeatureComponents}}.
#' @param replicateGroups The replicate groups name or number to filter the \code{xA} list,
#' which is named according to the experimental replicate groups as defined in \code{\link{setupProject}}.
#' @param features Feature identifier/s.
#' @param featData The feature data to get the feature identifiers.
#' Can be a \linkS4class{XCMSnExp} or \linkS4class{featureGroups} object from \pkg{xcms} or \pkg{patRoon}, respectively.
#' @param mz The \emph{m/z} of interest. If \code{features} are specified \code{mz is not used}.
#' @param ppm The expected mass deviation to search for features of a given \code{mz}.
#' @param mzWindowPlot The \emph{m/z} range for the plot.
#' @param rt The expected retention time of the \emph{m/z} of interest, only used if \code{features} are not speficied.
#' @param rtWindow The expected retention time deviation for searching.
#' @param rtUnit The time unit used. Default is minutes.
#' @param comp The component numbers to plot. Only used if both \code{features} and \code{mz} are \code{NULL} (i.e. not specified).
#' @param onlyAnnotated Plots only annotated features.
#' @param onlyRelated Plots only relevant features (i.e. isotopes and addcuts) of the specified features. 
#' @param intval The intensity value type. Default is "maxo", corresponding to the height.
#' @param log Logical, set to \code{TRUE} (the default) to plot the intensities in a log scale.
#'
#' @return A spectrum plot of given components or features.
#' 
#' @export
#' 
#' @importFrom patRoon as.data.table
#' @importFrom xcms featureSummary featureDefinitions
#' @importFrom CAMERA getPeaklist
#' @importFrom dplyr select everything between filter all_of
#' @importFrom stringr str_extract
#' @importFrom plotly plot_ly add_bars toRGB layout
#'
#' @examples
#' 
#' 
#' 
plotComponentSpectrum_Old <- function(xA = featComp, replicateGroups = NULL,
                                  features = NULL, featData = featData,
                                  mz = NULL, ppm = 5, mzWindowPlot = NULL,
                                  rt = NULL, rtWindow = 1, rtUnit = "min",
                                  comp = NULL,
                                  onlyAnnotated = FALSE, onlyRelated = TRUE,
                                  intval = "maxo",
                                  log = TRUE) {
  
  # featData = featData
  # xA = featComp
  # replicateGroups = 2
  # features = NULL #"FT0107"
  # comp = NULL
  # onlyAnnotated = FALSE
  # onlyRelated = TRUE
  # mz = 233.0249
  # ppm = 5
  # rt = 15.7
  # rtWindow = 1
  # mzWindowPlot = c(200,300)
  # rtUnit = "min"
  # intval = "maxo"
  # log = TRUE
  
  
  #library(magrittr)
  
  #filter for a replicate group or for all
  fileIndex <- 1:base::length(xA[[1]]@xcmsSet$sample_group)
  rGroups <- base::names(xA)
  rIndex <- 1:base::length(xA)
  if (!is.null(replicateGroups)){
    if (is.character(replicateGroups)) {
      fileIndex <- base::which(xA[[1]]@xcmsSet$sample_group %in% rGroups[rGroups == replicateGroups])
      sampleNames <-  xA[[1]]@xcmsSet$sample_name[xA[[1]]@xcmsSet$sample_group %in% rGroups[rGroups == replicateGroups]]
      groupNames <- xA[[1]]@xcmsSet$sample_group[xA[[1]]@xcmsSet$sample_group %in% rGroups[rGroups == replicateGroups]]
      rIndex <- base::which(rGroups %in% replicateGroups)
    } else {
      fileIndex <- base::which(xA[[1]]@xcmsSet$sample_group %in% rGroups[replicateGroups])
      sampleNames <-  xA[[1]]@xcmsSet$sample_name[xA[[1]]@xcmsSet$sample_group %in% rGroups[replicateGroups]]
      groupNames <- xA[[1]]@xcmsSet$sample_group[xA[[1]]@xcmsSet$sample_group %in% rGroups[replicateGroups]]
      rIndex <- replicateGroups
    }
  }
  
  #When feature ID is given to look for spectra
  FT <- NULL
  if (!base::is.null(features))
  {
    if (base::unique(base::grepl("M*_R", features, fixed = FALSE))) {
      x <- featData
      featData <- featData@xdata
      gKey <- base::cbind(patRoon::as.data.table(x, average = TRUE)[,.SD, .SDcols = "group"],
                          base::data.frame(FT = base::row.names(xcms::featureSummary(featData))))
      gKey <- base::as.data.frame(gKey)
      FT <- gKey[gKey$group %in% features, "FT", drop = T]
    } else {
      FT <- features
    }
    for (r in rIndex) {
      ft_temp <- dplyr::select(base::as.data.frame(xcms::featureDefinitions(featData)), "mzmed")
      ft_temp$FT <- base::row.names(ft_temp)
      ft_temp$rGroup <- rGroups[r]
      ano_temp <- CAMERA::getPeaklist(xA[[r]], intval = intval)
      ano_temp <- base::cbind(ft_temp[,c("rGroup","FT"), drop=F],ano_temp)
      if (r == 1) {
        ano <- ano_temp
      } else {
        ano <- base::rbind(ano,ano_temp)
      }
    }
    #get comp of features
    pcgroups <- as.numeric(ano[ano$FT %in% FT,"pcgroup", drop = T])
    
  } else {
    #When features are not given but a specific mz +/- ppm  
    if (!is.null(mz)) {
      
      if (!base::is.null(rt)) if (rtUnit == "min") rt <- rt*60
      if (!base::is.null(rtWindow)) if (rtUnit == "min") rtWindow <- rtWindow*60
      if (!base::is.null(rt) & !base::is.null(rtWindow) & (base::length(rtWindow) == 1)) { rtr <- c(rt-rtWindow, rt+rtWindow) }
      if (base::is.null(rt)) if (!base::is.null(rtWindow)) if (base::length(rtWindow) == 2) { rtr <- rtWindow }
      if (!base::is.null(rt) & base::is.null(rtWindow)) {rtr <- c(rt-10, rt+10)}
      
      if (length(mz) == 1) { mzr <- c(mz - ((ppm/1E6)*mz), mz + ((ppm/1E6)*mz)) }
      if (length(mz) == 2) { mzr <- c(mz[1], mz[2]) }
      
      for (r in rIndex) {
        ano_temp <- CAMERA::getPeaklist(xA[[r]], intval = intval)
        ano_temp$rGroup <- rGroups[r]
        ano_temp <- dplyr::select(ano_temp, rGroup, dplyr::everything())
        if (r == 1 | length(rIndex) == 1) {
          ano <- ano_temp
        } else {
          ano <- base::rbind(ano,ano_temp)
        }
      }
      
      pcgroups <- ano %>% dplyr::filter(dplyr::between(rt, rtr[1],rtr[2]))
      pcgroups <- pcgroups %>% dplyr::filter(dplyr::between(mz, mzr[1],mzr[2]))
      pcgroups <- as.numeric(pcgroups$pcgroup[drop=F])
      
    } else {
      
      #When only the comp number is given
      if (base::is.null(comp)) stop("At least one of the three, Feature ID, m/z and RT range or component number (comp) should be given.")
      
      for (r in rIndex) {
        ano_temp <- CAMERA::getPeaklist(xA[[r]], intval = intval)
        ano_temp$rGroup <- rGroups[r]
        ano_temp <- dplyr::select(ano_temp, rGroup, dplyr::everything())
        if (r == 1) {
          ano <- ano_temp
        } else {
          ano <- base::rbind(ano,ano_temp)
        }
      }
      
      pcgroups <- as.numeric(comp)
      
    }
  }
  
  #filter by relevent pcgroups
  ano <- ano[ano$pcgroup %in% pcgroups,]
  
  #filter by showing only related features (i.e. isotopes and adducts)
  if (onlyRelated) {
    #when features are given
    if (!base::is.null(FT)) {
      temp_FT <- ano[ano$FT %in% FT,, drop = F]
      iso_group <- base::gsub("\\D", "", temp_FT$isotopes[drop = T])
      iso_group <- iso_group[iso_group != ""]
      adduct_Mion <- stringr::str_extract(temp_FT$adduct[drop = T], "[0-9]+\\.[0-9]+")
      adduct_Mion <- adduct_Mion[!is.na(adduct_Mion)]
      adduct_temp <- ano[base::grepl(paste(adduct_Mion, collapse = "|"), ano$adduct[drop = T]),]
      iso_adduct <- base::gsub("\\D", "", adduct_temp$isotopes[drop = T])
      iso_adduct <- iso_adduct[iso_adduct != ""]
      iso_group <- base::unique(c(iso_group,iso_adduct))
      #leave only iso_group features
      rel_temp <- ano[base::grepl(base::paste(iso_group, collapse = "|"), ano$isotopes[drop = T])|
                        base::grepl(base::paste(adduct_Mion, collapse = "|"), ano$adduct[drop = T]),]
    } else {
      if (!is.null(mz)) {
        temp_FT <- ano %>% dplyr::filter(dplyr::between(rt, rtr[1],rtr[2]))
        temp_FT <- temp_FT %>% dplyr::filter(dplyr::between(mz, mzr[1],mzr[2]))
        iso_group <- base::gsub("\\D", "", temp_FT$isotopes[drop = T])
        iso_group <- iso_group[iso_group != ""]
        adduct_Mion <- stringr::str_extract(temp_FT$adduct[drop = T], "[0-9]+\\.[0-9]+")
        adduct_Mion <- adduct_Mion[!is.na(adduct_Mion)]
        adduct_temp <- ano[base::grepl(paste(adduct_Mion, collapse = "|"), ano$adduct[drop = T]),]
        iso_adduct <- base::gsub("\\D", "", adduct_temp$isotopes[drop = T])
        iso_adduct <- iso_adduct[iso_adduct != ""]
        iso_group <- base::unique(c(iso_group,iso_adduct))
        #leave only iso_group features
        rel_temp <- ano[base::grepl(base::paste(iso_group, collapse = "|"), ano$isotopes[drop = T])|
                          base::grepl(base::paste(adduct_Mion, collapse = "|"), ano$adduct[drop = T]),]
      }
    }
    
    if(base::exists("rel_temp")) {ano <- ano[ano$mz %in% rel_temp$mz,]}
  }
  
  #filter ano by selecting only annotated Features
  if (onlyAnnotated) {
    ano <- dplyr::filter(ano, isotopes != "" | ano$adduct != "")  
  }
  
  #add text for labels to ano
  ano$text <- paste0(round(ano$mz, digits = 2),ano$isotopes, ano$adduct)
  
  #check the col names for duplicates
  sampleNamesRepeats <- ifelse(xA[[1]]@xcmsSet$sample_name %in% xA[[1]]@xcmsSet$sample_group,
                               paste0(xA[[1]]@xcmsSet$sample_name,".1"),xA[[1]]@xcmsSet$sample_name)
  
  #calculate averages for the intensites
  anoS <- dplyr::select(ano, -all_of(xA[[1]]@xcmsSet$sample_name), -all_of(xA[[1]]@xcmsSet$sample_group), -all_of(sampleNamesRepeats))
  
  for (r in rIndex) {
    anoS[,rGroups[r]] <- apply(ano[,colnames(ano) %in% sampleNamesRepeats[xA[[1]]@xcmsSet$sample_group %in% rGroups[r]], drop = F],
                               MARGIN = 1, FUN = function(x) base::ifelse(log, base::log(base::mean(x)), base::mean(x))) 
  }
  
  if (!base::is.null(mzWindowPlot)) {
    if (base::length(mzWindowPlot) == 2) mzrange <- mzWindowPlot
  } else {
    mzrange <- c(base::min(anoS$mz)*0.9, base::max(anoS$mz)*1.1)
  }
  intrange <- c(0, base::max(anoS[,colnames(anoS) %in% rGroups], na.rm = TRUE)*2) 
  
  colors <- ntsIUTA::getColors(base::length(rGroups))
  
  fig <- plotly::plot_ly(anoS, type = "bar")
  
  for (s in rIndex) { 
    
    fig <- fig %>% plotly::add_bars(x = anoS$mz[anoS$rGroup %in% rGroups[s]],
                                    y = anoS[anoS$rGroup %in% rGroups[s],colnames(anoS) %in% rGroups[s], drop = T],
                                    marker = list(color = colors[s]),  width = 0.05,
                                    text = anoS$text[anoS$rGroup %in% rGroups[s]],
                                    textposition = 'outside', textangle = -90, insidetextanchor = "middle",
                                    textfont = list(size = 10, color = colors[s]),
                                    name = rGroups[s])
  }
  
  #title <- base::list(text = "Coisas", x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))
  
  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "m/z",
                      titlefont = base::list(size = 12, color = "black"),
                      range = mzrange)
  
  yaxis = list(linecolor = plotly::toRGB("black"), linewidth = 2, title = base::ifelse(log, "log(Intensity)", "Intensity"),
               titlefont = list(size = 12, color = "black"), range = intrange) #, range = intrange
  
  fig <- fig %>% plotly::layout(xaxis = xaxis, yaxis = yaxis, barmode = "overlay", uniformtext = list(minsize = 4, mode = "show"))
  
  return(fig)
  
}




#' @title plotFeatures_Old
#' @description Plot features from a \linkS4class{featureGroups}.
#'
#' @param features A \linkS4class{featureGroups} object with one or more files and grouped peaks (i.e., features).
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param ID The identifier of the features of interest. When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.
#' @param mz Optional target \emph{m/z} to find features using the mass deviation specified by the \code{ppm} parameter.
#' @param ppm The mass deviation to extract the features when \code{mz} is specified.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below. Only used when \code{mz} is specified.
#' @param rtWindow The time deviation to collect features. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{min}.
#' @param title The title for the plot, optional.
#' @param plotBy The grouping for the plot. Possible are "samples", "replicates" and "features"
#' for plotting by individual samples, replicate samples or features, respectively. The default is "features".
#'
#' @return A plot produced through pkg{gglpot2}.
#'
#' @note If \code{ID} and \code{mz} are \code{NULL} all the features in the \linkS4class{featureGroups} object are plotted.
#'
#' @export
#'
#' @importClassesFrom patRoon featureGroups
#' @importClassesFrom xcms XCMSnExp
#' @importMethodsFrom patRoon algorithm analyses
#' @importMethodsFrom xcms featureDefinitions chromPeaks
#' @importFrom xcms featureChromatograms
#' @importFrom BiocGenerics as.data.frame
#' @importFrom plotly toRGB plot_ly add_trace layout
#'
#' @examples
#'
plotFeatures_Old <- function(features = features,
                            fileIndex = NULL,
                            ID = NULL,
                            mz = NULL, ppm = 5,
                            rt = NULL, rtWindow = 1,
                            rtUnit = "min",
                            title = NULL, plotBy = "features") {
  
  fileIndex <- NULL
  ID <- NULL #"M748_R891_1274"
  mz <- 748.4842
  rt <- 14.9
  rtUnit <- "min"
  ppm <- NULL
  rtWindow <- NULL
  plotBy <- "features"
  title <- NULL
  
  x <- featuresOpenms
  
  if (!is.null(fileIndex)) x <- x[fileIndex, ]
  
  if (!is.null(ID)) {
    x <- x[, ID]
  } else {
    if (!is.null(mz)) {
      mzr <- mzrBuilder(mz = mz, ppm = ppm)
      rtr <- rtrBuilder(rt = rt, rtWindow = rtWindow, rtUnit = rtUnit)
      x <- patRoon::filter(obj = x, retentionRange = rtr, mzRange = mzr)
    }
  }
  
  ft <- names(x)
  
  sInfo <- analyses(x)
  
  pks <- patRoon::as.data.frame(x,
                                average = FALSE,
                                areas = FALSE,
                                features = TRUE,
                                regression = FALSE
  )
  
  EICs <- list()
  
  for (i in seq_len(nrow(pks))) {
    flidx <- which(sInfo == pks$analysis[i])
    EICs[[rownames(pks)[i]]] <- extractEIC(x, fileIndex = flidx,
                                           mz = c(pks$mzmin[i], pks$mzmax[i]),
                                           rtWindow = c(pks$retmin[i], pks$retmax[i]),
                                           rtUnit = "sec")
  }
  
  
  
  
  
  
  if (algorithm(x) == "openms") {
    x <- getXCMSnExp(x)
    y <- temp
    sn <- matrix(rep(0, nrow(chromPeaks(y))), ncol = 1)
    colnames(sn) <- "sn"
    chromPeaks(y) <- cbind(chromPeaks(y), sn)
    featureDefinitions(y) <- featureDefinitions(temp)
  }
  
  
  
  def <- xcms::featureDefinitions(y, type = "within", msLevel = 1)
  
  chrom <- xcms::featureChromatograms(y, aggregationFun = "sum", expandRt = 60, include = "feature_only", filled = TRUE, missing = 0)
  
  head(xcms::chromPeaks(y, isFilledColumn = TRUE))$sn <- 0
  
  sn <- matrix(rep(0, nrow(xcms::chromPeaks(y))), ncol = 1)
  colnames(sn) <- "sn"
  
  head(sn)
  
  xcms::chromPeaks(y) <- cbind(xcms::chromPeaks(y), sn)
  
  xcms::featureDefinitions(y) <- xcms::featureDefinitions(x@xdata)
  
  
  
  features <- base::character()
  for (f in 1:base::length(FT)) {
    features <- c(features, paste("M", base::round(defFT$mzmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                  "_R", base::round(defFT$rtmed[base::row.names(defFT) %in% FT[f],drop = T], digits = 0),
                                  "_", FT[f],sep=""))
  }
  
  # temp_rtFT <- as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
  # chrom <- MSnbase::chromatogram(y, aggregationFun = "sum", rt = rtr, msLevel = 1, missing = 0,
  #                                mz = c(base::min(temp_rtFT$mzmin-(ppm/1E6*temp_rtFT$mzmin), na.rm = T),
  #                                       base::max(temp_rtFT$mzmax+(ppm/1E6*temp_rtFT$mzmax), na.rm = T)))
  chrom <- xcms::featureChromatograms(y, features = FT, aggregationFun = "sum", expandRt = 60, include = "any", filled = TRUE, missing = 0)
  
  df <- base::data.frame(rtime = base::as.numeric(), intensity = base::as.numeric(), sample = base::as.character(), FT = base::as.character())
  for (f in 1:base::nrow(def)) {
    for (s in 1:base::ncol(chrom)) {
      temp <- chrom[f, s]
      temp <- BiocGenerics::as.data.frame(temp)
      temp$sample <- y$sample_name[s]
      temp$ft <- ft[f]
      df <- base::rbind(df,temp)
    }
  }
  
  #colors <- ntsIUTA::getColors(y, "samples")
  
  if (base::is.null(title)) title <- if (plotBy == "features") base::ifelse(base::length(FT) == 1, features, "")
  
  title <- base::list(text = title, x = 0.1, y = 0.98, font = base::list(size = 14, color = "black"))
  
  xaxis <- base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Retention Time (sec.)",
                      titlefont = base::list(size = 12, color = "black"),
                      range = c(rtr[1], rtr[2]), autotick = T, ticks = "outside")
  
  yaxis = base::list(linecolor = plotly::toRGB("black"), linewidth = 2, title = "Intensity",
                     titlefont = list(size = 12, color = "black"))
  
  
  if (plotBy == "samples") {
    colors <- ntsIUTA::getColors(y, "samples")
    legG <- y$sample_name
  } else {
    if (plotBy == "features") {
      colors <- ntsIUTA::getColors(base::length(features))
      legG <- features
    } else {
      colors = ntsIUTA::getColors(y, "groups")
      legG <- y$sample_name
    }
  }
  
  
  # if (plotBy == "replicates") legG <- y$sample_group
  # if (plotBy == "samples") legG <- y$sample_name
  
  #By feature the eic followed by the integrated are, legend grouped as function input plotBy
  plot <- plotly::plot_ly(df)
  showlegend = base::rep(0,base::length(features))
  for (s in 1:base::ncol(chrom)) { #base::ncol(feat)
    for (f in 1:base::length(FT)) { #base::nrow(feat)
      showlegend[f] = 1 + showlegend[f]
      plot <- plot %>% plotly::add_trace(df,
                                         x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "rtime"],
                                         y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s], "intensity"],
                                         type = "scatter", mode = "lines",
                                         line = base::list(width = 0.5, color = base::unname(colors[base::ifelse(plotBy == "features", f, s)])),
                                         connectgaps = TRUE,
                                         name = legG[base::ifelse(plotBy == "features", f, s)],
                                         legendgroup = legG[base::ifelse(plotBy == "features", f, s)],
                                         showlegend = base::ifelse(plotBy == "samples", base::ifelse(f == 1, T, F),
                                                                   base::ifelse(plotBy == "features", base::ifelse(showlegend[f] == 1, T, F),
                                                                                base::ifelse(f == 1, T, F))))
      
      rtFT <- base::as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
      if(unique(TRUE %in% (rtFT$sample == s))) {
        rtFT <- rtFT[rtFT$sample == s, c("rtmin","rtmax")]
        plot <- plot %>%  plotly::add_trace(df,
                                            x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2], "rtime"],
                                            y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2],"intensity"],
                                            type = "scatter", mode =  "lines+markers", fill = 'tozeroy', connectgaps = TRUE,
                                            fillcolor = base::paste(color = base::unname(colors[base::ifelse(plotBy == "features", f, s)]),50, sep = ""),
                                            line = base::list(width = 0.1, color = base::unname(colors[base::ifelse(plotBy == "features", f, s)])),
                                            marker = base::list(size = 3, color = base::unname(colors[base::ifelse(plotBy == "features", f, s)])),
                                            name = legG[base::ifelse(plotBy == "features", f, s)],
                                            legendgroup = legG[base::ifelse(plotBy == "features", f, s)],
                                            showlegend = F,
                                            hoverinfo = 'text', text = base::paste('</br> feature: ', features[f],
                                                                                   '</br> sample: ', y$sample_name[s]))
      }
    }
  }
  # showlegend = 0
  # for (s in 1:base::ncol(chrom)) { #base::ncol(feat)
  #   for (f in 1:base::nrow(chrom)) { #base::nrow(feat)
  #     rtFT <- base::as.data.frame(xcms::chromPeaks(y)[base::unlist(defFT[base::row.names(defFT) %in% FT[f] ,"peakidx"]), ])
  #     if(base::unique(TRUE %in% (rtFT$sample == s))) {
  #       showlegend <- 1 + showlegend
  #       rtFT <- rtFT[rtFT$sample == s, c("rtmin","rtmax")]
  #       plot <- plot %>%  plotly::add_trace(df,
  #                                           x = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2], "rtime"],
  #                                           y = df[df$FT == FT[f] & df$sample ==  y$sample_name[s] & df$rtime >= rtFT[1,1] & df$rtime <= rtFT[1,2],"intensity"],
  #                                           type = "scatter", mode = "lines+markers", fill = 'tozeroy', connectgaps = TRUE, fillcolor = paste(color = base::unname(colors[s]),50, sep = ""),
  #                                           line = base::list(width = 0.1, color = base::unname(colors[s])),
  #                                           marker = base::list(size = 3, color = base::unname(colors[s])),
  #                                           name = features[f], legendgroup = FT[f], showlegend = base::ifelse(showlegend == 1, T, F))
  #     }
  #   }
  # }
  
  plot <- plot %>% plotly::layout(xaxis = xaxis,yaxis = yaxis, title = title) #legend = list(title = list(text='<b> Sample: </b>'))
  
  return(plot)
  
}




#' @title plotlyRawChrom
#' @description Plots total, base and extracted ion chromatograms
#' (TIC, BPC and EIC, respectively) of an \linkS4class{ntsData} object,
#' using the \pkg{plotly} package to produce an iterative plot.
#'
#' @param obj An \linkS4class{ntsData} object with one or more files.
#' @param fileIndex The index of the file/s to extract the data.
#' @param mz Optional target \emph{m/z} to obtain an EIC.
#' Note that when not \code{NULL} (the default), EIC is always returned.
#' @param ppm The mass deviation to extract the data for the EIC, in \code{ppm}.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time window or deviation to collect the data.
#' The time unit is defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' A vector of length 1 is assumed as a deviation of a given \code{rt}.
#' @param rtUnit Possible entries are \code{sec} (the default) or \code{min}.
#' @param msLevel The MS level to extract the data.
#' For the moment, only 1 is possible.
#' @param type The type of chromatogram.
#' Possible entries are "bpc" for base peak chromatogram
#' or "tic" for total ion chromatogram.
#' The default is "tic". If \code{mz} is specified (not \code{NULL}),
#' the type is set automatically to EIC.
#' @param colorBy Possible values are \code{"samples"} or \code{samplegroups}
#' (the default), for colouring by samples or sample replicate groups respectively.
#'
#' @return An iterative plot for inspection of the raw data
#' in the given \linkS4class{ntsData} object.
#'
#' @export
#'
#' @importMethodsFrom MSnbase filterFile
#' @importFrom plotly toRGB plot_ly add_trace layout
#' @importFrom dplyr group_by arrange top_n summarize
#'
#' @examples
#'
plotlyRawChrom <- function(obj = NULL, fileIndex = NULL,
                           mz = NULL, ppm = 20,
                           rt = NULL, rtWindow = NULL,
                           rtUnit = "sec",
                           msLevel = 1,
                           type = "tic",
                           colorBy = "samplegroups") {
  
  # raw = rawdata
  # fileIndex = NULL
  # mz <- c(233.0243)
  # rt <- NULL
  # rtUnit = "min"
  # ppm <- 20
  # rtWindow = NULL
  # msLevel = 1
  
  if (!is.null(mz) && length(mz) == 1) {
    if (!is.null(ppm)) ppm <- 20
    main <- paste0("EIC of ", round(mz, digits = 4), " +/- ", round(ppm, digits = 0), " ppm")
    type <- "eic"
  } else {
    if (!is.null(mz) && length(mz) == 2) {
      main <- paste0("EIC for ", round(mz[1], digits = 4), "to ", round(mz[2], digits = 4))
      type <- "eic"
    } else {
      main <- toupper(type)
    }
  }
  
  if (!is.null(fileIndex)) obj <- obj[fileIndex]
  
  df <- extractEIC(obj = obj,
                   fileIndex = NULL,
                   mz = mz, ppm = ppm,
                   rt = rt, rtWindow = rtWindow,
                   rtUnit = rtUnit, msLevel = 1,
                   normIntensity = FALSE)
  
  if (type == "tic") {
    mz <- df %>% group_by(rt) %>% top_n(1, i)
    mz <- as.data.frame(mz)
    mz <- arrange(mz, rt)
    y <- df %>% group_by(rt) %>% summarize(i = sum(i))
    y <- as.data.frame(y)
    y <- arrange(y, rt)
    mz$i <- y[, "i", drop = TRUE]
    df <- mz
  }
  
  if (type == "bpc") {
    y <- df %>% group_by(rt) %>% top_n(1, i)
    y <- as.data.frame(y)
    df <- arrange(y, rt)
  }
  
  for (i in seq_len(length(unique(df$file)))) {
    df[df$file == i, "file"] <- samples(obj)[i]
  }
  
  cl <- getColors(obj, which = colorBy)
  
  title <- list(text = main, x = 0.1, y = 0.98, font = list(size = 14, color = "black"))
  
  xaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Retention Time (sec.)",
                titlefont = list(size = 12, color = "black"))
  
  yaxis <- list(linecolor = toRGB("black"),
                linewidth = 2, title = "Intensity",
                titlefont = list(size = 12, color = "black"))
  
  plot <- plot_ly(df,
                  x = df[df$file == samples(obj)[1], "rt"],
                  y = df[df$file == samples(obj)[1], "i"],
                  type = "scatter", mode = "lines+markers",
                  line = list(width = 0.5, color = unname(cl[1])),
                  marker = list(size = 2, color = unname(cl[1])),
                  name = samples(obj)[1])
  
  if (length(unique(df$file)) > 1) {
    for (i in 2:length(unique(df$file))) {
      plot  <- plot %>% add_trace(df,
                                  x = df[df$file == samples(obj)[i], "rt"],
                                  y = df[df$file == samples(obj)[i], "i"],
                                  type = "scatter", mode = "lines+markers",
                                  line = list(width = 0.5, color = unname(cl[i])),
                                  marker = list(size = 2, color = unname(cl[i])),
                                  name = samples(obj)[i])
    }
  }
  
  plot <- plot %>% layout(legend = list(title = list(text = "<b> Sample: </b>")),
                          xaxis = xaxis, yaxis = yaxis, title = title)
  
  return(plot)
  
}



#' @title plotRawChromOld
#' @description Plots total, base and extracted ion chromatograms
#' (TIC, BPC and EIC, respectively) of an \linkS4class{ntsData} object.
#'
#' @param obj An \linkS4class{ntsData} object with one or more files.
#' @param fileIndex The index of the file/s to extract the data.
#' @param mz Optional target \emph{m/z} to obtain an EIC.
#' Note that when not \code{NULL} (the default), EIC is always returned.
#' @param ppm The mass deviation to extract the data for the EIC, in \code{ppm}.
#' @param rt The retention time in minutes or seconds,
#' depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time window or deviation to collect the data.
#' The time unit is defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector,
#' defining the minimum and maximum retention time.
#' A vector of length 1 is assumed as a deviation of a given \code{rt}.
#' @param rtUnit Possible entries are \code{sec} (the default) or \code{min}.
#' @param msLevel The MS level to extract the data.
#' For the moment, only 1 is possible.
#' @param type The type of chromatogram.
#' Possible entries are "bpc" for base peak chromatogram
#' or "tic" for total ion chromatogram.
#' The default is "tic". If \code{mz} is specified (not \code{NULL}),
#' the type is set automatically to EIC.
#' @param colorBy Possible values are \code{"samples"} or \code{samplegroups}
#' (the default), for colouring by samples or sample replicate groups respectively.
#'
#' @return A plot for inspection of the raw data
#' in the given \linkS4class{ntsData} object.
#'
#' @export
#'
#' @importMethodsFrom ProtGenerics rtime
#' @importMethodsFrom MSnbase filterFile rtime
#' @importFrom plotly toRGB plot_ly add_trace layout
#' @importFrom dplyr group_by arrange top_n summarize
#'
#' @examples
#'
plotRawChromOld <- function(obj = NULL, fileIndex = NULL,
                          mz = NULL, ppm = 20,
                          rt = NULL, rtWindow = NULL,
                          rtUnit = "sec",
                          msLevel = 1,
                          type = "tic",
                          colorBy = "samplegroups") {
  
  # raw = rawdata
  # fileIndex = NULL
  # mz <- c(233.0243)
  # rt <- NULL
  # rtUnit = "min"
  # ppm <- 20
  # rtWindow = NULL
  # msLevel = 1
  
  if (!is.null(mz) && length(mz) == 1) {
    if (!is.null(ppm)) ppm <- 20
    main <- paste0("EIC of ", round(mz, digits = 4), " +/- ", round(ppm, digits = 0), " ppm")
    type <- "eic"
  } else {
    if (!is.null(mz) && length(mz) == 2) {
      main <- paste0("EIC for ", round(mz[1], digits = 4), "to ", round(mz[2], digits = 4))
      type <- "eic"
    } else {
      main <- toupper(type)
    }
  }
  
  if (!is.null(fileIndex)) obj <- obj[fileIndex]
  
  df <- extractEIC(obj = obj,
                   fileIndex = NULL,
                   mz = mz, ppm = ppm,
                   rt = rt, rtWindow = rtWindow,
                   rtUnit = rtUnit, msLevel = 1,
                   normIntensity = FALSE)
  
  if (type == "tic") {
    mz <- df %>% group_by(rt) %>% top_n(1, i)
    mz <- as.data.frame(mz)
    mz <- arrange(mz, rt)
    y <- df %>% group_by(rt) %>% summarize(i = sum(i))
    y <- as.data.frame(y)
    y <- arrange(y, rt)
    mz$i <- y[, "i", drop = TRUE]
    df <- mz
  }
  
  if (type == "bpc") {
    y <- df %>% group_by(rt) %>% top_n(1, i)
    y <- as.data.frame(y)
    df <- arrange(y, rt)
  }
  
  for (i in seq_len(length(unique(df$file)))) {
    df[df$file == i, "file"] <- samples(obj)[i]
  }
  
  cl <- getColors(obj, which = colorBy)
  
  plot <- ggplot(data = df, aes(x = rt, y = i, color = file)) +
    geom_line(size = 0.5) +
    scale_color_manual(values = cl) +
    ggtitle(main) +
    theme_bw() +
    ylab("Intensity") +
    xlab("Retention Time") +
    theme(legend.title = element_blank())
  
  return(plot)
  
}



#' @param ppmForFillingGroups The mass (in ppm) to expand the \emph{mz} for filling missing peaks in incomplete features.
#' See \code{\link[xcms]{fillChromPeaks}} for more information.

#Fill missing peaks in feature groups
  #Old parameter used: FillChromPeaksParam(ppm = ppmForFillingGroups)
  featData <- base::suppressWarnings(xcms::fillChromPeaks(featData, param = xcms::ChromPeakAreaParam(),
                                                          BPPARAM = BiocParallel::bpparam("SnowParam")))
  
  # add optinally to apply correction of retention time, avoiding another function from xcms
  # if (length(peaksDataUnified$sample_name) > 1) {
  #   featData <- xcms::applyAdjustedRtime(featData)
  # }



#' @title refinePeaks
#' @description Refines chromatographic peaks
#'
#' @param peaksData The \linkS4class{OnDiskMSnExp} object generated by the \code{\link{importRawData}} function.
#' @param expandRt Time (in seconds) to expand the peak width to merge neighboring peaks.
#' @param minProp The proportion (between 0 and 1) representing the proporion of intensity to be required for peaks to be joined.
#' @param expandMz Decreases the minimum and increases the maximum \emph{m/z} by a fixed given value.
#' @param ppm Decreases the minimum and increases the maximum \emph{m/z} by a fixed given value, in ppm.
#' @param maxPeakwidth The maximum width allowed for a peak, in seconds.
#' @param save Logical, set to \code{TRUE} to save the generated list of \code{XCMSnExp} objects in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#' @param maxMultiProcess Logical, set to \code{TRUE} to enable max parallel processing. Changes the number of workers to the maximum available.
#'
#' @return A \code{list} of \linkS4class{XCMSnExp} objects with length equal to the number of replicate groups defined in the \code{setup}.
#' The \code{list} is named with the name given to the sample replicate groups.
#'
#' @note See parameters in \code{?refineChromPeaks} from \code{xcms} for more information when \code{refinePeaks} is set to \code{TRUE}.
#'
#' @references
#' \insertRef{MSnbase1}{ntsIUTA}
#' \insertRef{MSnbase2}{ntsIUTA}
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#'
#' @export
#'
#' @importFrom BiocParallel registered register bpparam
#' @importClassesFrom BiocParallel SnowParam
#' @importFrom parallel detectCores
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @importClassesFrom xcms MergeNeighboringPeaksParam CleanPeaksParam
#' @importMethodsFrom xcms findChromPeaks refineChromPeaks
#' @importMethodsFrom MSnbase filterFile
#'
#' @examples
#'
refinePeaks <- function(peaksData,
                        expandRt = 1, minProp = 0.9, expandMz = 0, ppm = 0, maxPeakwidth = NULL,
                        save = TRUE, projPath = setup$projPath, maxMultiProcess = TRUE) {

  mpp <- xcms::MergeNeighboringPeaksParam(expandRt = expandRt, minProp = minProp, expandMz = expandMz, ppm = ppm)
  peaksData[[groups[rgidx]]] <- xcms::refineChromPeaks(peaksData[[groups[rgidx]]],
                                                        mpp, msLevel = 1, BPPARAM = BiocParallel::bpparam("SnowParam"))
  if (!base::is.null(maxPeakwidth)) {
    peaksData[[groups[rgidx]]] <- xcms::refineChromPeaks(peaksData[[groups[rgidx]]],
                                                    param = xcms::CleanPeaksParam(maxPeakwidth = maxPeakwidth),
                                                    msLevel = 1L, BPPARAM = BiocParallel::bpparam("SnowParam"))
  }
}



#' @title makeFeatures
#' @description The \code{makeFeatures} consists of alignment and grouping of chromatographic peaks across samples.
#' The function uses methods from the package \pkg{xcms}, specifically \code{\link[xcms]{adjustRtime}} and \code{\link[xcms]{groupChromPeaks}}
#' methods. The input is a \linkS4class{XCMSnExp} object or a \code{list} of \linkS4class{XCMSnExp} objects
#' named with the given sample replicate group during the \code{setup}. If the input is a \code{list}, \linkS4class{XCMSnExp} objects
#' are concatenated to produce a merged \linkS4class{XCMSnExp} object,
#' containing features (\emph{i.e.} grouped chromatographic peaks across samples).
#'
#' @param peaksData A \linkS4class{XCMSnExp} object or a \code{list} of \linkS4class{XCMSnExp} objects named according to the 
#' sample replicate group represented by each \linkS4class{XCMSnExp} object as obtained by \code{\link{peakPicking}}.
#' @param paramPreGrouping If \code{\link[xcms]{adjustRtime}} is preformed with the method \code{PeakGroups},
#' a pre-grouping of peaks is required for alignment. The \code{paramPreGrouping} is the parameters obtained by the selected grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} for more information.
#' @param paramAlignment The parameters for the chosen alignment method. See documentation of \code{\link[xcms]{adjustRtime}} for more information.
#' @param paramGrouping The parameters for the chosen grouping method.
#' See documentation of \code{\link[xcms]{groupChromPeaks}} for more information.
#' @param ppmForFillingGroups The mass (in ppm) to expand the \emph{mz} for filling missing peaks in incomplete features.
#' See \code{\link[xcms]{fillChromPeaks}} for more information.
#' @param save Logical, set to \code{TRUE} to save the generated \code{XCMSnExp} object in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#' @param maxMultiProcess Logical, set to \code{TRUE} to enable max parallel processing. Changes the number of workers to the maximum available.
#'
#' @return An \linkS4class{XCMSnExp} containing features which are grouped and aligned peaks across samples.
#' 
#'  @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#' 
#' @importFrom BiocParallel registered SnowParam register bpparam
#' @importFrom parallel detectCores
#' @importFrom xcms sampleGroups groupChromPeaks adjustRtime fillChromPeaks ChromPeakAreaParam
#'
#' @examples
#' 
#' 
#' 
makeFeatures <- function(peaksData = peaksData,
                         paramPreGrouping = instParam$preGrouping,
                         paramAlignment = instParam$alignment,
                         paramGrouping = instParam$grouping,
                         ppmForFillingGroups = 5,
                         save = TRUE, projPath = setup$projPath, maxMultiProcess = TRUE) {
  
  #Examples
  #setup <- ntsIUTA::makeSetup(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)
  # setup$sampleInfo[1:3,"group"] <- "Sample"
  # rawDataExample <- ntsIUTA::importRawData(setup$sampleInfo[1:3,], save = FALSE, centroidedData = TRUE)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = param, save = FALSE)
  # paramList <- paramListExample
  # instParam <- getInstParam(paramList)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = instParam$PP, save = FALSE)
  # featDataExample <- ntsIUTA::makeFeatures(peaksDataExample, paramPreGrouping = instParam$preGrouping, paramAlignment = instParam$alignment, paramGrouping = instParam$grouping, save = FALSE)
  # featDataExample
  
  
  
  # require(xcms)
  # require(BiocParallel)
  # require(parallel)
  # require(grDevices)
  # require(RColorBrewer)
  # require(dplyr)
  # require(graphics)
  
  #Enable full parallel processing
  if (maxMultiProcess)
  {
    snow <- BiocParallel::registered("SnowParam")
    if (snow$workers < parallel::detectCores())
    {
      snow <- BiocParallel::SnowParam(workers = parallel::detectCores(), type = "SOCK", exportglobals = FALSE)
      BiocParallel::register(snow, default = TRUE)
    }
  }
  
  
  #Concatenation for peaksData objects before grouping
  if (base::class(peaksData) == "XCMSnExp")
  {
    peaksDataUnified <- peaksData
  } else {
    
    if (base::length(peaksData) == 1)
    {peaksDataUnified <- peaksData[[1]]}
    else {peaksDataUnified <- base::do.call(c, base::unlist(peaksData, recursive = FALSE))}
    
  }
  
  if (base::length(peaksDataUnified$sample_name) > 1) {
    
    #Pregrouping necessary for alignment with PeakGroups method from xcms
    if(base::class(paramAlignment) == "PeakGroupsParam")
    {
      xcms::sampleGroups(paramPreGrouping) <- peaksDataUnified$sample_group
      featData <- xcms::groupChromPeaks(peaksDataUnified, param = paramPreGrouping)
    }
    
    
    featData <- xcms::adjustRtime(featData, msLevel = 1, param = paramAlignment)
    
    xcms::sampleGroups(paramGrouping) <- peaksDataUnified$sample_group
    featData <- xcms::groupChromPeaks(featData, param = paramGrouping)
    
  } else {
    
    xcms::sampleGroups(paramGrouping) <- peaksDataUnified$sample_group
    featData <- xcms::groupChromPeaks(peaksDataUnified, param = paramGrouping)
    
  }
  
  #Fill missing peaks in feature groups
  #Old parameter used: FillChromPeaksParam(ppm = ppmForFillingGroups)
  featData <- base::suppressWarnings(xcms::fillChromPeaks(featData, param = xcms::ChromPeakAreaParam(),
                                                          BPPARAM = BiocParallel::bpparam("SnowParam")))
  
  # add optinally to apply correction of retention time, avoiding another function from xcms
  # if (length(peaksDataUnified$sample_name) > 1) {
  #   featData <- xcms::applyAdjustedRtime(featData)
  # }
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(featData, file = base::paste0(rData,"\\featData.rds"))
  }
  
  return(featData)
}







#' @title getInstParam
#' @description Tranforms the list of parameters for peak picking (PP), alignment (\emph{i.e.} retention time adjustment across samples)
#' and grouping into the respetive class objects.
#' As input a list as follows should be given where for each data processing step the input should be given as described in
#' the documentation of the aimed function of the \pkg{xcms} package.
#' \itemize{
#'   \item \code{paramList <- list(PP = list("..."), preGrouping = list("..."), alignment = list("..."), grouping = list("..."))}
#' }
#' \code{preGrouping} is used for \code{\link[xcms]{adjustRtime}} function
#' when using the method \code{PeakGroups} as it requires peaks groups (\emph{i.e.} features) for alignment.
#' Note that currently for PP only \code{MassifquantParam} is possible, for alignment only \code{PeakGroupsParam} is possible
#' and for grouping only \code{PeakDensityParam}. Other methods will be integrated in the future.
#' 
#' @param paramList A list of parameters for PP, alignment and grouping as described above.
#' See \code{\link{paramListExample}} object for a structure example.
#'
#' @return A \code{list} with parameters for PP, aligment and grouping to be passed to \code{\link{peakPicking}} and \code{\link{makeFeatures}}. 
#' 
#'  @references
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#'
#' @examples
#' paramList <- ntsIUTA::paramListExample
#' getInstParam(paramList)
#' 
getInstParam <- function(paramList = paramList) {
  
  instParam <- list()
  
  instParam$PP <- NULL
  
  # Make instParam for PP
  if(!is.null(paramList$PP)) {
    
    if(paramList$PP$funcName == class(xcms::MassifquantParam())) {
      
      instParam$PP <- xcms::MassifquantParam(ppm = paramList$PP$ppm,
                                             peakwidth = paramList$PP$peakwidth,
                                             snthresh = paramList$PP$snthresh,
                                             prefilter = paramList$PP$prefilter,
                                             mzCenterFun = paramList$PP$mzCenterFun,
                                             integrate = paramList$PP$integrate,
                                             mzdiff = paramList$PP$mzdiff,
                                             fitgauss = paramList$PP$fitgauss,
                                             noise = paramList$PP$noise,
                                             criticalValue = paramList$PP$criticalValue,
                                             consecMissedLimit = paramList$PP$consecMissedLimit,
                                             unions = paramList$PP$unions,
                                             checkBack = paramList$PP$checkBack,
                                             withWave = paramList$PP$verboseColumns,
                                             verboseColumns = paramList$PP$verboseColumns)
      
    }
    
    #Add parameters for other methods
    
  }
  
  instParam$preGrouping <- NULL
  instParam$alignment <- NULL
  
  # Make instParam for alignment, including preGrouping for PeakGroupsParam method
  if(!is.null(paramList$alignment)) {
    
    if(paramList$alignment$funcName == class(xcms::PeakGroupsParam())) {
      
      instParam$preGrouping <-  xcms::PeakDensityParam(sampleGroups = "placeHolder",
                                                       bw = paramList$preGrouping$bw,
                                                       minFraction = paramList$preGrouping$minFraction,
                                                       minSamples = paramList$preGrouping$minSamples,
                                                       binSize = paramList$preGrouping$binSize,
                                                       maxFeatures = paramList$preGrouping$maxFeatures)
      
      instParam$alignment <- xcms::PeakGroupsParam(minFraction = paramList$alignment$minFraction,
                                                   extraPeaks = paramList$alignment$extraPeaks,
                                                   smooth = paramList$alignment$smooth,
                                                   span = paramList$alignment$span,
                                                   family = paramList$alignment$family,
                                                   peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
                                                   subset = integer(),
                                                   subsetAdjust = "average")
      
    }
    
    #Add parameters for other methods
    
  }
  
  instParam$grouping <- NULL
  
  #Make instParam for grouping.
  if(!is.null(paramList$grouping)) {
    
    if(paramList$grouping$funcName == class(xcms::PeakDensityParam(sampleGroups = "placeHolder"))) {
      
      instParam$grouping <- PeakDensityParam(sampleGroups = "placeHolder",
                                             bw = paramList$grouping$bw,
                                             minFraction = paramList$grouping$minFraction,
                                             minSamples = paramList$grouping$minSamples,
                                             binSize = paramList$grouping$binSize,
                                             maxFeatures = paramList$grouping$maxFeatures)
      
    }
    
    #Add parameters for other methods
    
  }
  
  return(instParam)
  
}



### paramList -----

# paramListExample <- list(
#   
#   instName = "Generic",
#   
#   PP = list(
#     funcName = "MassifquantParam",
#     ppm = 20,
#     peakwidth = c(8, 60),
#     snthresh = 3,
#     prefilter = c(6, 500),
#     mzCenterFun = "wMean",
#     integrate = 1,
#     mzdiff = -0.001,
#     fitgauss = TRUE,
#     noise = 200,
#     criticalValue = 1.5,
#     consecMissedLimit = 2,
#     unions = 1,
#     checkBack = 1,
#     withWave = TRUE,
#     verboseColumns = TRUE
#   ),
#   
#   preGrouping = list(
#     funcName = "PeakDensityParam",
#     bw = 5,
#     minFraction = 0.5,
#     minSamples = 1,
#     binSize = 0.003,
#     maxFeatures = 100
#   ),
#   alignment = list(
#     funcName = "PeakGroupsParam",
#     minFraction = 1,
#     extraPeaks = 0,
#     smooth = "loess",
#     span = 0.2,
#     family = "gaussian"
#   ),
#   grouping = list(
#     funcName = "PeakDensityParam",
#     bw = 3,
#     minFraction = 0.5,
#     minSamples = 1,
#     binSize = 0.003,
#     maxFeatures = 100
#   ),
#   componentization = list(
#     
#   )
# )


#additional function for ploting adepted XIC from MSnbase 
.vertical_sub_layout_ex <- function(x, sub_plot = 2) {
  sqrt_x <- base::sqrt(base::length(x))
  ncol <- base::ceiling(sqrt_x)
  nrow <- base::round(sqrt_x)
  rws <- base::split(1:(ncol * nrow * sub_plot), f = base::rep(1:nrow, each = sub_plot * ncol))
  base::do.call(rbind, base::lapply(rws, matrix, ncol = ncol))
}



#additional function for ploting adepted XIC from MSnbase
.plotXIC_ex <- function(x, main = "", col = "white", colramp = colorRampPalette(c("#383E47", "#5E8CAA", "#16B9E5", "#16E5C9", "#16E54C")), #topo.colors colorRampPalette(rev(RColorBrewer::brewer.pal(7,"RdGy")))
                        grid.color = "lightgrey", pch = 21,
                        layout = base::matrix(1:2, ncol = 1), plotTargetMark = plotTargetMark,...) {
  # start edit
  plot_strip <- function(..., v, h) base::plot(...)
  dots <- base::list(...)
  #print(list(...))
  # end edit
  
  if (!plotTargetMark)
  {
    dots$h <- NULL
    dots$v <- NULL
  }
  
  targetColor <- "#0644E9"
  
  if (base::is.matrix(layout))
    graphics::layout(layout)
  
  ## Chromatogram.
  bpi <- base::unlist(base::lapply(base::split(x$i, x$rt), max, na.rm = TRUE))
  brks <- lattice::do.breaks(base::range(x$i), nint = 256)
  graphics::par(mar = c(0, 4, 2, 1))
  
  # start edit
  plot_strip(base::as.numeric(base::names(bpi)), bpi, xaxt = "n", col = col, main = main,
             bg = lattice::level.colors(bpi, at = brks, col.regions = colramp), xlab = "",
             pch = pch, cex = 1.5, ylab = "", las = 2, ...)
  # end edit
  
  if(!is.null(dots$v)) graphics::abline(v=dots$v,col=targetColor,lty=3)
  
  graphics::mtext(side = 4, line = 0, "Intensity", cex = graphics::par("cex.lab"))
  graphics::grid(col = grid.color)
  graphics::par(mar = c(3.5, 4, 0, 1))
  
  # start edit
  plot_strip(x$rt, x$mz, main = "", pch = pch, col = col, xlab = "", ylab = "", cex = 1.5,
             yaxt = "n", bg = lattice::level.colors(x$i, at = brks, col.regions = colramp),
             ...)
  
  if(!is.null(dots$h)) graphics::abline(h=dots$h,col=targetColor,lty=3)
  if(!is.null(dots$v)) graphics::abline(v=dots$v,col=targetColor,lty=3)
  
  if(!is.null(dots$v) & !is.null(dots$h))
  {
    graphics::rect(dots$v-10, dots$h-((5/1E6)*dots$h), dots$v+10, dots$h+((5/1E6)*dots$h),
                   col = NA, lty = 2, border = targetColor)
  }
  # end edit
  
  graphics::axis(side = 2, las = 2)
  graphics::grid(col = grid.color)
  graphics::mtext(side = 1, line = 2.5, "Retention time (sec.)", cex = par("cex.lab"))
  graphics::mtext(side = 4, line = 0, "m/z", cex = par("cex.lab"))
  
}


#' @title plotRawChrom_old
#' @description Plot BPC or TIC chromatograms of an \linkS4class{OnDiskMSnExp} object.
#'
#' @param x An \linkS4class{OnDiskMSnExp} object with one or more files.
#' @param fileIndex The index of the file/s to extract the centroids or profile data.
#' @param mz Optional target \emph{m/z} to obtain an extracted ion chromatogram (EIC).
#' @param ppm The mass deviation to extract the data for the EIC in \code{ppm}.
#' @param rt The retention time in minutes or seconds, depending on the defined \code{rtUnit}, see below.
#' @param rtWindow The time deviation to collect centroids or profile data. The time unit is the defined by \code{rtUnit}.
#' A time interval can be given with a length 2 vector, defining the minimum and maximum retention time.
#' @param rtUnit Possible entries are \code{min} or \code{sec}. The default is \code{sec}.
#' @param msLevel The MS level to extract the data. POssible values are 1 or 2.
#' @param type The type of chromatogram. Possible entries are "bpc" of base peak chromatogram or "tic" for total ion chromatogram.
#' The default is "tic".
#'
#' @return A plot.
#'
#' @examples
#' 
#' plotRawChrom_old(ntsIUTA::rawDataExample)
#' 
#' 
plotRawChrom_old <- function(x = rawData, fileIndex = NULL,
                             mz = NULL, ppm = 20,
                             rt = NULL, rtWindow = NULL,
                             rtUnit = "min",
                             msLevel = 1, type = "tic") {
  
  require(magrittr)
  
  if (rtUnit == "min") if (!is.null(rt)) rt <- rt*60
  if (rtUnit == "min") if (!is.null(rtWindow)) rtWindow <- rtWindow*60
  
  if (!base::is.null(fileIndex)) {x <- MSnbase::filterFile(x, fileIndex)}
  
  rtr <- c(base::min(MSnbase::rtime(x)), base::max(MSnbase::rtime(x)))
  mzr <- NULL
  if (!is.null(mz)) {
    
    if (length(mz) == 1) { mzr <- c(mz - ((ppm/1E6)*mz), mz + ((ppm/1E6)*mz)) }
    if (length(mz) == 2) { mzr <- c(mz[1], mz[2]) }
    
  }
  
  if (!is.null(rt)) { rtr <- c((rt) - ifelse(!is.null(rtWindow), rtWindow, 1), (rt) + ifelse(!is.null(rtWindow), rtWindow, 1)) }
  if (is.null(rt)) if (unique(!is.null(rtWindow))) if (length(rtWindow) == 2) { rtr <- c(rtWindow[1], rtWindow[2]) }
  
  colors <- getColors(x, "samples")
  
  if (is.null(mzr)) {
    chrom <- MSnbase::chromatogram(x, aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))), rt = rtr, msLevel = msLevel)
  } else {
    chrom <- MSnbase::chromatogram(x, aggregationFun = base::ifelse(type == "bpc", "max", base::ifelse(type == "tic", "sum", stop("Unkown type, please use bpc or tic."))), rt = rtr, mz = mzr, msLevel = msLevel)
  }
  
  if (!base::is.null(mzr))
  {
    main <- base::paste0("EIC of ", round(mz, digits = 4)," +/- ", round(ppm, digits = 0)," ppm")
  } else {
    main <- base::toupper(type)
  }
  
  # graphics::layout(base::matrix(1:1))
  # graphics::par(mar=c(5,4,4,2)+0.1)
  
  plot <- MSnbase::plot(chrom, col = colors, lwd = 0.5, main = main, xlab = "Retention time (sec.)", ylab = "Intensity",)
  
  plot <- plot %>% graphics::legend("topleft" , legend = x$sample_name, col = colors, lty = 1, lwd = 2, cex = 0.6, bty = "n", xjust = 0)
  
  
  return(plot)
  
}




#' @title makeSetup
#' @description Project setup by screening a given folder for MS files or
#' by creating a new folder where the MS files can be added using \code{addFiles}.
#' Setting \code{createRproject} to \code{TRUE} will create and open an R project in the selected or given folder.
#' When \code{convertFiles} is \code{TRUE} the function uses \code{\link{mzMLconverter}} to automatically convert specified raw MS files to mzML.
#' 
#'
#' @param projPath The directory for the project. The default is the \code{utils::choose.dir()} to select or create a folder.
#' Note that the function will look into subfolders for mzML files or specified raw files when \code{convertFiles} is set to \code{TRUE}.
#' @param date The project date. The default is the system date generated by \code{\link[base]{Sys.time}}.
#' @param groups A vector with the identifier/name of each sample replicate group. If \code{NULL},
#' a grouping tentative will be made using the name of the MS files found in the given \code{projPath}.
#' @param blanks The name of the sample replicate group to be used for blank subtraction.
#' Different blank groups can be assigned to replicate sample groups by editing the \code{sampleInfo} in the output \code{setup} list.
#' @param polarity A vector specifying the polarity mode used in each MS file.
#' Possible values are \code{positive} or \code{negative} for positive and negative mode, respectively.
#' @param convertFiles Logical, set to \code{TRUE} for non mzML MS files being converted to centroided mzML.
#' The default is \code{FALSE} for no conversion even if files are present.
#' @param convertFrom The name of the file format or vendor format to be found and converted.
#' Possible entries are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.
#' @param convertToCentroid Logical, set to \code{TRUE} to convert profile data to centroided. The default and recommended is \code{TRUE}.
#' @param save Logical, set to \code{TRUE} to save the object setup in the \strong{rData} folder.
#' If \strong{rData} folder does not exist in given \code{projPath}, the folder will be created.
#' @param makeNewProject Logical, set to TRUE to create an R project in the given \code{projPath} and open a new R session.
#'
#' @return The output is a simple list, containing: (1) \code{projPath}, the chosen project path; (2) the \code{date}, the given project date; and 
#' (3) the \code{sampleInfo}, a \code{\link[base]{data.frame}} with the following columns:
#' (1) \code{filePath}, the path of each mzML file, including any converted MS files;
#' (2) \code{sample}, the retrieved file name to be used as sample identifier;
#' (3) \code{group}, the name of each replicate sample group
#' (e.g. samples with the same name but with different numbering, such as in the case of replicate samples);
#' (4) \code{blank}, the specified blank group names or the blank group detected by file name (\emph{i.e.} blank, Blank and/or B as file name);
#' and (5) \code{polarity}, the polarity mode used for each sample, possible entries are \code{positive} and \code{negative};
#' 
#' 
#' @details The data.frame \code{sampleInfo} in the output list can always be edited for correction of
#' sample replicate group names, blank replicate groups and assigned polarity, which is set to \code{positive} by default when not given as argument.
#' 
#' @export
#'
#' @importFrom tools file_path_sans_ext
#' @importFrom dplyr mutate
#' @importFrom rstudioapi initializeProject openProject
#'
#' @examples
#' 
#' 
#' 
setupProject <- function(projPath = utils::choose.dir(base::getwd(), "Select or create a project folder"),
                         date = base::Sys.Date(),
                         groups = NULL,
                         blanks = NULL,
                         polarity = "positive",
                         convertFiles = FALSE,
                         convertFrom = NULL,
                         convertToCentroid = TRUE,
                         save = TRUE,
                         makeNewProject = FALSE) {
  
  #Examples
  # projPath <-  system.file(package = "ntsIUTA", dir = "extdata")
  # setup <- setupProject(projPath = projPath, save = FALSE)
  # setup
  
  
  setup <- base::list()
  setup$projPath  <- projPath
  setup$date <- date
  
  #Create holder for list of samples
  sampleInfo <- base::data.frame(filePath = base::as.character(), sample = base::as.character(),
                                 group = base::as.character(), blank = base::as.character())
  
  
  if (convertFiles)
  {
    if (!base::is.null(convertFrom))
      ntsIUTA::mzMLconverter(path = projPath, convertFrom = convertFrom, centroidData = centroidData)
    else { stop("Vendors should be specified for file recognition. Possible entries are: thermo, bruker, agilent,
                 ab (from AB Sciex) and waters.") }
  }
  
  
  #Screen for MS files in project folder and add info to sampleInfo
  msFiles <- base::list.files(path = projPath, pattern = ".mzML|.mzXML", recursive = TRUE, full.names = TRUE, no.. = TRUE)
  
  if (base::length(msFiles) == 0)
  {
    warning("No mzML or mzXML files were found in selected project path. Use addFiles() to add files to project.")
    
    #Adds files info to sample data frame
  } else {
    sampleInfo[1:base::length(msFiles),] <- NA
    sampleInfo$filePath <- msFiles #base::dirname(msFiles)
    sampleInfo$sample <- tools::file_path_sans_ext(base::basename(msFiles))
    sampleInfo$blank <- blanks
    
    if (base::length(groups) < 2 & base::is.null(groups))
    {
      sampleInfo$group <- tools::file_path_sans_ext(base::basename(msFiles))
      #tentative to group samples and add blank group
      sampleInfo <- dplyr::mutate(sampleInfo, group = base::ifelse(base::grepl("A-r|qc|QC", sampleInfo$sample), "QC", group))
      sampleInfo <- dplyr::mutate(sampleInfo, group = base::ifelse(base::grepl("Blank|blank|B-r", sampleInfo$sample), "Blank", group))
      
    } else sampleInfo$group <- groups
    
    if (base::length(blanks) < 2 & base::is.null(blanks)) sampleInfo$blank <- base::ifelse("Blank" %in% sampleInfo$group, "Blank", NA)
    
    sampleInfo$polarity <- polarity
    
    setup$sampleInfo <- sampleInfo
  }
  
  
  if (save) # Add possibility to move the R project file 
  {
    rData <- base::paste0(projPath,"\\rData")
    
    if (base::dir.exists(rData))
    {
      base::saveRDS(setup, file = base::paste0(rData,"\\setup.rds"))
    } else {
      base::dir.create(rData)
      base::saveRDS(setup, file = base::paste0(rData,"\\setup.rds"))
    }
  }
  
  
  if (makeNewProject)
  {
    # add template for main script file
    sp <- base::file.path(projPath,"mainScript.R")
    base::cat(
      
      "#Copy template from template.R using x and/or use ?ntsIUTA for a tutorial.\n
#Run the following code to load the project setup:\n
library(ntsIUTA) #/n
setup <- readRDS('rData/setup.rds')",
      
      file = sp, sep = "")
    rstudioapi::initializeProject(projPath)
    rstudioapi::openProject(projPath, newSession = TRUE)
    base::print("Run  rstudioapi::navigateToFile('mainScript.R')  in the new project to open the mainScript.R file.")
  }
  
  return(setup)
  
}

