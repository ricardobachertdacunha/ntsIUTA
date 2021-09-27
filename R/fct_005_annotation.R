

#' annotateFeatures
#'
#' @description Group features into components according to co-elution and
#' EIC similarities. Then, it annotates isotopes and adducts in each component
#' per sample replicate group.
#'
#' @param obj An \linkS4class{ntsData} object with features.
#' @param algorithm The algorithm for finding isotopes.
#' Currently, "alteredcamera" in the only possible value.
#' @param param The param used for annotation of isotopes and adducts.
#' See ?AlteredCameraParam for more information.
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
#' @importFrom stringr str_extract
#'
annotateFeatures <- function(obj = NULL,
                             algorithm = NULL,
                             param = NULL,
                             excludeBlanks = FALSE,
                             save = FALSE
                             ) {

  assertClass(obj, "ntsData")

  if (is.null(algorithm)) algorithm <- annotationParameters(obj)@algorithm

  if (is.na(algorithm)) {
    warning("Annotation algorihtm not defined!")
    return(obj)
  }

  if (is.null(param)) param <- annotationParameters(obj)@param

  if (length(param) == 0) {
    warning("Parameters for annotation not found!
            Use the function AlteredCameraParam to obtain a deafult list of parameters.")
    return(obj)
  }

  if (excludeBlanks) {
    rg <- unique(obj@samples$group[!(obj@samples$group %in% obj@samples$blank)])
  } else {
    rg <- unique(obj@samples$group)
  }

  if (algorithm == "alteredcamera") {
    
    if (is.list(param)) param <- param[[1]]

    xs <- obj@patdata

    xs <- getXCMSSet(xs, verbose = TRUE, loadRawData = TRUE)

    xAL <- list()

    pb <- txtProgressBar(min = 0, max = 100, initial = 0, char = "=", width = 80, style = 3)

    setTxtProgressBar(pb, 0)

    for (rgidx in seq_len(length(rg))) {

      sampleidxs <- which(obj@samples$group == rg[rgidx])

      xA_temp <- xsAnnotate(xs = xs[, sampleidxs],
                            sample = c(1:length(sampleidxs)),
                            polarity = obj@polarity)

      xA_temp <- groupFWHM(xA_temp, sigma = param@sigma,
                           perfwhm = param@perfwhm,
                           intval = "maxo")
      
      xA_temp <- groupCorr(xA_temp,
                           cor_eic_th = param@cor_eic_th,
                           cor_exp_th = param@cor_exp_th,
                           pval = param@pval,
                           graphMethod = "hcs",
                           calcIso = FALSE,
                           calcCiS = TRUE,
                           calcCaS = FALSE,
                           psg_list = NULL,
                           xraw = NULL,
                           intval = "maxo")
      
      # TODO check and test isotopic detection, specially for one sample
      xA_temp <- FindIsotopesWithValidationAltered(
        xA = xA_temp,
        obj = obj,
        sampleidxs = sampleidxs,
        ppm = param@ppmIsotopes,
        noise = param@noise,
        maxcharge = 3,
        intval = "maxo",
        validateIsotopePatterns = param@validateIsotopePatterns
      )

      # df_temp <- CAMERA::getPeaklist(xA_temp2, intval = "maxo")
      # df_temp$group <- rg[rgidx]
      # df_temp$ID <- unique(obj@peaks$feature[obj@peaks$sample %in% samples(obj)[sampleidxs]])
      # #df_temp$ID <- obj@features$ID
      # df_temp <- df_temp[df_temp$pcgroup == 18, ]
      # View(df_temp)
      
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
          # TODO add deprotonated adduct directly in the rule table
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

      setTxtProgressBar(pb, ((rgidx / length(rg)) * 100))

    }

    
    #prepare data frame
    for (r in seq_len(length(xAL))) {
      df_temp <- CAMERA::getPeaklist(xAL[[r]], intval = "maxo")
      df_temp$group <- names(xAL)[r]
      
      #for one sample, peaks might be missing
      sp <- samples(obj)[sampleGroups(obj) == names(xAL)[r]]
      
      if (length(sp) > 1) {
        
        df_temp$ID <- obj@features$ID
        df_temp$intensity <- apply(df_temp[, gsub("[-.]", "", colnames(df_temp)) %in% gsub("[-.]", "", samples(obj))], MARGIN = 1, function(x) mean(x, na.rm = TRUE))
        df_temp$intensity_sd <- apply(df_temp[, gsub("[-.]", "", colnames(df_temp)) %in% gsub("[-.]", "", samples(obj))], MARGIN = 1, function(x) sd(x, na.rm = TRUE))
        df_temp$intensity[is.na(df_temp$intensity)] <- 0
        df_temp$intensity_sd[is.na(df_temp$intensity_sd)] <- 0
        df_temp <- df_temp[, !gsub("[-.]", "", colnames(df_temp)) %in% gsub("[-.]", "", samples(obj))]
        df_temp <- df_temp[, !colnames(df_temp) %in% rg]
        
      } else {
        
        df_temp$ID <- obj@features$ID[obj@features[, names(xAL)[r], drop = TRUE] > 0]
        df_temp <- rename(df_temp, intensity = maxo)
        df_temp$intensity_sd <- 0
        df_temp$intensity[is.na(df_temp$intensity)] <- 0
        df_temp <- select(df_temp, -into, -sample)
        df_temp$npeaks <- 1
        df_temp <- select(df_temp, mz, mzmin, mzmax, rt, rtmin, rtmax, npeaks, isotopes, adduct, pcgroup, group,  ID, intensity, intensity_sd)

      }
      
      if (r == 1 | length(xAL) == 1) {
        df <- df_temp
      } else {
        df <- rbind(df, df_temp)
      }

    }
    df$pcgroup <- paste0(df$pcgroup, "_", df$group)
    df <- rename(df, comp = pcgroup)


  } else {
    # TODO Addapt to other algorithms
    ag <- list(fGroups = obj@patdata, algorithm = algorithm)
    pat <- do.call(generateComponents, c(ag, param))
  }


  ## Isotopic information for each feature per sample replicate group
  isotopes <- df$isotopes
  isotopes <- strsplit(isotopes, split = "\\]\\[")

  isogroup <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][1])
  isogroup <- str_extract(isogroup, pattern = "([0-9]+)")
  noIsogroup <- is.na(isogroup)
  isogroup <- paste0(isogroup, "_", df$group)
  isogroup[noIsogroup] <- NA

  isoclass <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][2])
  isoclass <- as.vector(ifelse(!is.na(isoclass), paste0("[",isoclass), NA))

  isonr <- gsub("\\[M\\]\\+|\\[M\\]\\-", 0, isoclass)
  isonr <- as.numeric(str_extract(isonr, "[0-9]"))

  df$isogroup <- isogroup
  df$isoclass <- isoclass
  df$isonr <- isonr


  ## Make Mions
  # When not annotated, charge is 1 by default and respective Mion is calculated
  Mion <- df$mz
  index <- seq_len(length(Mion))
  isMion <- grepl(pattern = "[M]", isoclass, fixed = TRUE)

  charge <- sapply(index, FUN =  function(x) {
   ifelse(isMion[x], as.numeric(str_extract(isoclass[x], pattern = "([0-9]+)")), NA)
  })

  charge[is.na(charge)] <- 1
  if (obj@polarity == "positive") {
    baseadduct <- "[M+H]+"
    Mion <- (Mion - 1.007276) * charge
  } else {
    baseadduct <- "[M-H]-"
    Mion <- (Mion + 1.007276) * charge
  }
  Mion <- round(Mion, digits = 3)


  ## Collect possible adducts and add conflicts information
  adducts <- strsplit(df$adduct, split = " ")
  conflicts <- unlist(lapply(X = seq_len(length(adducts)), function(x) ifelse(length(adducts[[x]]) > 2, TRUE, FALSE)))

  adductclass <- lapply(X = seq_len(length(adducts)), function(x) grep("M", adducts[[x]], value = TRUE))
  adductMion <- lapply(X = seq_len(length(adducts)), function(x) as.numeric(grep("M", adducts[[x]], value = TRUE, invert = TRUE)))

  possibleAdducts <- unlist(lapply(adductMion, length))
  nAdducts <- possibleAdducts
  isAdduct <- possibleAdducts != 0


  ## Updates Mion with non-conflicting adducts
  Mion <- unlist(lapply(index, FUN = function(x) ifelse(isAdduct[x] & !conflicts[x], as.numeric(adductMion[[x]]), Mion[x])))

  df$conflicts <- conflicts
  df$nAdducts <- nAdducts

  adductclass <- lapply(adductclass, function(x) {
    y <- x
    if (length(x) == 0) y <- NA
    return(y)
  })

  adductMion <- lapply(adductMion, function(x) {
    y <- x
    if (length(x) == 0) y <- NA
    return(y)
  })

  df$adductMion <- I(adductMion)

  df$adductclass <- I(adductclass)


  ## solve conflicts
  # Rules
  # - Give priority to [M+H]+
  # - When ties between adducts return the one with lowest mass error

  consolidate <- lapply(index, FUN =  function(x, dt, adductclass) {

    hasH <- numeric()

    hasP <- numeric()

    y <- adductclass[[x]]

    if (length(y) > 1) {

      hold <- dt[x, ]

      #priority
      hasH <- grep("H]", y, fixed = TRUE)

      if (length(hasH) > 0) {

        df$adductMion[x] <<- unlist(hold$adductMion)[hasH]
        df$adductclass[x] <<- y[hasH]
        Mion[x] <<- unlist(hold$adductMion)[hasH]

      } else {

        #when priority adduct not present, calculate mass error and return the lowest
        massrules <- rules[rules$name %in% y, ]
        massrules <- massrules[match(y, massrules$name), ]
        errors <- data.frame(mions = unlist(hold$adductMion),
                             adu = y,
                             mz = rep(hold$mz, length(y)),
                             massdiff = massrules$massdiff)
        errors$expect <- errors$mz - errors$massdiff
        errors$error <- abs(errors$mions - errors$expect)

        hasP <- min(errors$error, na.rm = TRUE)
        hasP <- which(errors$error == hasP)

        if (length(hasP) < 2) {
          df$adductMion[x] <<- unlist(hold$adductMion)[hasP]
          df$adductclass[x] <<- y[hasP]
          Mion[x] <<- unlist(hold$adductMion)[hasP]
        } else {
          #tied case for same mass error, leave it NA with default Mion
          df$adductMion[x] <<- NA
          df$adductclass[x] <<- NA
        }
      }
    }
  }, adductclass = adductclass, dt = df)


  ## Update Mion with isotopologues after settling adducts priority when multiple hits
  isopolog <- !isMion
  isopolog[is.na(isoclass)] <- FALSE

  names(Mion) <- isogroup

  Mion <- unlist(lapply(index, function(x) ifelse(isopolog[x], Mion[names(Mion) %in% isogroup[x]], Mion[x])))

  df$Mion <- as.numeric(Mion)


  ## finishing data.frame and adding info to obj
  df <- select(df, ID, group, comp, Mion, isonr, isoclass, isogroup, adductclass, adductMion, conflicts, nAdducts, everything())

  obj@annotation$comp <- df

  ## TODO make componentsFeatures to integrate in patRoon workflow

  obj@annotation$raw <- xAL

  obj <- annotationParameters(obj, algorithm = algorithm, param = param)
  
  if (save) saveObject(obj = obj)

  if (is.character(save)) saveObject(obj = obj, filename = save)

  return(obj)

}
