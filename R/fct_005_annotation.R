

#' @title AlteredCameraParam
#'
#' @slot sigma The multiplier of the standard deviation for grouping features
#' by retention time. 
#' @slot perfwhm Percentage of the overlapping FWHM to group features. 
#' @slot cor_eic_th Threshold for feature EIC correlation in each sample.
#' @slot calcCaS Calculate correlation across samples. 
#' @slot cor_exp_th Threshold for intensity correlation across samples.
#' @slot pval p-value threshold for testing correlation of significance.
#' @slot validateIsotopePatterns Logical, set to \code{TRUE} for validating
#' the annotated isotopes with the \emph{kegg} database.
#' @slot ppmIsotopes The expected mass deviation (in ppm) to find isotopes. 
#' @slot mzabs numeric. 
#' @slot noise numeric. 
#' @slot searchAdducts Logical, set to \code{TRUE} to screen for adducts after finding isotopes. 
#' @slot ppmAdducts The expected mass deviation (in ppm) to find adducts. 
#' @slot extendedList Logical, set to \code{TRUE} to use the extended list of
#' adducts. The default is \code{FALSE}.
#'
#' @return An \linkS4class{AlteredCameraParam} object containing parameters for
#' annotation of features in an \linkS4class{ntsData} object.
#' 
#' @export
#'
setClass("AlteredCameraParam",
  slots = c(
    sigma = "numeric",
    perfwhm = "numeric",
    cor_eic_th = "numeric",
    calcCaS = "logical",
    cor_exp_th = "numeric",
    calcIso = "logical",
    pval = "numeric",
    validateIsotopePatterns = "logical",
    ppmIsotopes = "numeric",
    mzabs = "numeric",
    noise = "numeric",
    searchAdducts = "logical",
    ppmAdducts = "numeric",
    extendedList = "logical"
  ))

AlteredCameraParam <- function(
  sigma = 5,
  perfwhm = 0.45,
  cor_eic_th = 0.85,
  calcCaS = TRUE,
  cor_exp_th = 0.85,
  pval = 0.05,
  validateIsotopePatterns = TRUE,
  ppmIsotopes = 50,
  mzabs = 0.01,
  noise = 350,
  searchAdducts = TRUE,
  ppmAdducts = 5,
  extendedList = FALSE) {
  
  paramobj <- do.call(new, c("AlteredCameraParam", as.list(environment())))
  
  paramobj
  
  return(paramobj)
  
}




#' @title annotateFeatures
#' 
#' @description Write
#'
#' @param obj An \linkS4class{ntsData} object with features.
#' @param algorithm The algorithm for finding isotopes. Possible values are
#' "alteredcamera", "camera" or "ramclustr".
#' @param param The param used for annotation of isotopes and adducts.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return
#' 
#' @export
#'
#' @importFrom checkmate checkSubset
#' @importFrom patRoon getXCMSSet
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr findAdducts getPeaklist
#' @importFrom utils txtProgressBar setTxtProgressBar read.table
#' @importFrom dplyr rename select everything
#'
annotateFeatures <- function(obj = NULL,
                             algorithm = NULL,
                             param = NULL,
                             save = FALSE
                             ) {
  
  if (is.null(obj)) return(cat("An ntsData object must be given!"))
  
  if (is.null(algorithm)) if (!is.na(obj@algorithms$annotation)) algorithm = obj@algorithms$annotation
  
  if (is.null(param)) if (length(obj@parameters$annotation) > 0) param = obj@parameters$annotation
  
  checkmate::checkSubset(algorithm, c("alteredcamera",
                                      "camera",
                                      "ramclustr",
                                      "nontarget",
                                      "intclust"))
  
  if (is.null(param)) return(cat("Parameters for annotation must be given!"))
  
  rg <- unique(obj@samples$group[!(obj@samples$group %in% obj@samples$blank)])
  
  if (!(algorithm == "alteredcamera")) {
    
    ag <- list(fGroups = obj@patdata, algorithm = algorithm)
    
    pat <- do.call(generateComponents, c(ag, param))
    
    # TODO convert to data frame from other functions and integrate with df scheme
    
  } else {
    
    xs <- obj@patdata
    
    xs <- patRoon::getXCMSSet(xs, verbose = TRUE, exportedData = TRUE)
    
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
      
      setTxtProgressBar(pb, ((rgidx/length(rg))*100))
      
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
        validateIsotopePatterns = validateIsotopePatterns
      )
      
      if (param@searchAdducts) {
        
        if (param@extendedList) {
          rules_pos <- system.file('rules/extended_adducts_pos.csv', package = "CAMERA")
          rules_neg <- system.file('rules/extended_adducts_neg.csv', package = "CAMERA")
        } else {
          rules_pos <- system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
          rules_neg <- system.file('rules/primary_adducts_neg.csv', package = "CAMERA")
        }
        
        if (obj@polarity == "positive") {
          rules <- utils::read.table(rules_pos, header = TRUE, sep = ",")
        }
        
        if (obj@polarity == "negative") {
          rules <- utils::read.table(rules_neg, header = TRUE, sep = ",")
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
  
  #retrieve isotopologues
  isotopes <- df$isotopes
  isotopes <- strsplit(isotopes, split = "\\]\\[")
  isogroup <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][1])
  isogroup <- str_extract(isogroup, pattern = "([0-9]+)")
  isogroup <- as.numeric(isogroup)
  isoclass <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][2])
  isoclass <- as.vector(ifelse(!is.na(isoclass), paste0("[",isoclass), NA))
  df$isogroup <- isogroup
  df$isoclass <- isoclass

  
  #assign Mions
  Mion <- df$mz
  index <- seq_len(length(Mion))
  isMion <- grepl(pattern = "[M]", isoclass, fixed = TRUE)
  charge <- sapply(index, FUN =  function(x) if (isMion[x]) as.numeric(str_extract(isoclass[x], pattern = "([0-9]+)")))
  charge[is.na(charge)] <- 1
  charge <- unlist(lapply(charge, function(x) ifelse(is.null(x), NA, x)))
  Mion <- Mion - 1.007276 / charge
  names(Mion) <- isogroup
  
  #assign Mions to isotopologues
  isopolog <- !isMion
  isopolog[is.na(isoclass)] <- FALSE
  Mion <- unlist(lapply(index, function(x) ifelse(isopolog[x], Mion[names(Mion) %in% isogroup[x]], Mion[x])))
  
  df$Mion <- Mion
  
  #assign adducts
  adducts <- strsplit(df$adduct, split = " ")
  df$adductclass <- as.character(lapply(X = seq_len(length(adducts)), function(x) adducts[[x]][1]))
  df$adductMion <- as.numeric(lapply(X = seq_len(length(adducts)), function(x) adducts[[x]][2]))
  
  #organise to end
  df <- rename(df, comp = pcgroup)
  df <- select(df, ID, group, comp, Mion, isoclass, isogroup, adductclass, adductMion, everything())
  
  obj@annotation$df <- df
  
  obj@annotation$algorithm <- algorithm
  
  obj@annotation$param <- param
  
  obj@annotation$raw <- xAL
  
  if (save) saveObject(obj = obj)
  
  if (is.character(save)) saveObject(obj = obj, filename = save)
  
  return(obj)
  
}
