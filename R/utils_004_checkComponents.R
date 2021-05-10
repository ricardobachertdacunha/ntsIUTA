


#' @title checkComponents
#' @description Inspect components by of given features, mz and retention time or component numbers of specified replicate groups
#' as defined in \code{\link{setupProject}}.
#'
#' @param xA The list of \linkS4class{xsAnnotate} objects obtained via \code{\link{makeFeatureComponents}}.
#' @param replicateGroups The replicate groups name or number to filter the \code{xA} list,
#' which is named according to the experimental replicate groups as defined in \code{\link{setupProject}}.
#' @param features Feature identifier/s.
#' @param featData The feature data to get the feature identifiers.
#' Can be a \linkS4class{XCMSnExp} or \linkS4class{featureGroups} object from \pkg{xcms} or \pkg{patRoon}, respectively.
#' @param mz The \emph{m/z} of interest. If \code{features} are specified \code{mz is not used}.
#' @param ppm The expected mass deviation to search for features of a given \code{mz}.
#' @param rt The expected retention time of the \emph{m/z} of interest, only used if \code{features} are not speficied.
#' @param rtWindow The expected retention time deviation for searching.
#' @param rtUnit The time unit used. Default is minutes.
#' @param comp The component numbers to plot. Only used if both \code{features} and \code{mz} are \code{NULL} (i.e. not specified).
#' @param onlyAnnotated Plots only annotated features.
#' @param onlyRelated Plots only relevant features (i.e. isotopes and addcuts) of the specified features. 
#' @param rGroupsAverage Logical, set to \code{TRUE} for averaging the intensities of each sample acording the defined sample replicate groups.
#' @param intval The intensity value type. Default is "maxo", corresponding to the height.
#'
#' @return A data.frame with detailed information about a component of interest.
#' 
#' @export
#'
#' @import magrittr
#' @importFrom patRoon as.data.table
#' @importFrom xcms featureSummary featureDefinitions
#' @importFrom CAMERA getPeaklist
#' @importFrom dplyr select everything filter between all_of
#' @importFrom stringr str_extract
#'
#' @examples
#' 
#' 
#' 
checkComponents <- function(xA = featComp, replicateGroups = NULL,
                            features = NULL, featData = featData,
                            mz = NULL, ppm = 5,
                            rt = NULL, rtWindow = 1, rtUnit = "min",
                            comp = NULL,
                            onlyAnnotated = FALSE, onlyRelated = TRUE,
                            rGroupsAverage = FALSE,
                            intval = "maxo") {
  
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
  # rtUnit = "min"
  # intval = "maxo"
  #rGroupsAverage = FALSE
  
  
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
  
  if (rGroupsAverage)
  {
    #check the col names for duplicates
    sampleNamesRepeats <- ifelse(xA[[1]]@xcmsSet$sample_name %in% xA[[1]]@xcmsSet$sample_group,
                                 paste0(xA[[1]]@xcmsSet$sample_name,".1"),xA[[1]]@xcmsSet$sample_name)
    
    #calculate averages for the intensites
    anoS <- dplyr::select(ano, -all_of(xA[[1]]@xcmsSet$sample_name), -all_of(xA[[1]]@xcmsSet$sample_group), -all_of(sampleNamesRepeats))
    
    for (r in rIndex) {
      anoS[,rGroups[r]] <- apply(ano[,colnames(ano) %in% sampleNamesRepeats[xA[[1]]@xcmsSet$sample_group %in% rGroups[r]], drop = F],
                                 MARGIN = 1, FUN = function(x) base::ifelse(log, log(mean(x)), mean(x))) 
    }
    
    return (base::as.data.frame(anoS))
    
  } else {
    
    return(base::as.data.frame(ano))
    
  }
}


