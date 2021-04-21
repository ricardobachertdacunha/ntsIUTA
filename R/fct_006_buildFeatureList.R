


#' @title buildFeatureList
#' @description Function to create a \code{data.frame} with detailed information for each feature in the given a \linkS4class{XCMSnExp} object object.
#' Optionally, a list of \linkS4class{xsAnnotate} objects per replicate group as obtained by \code{\link{makeFeatureComponents}} can be given.
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
#' @importFrom MSnbase fileNames
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
buildFeatureList <-  function(x = featData, xA = featComp, xPat = patData,
                              snWindow = 240,
                              save = TRUE, projPath = setup$projPath){
  
  maxMultiProcess = TRUE
  if (maxMultiProcess)
  {
    snow <- BiocParallel::registered("SnowParam")
    if (snow$workers < parallel::detectCores())
    {
      snow <- BiocParallel::SnowParam(workers = parallel::detectCores()-1, type = "SOCK", exportglobals = FALSE, progressbar = TRUE)
      BiocParallel::register(snow, default = TRUE)
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







