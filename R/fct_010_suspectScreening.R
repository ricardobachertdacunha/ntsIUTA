


#' @title screenSuspectsFromCSV
#'
#' @param patData A \linkS4class{featureGroups} object with the feature groups of interest
#' @param polarity The polatiry of the MS files in \code{patData}. Possible values are 'positive' or 'negative'.
#' @param ppm The allowed mass deviation, in ppm.
#' @param rtWindow The allowed retention time window, in seconds, when 'rt' is given in the screening list (\code{DB}).
#' @param DB A screeningList.
#' @param withMS2 Logical, set to \code{TRUE} for using MS2 for confirming suspects.
#' @param listMS2 An MS2 list, if not given it will be calculated.
#'
#' @return
#' 
#' @references \insertRef{Helmus2021}{ntsIUTA}
#' 
#' @export
#' 
#' @importFrom dplyr filter arrange left_join select everything mutate distinct desc
#' @importFrom utils read.csv head
#' @importFrom patRoon screenSuspects as.data.table screenInfo getDefAvgPListParams replicateGroups generateMSPeakLists generateFormulasSIRIUS
#' @importFrom fuzzyjoin difference_inner_join
#'
#' @examples
#' 
#' 
#' 
screenSuspectsFromCSV <- function(patData = patData, polarity = "positive",
                                  ppm = 10, rtWindow = 15, DB = utils::read.csv(base::file.choose()),
                                  withMS2 = FALSE, listMS2 = NULL) {
  
  #TODO Improve description of the function
  
  sDB <- DB
  sDB <- dplyr::filter(sDB, name != "")
  has_rt = FALSE
  if("rt" %in% base::colnames(sDB)) {has_rt <- TRUE}
  if(has_rt) {sDB$rt <- sDB$rt*60}
  
  if(base::length(polarity) > 1) stop("Multiple polarities are not supported! Select either pos for positive or neg for negative ionization.")
  if (polarity == "positive") adduct <- "[M+H]+"
  if (polarity == "negative") adduct <- "[M-H]-"
  
  patSuspects <- patRoon::screenSuspects(patData, sDB, rtWindow = rtWindow, mzWindow = 0.03, adduct = adduct, onlyHits = TRUE)
  
  suspectsDT <- dplyr::arrange(patRoon::as.data.table(patSuspects, average = T), group)
  suspectsDT <- dplyr::left_join(suspectsDT, sDB[,c("name", "formula", "comment", "int10")], by = "name")
  suspectsDT  <- dplyr::left_join(suspectsDT, dplyr::select(dplyr::arrange(patRoon::screenInfo(patSuspects), group), group, d_mz, d_rt), by = "group")
  suspectsDT$d_ppm <- (base::abs(suspectsDT$d_mz)/suspectsDT$mz)*1E6
  suspectsDT <- dplyr::select(suspectsDT, name, formula, d_ppm, d_rt, mz, ret, group, dplyr::everything(), -d_mz)
  suspectsDT <- dplyr::filter(suspectsDT, d_ppm <= ppm)
  suspectsDT <- base::as.data.frame(suspectsDT)
  
  if (withMS2)
  {
  
    if(base::is.null(listMS2)) listMS2 <- base::list()
    MS2_avgPListParams <- patRoon::getDefAvgPListParams(clusterMzWindow = 0.005,
                                               topMost = 50,
                                               minIntensityPre = 10,
                                               minIntensityPost = 10)
    
    # Evaluation and categorization of hits per replicate group
    for (i in 1:base::length(patRoon::replicateGroups(patData))) {
      
      #Collecting name and prepare table
      colMS2 <- base::paste0("ms2_",patRoon::replicateGroups(patData)[i])
      colCat <- base::paste0("cat_",patRoon::replicateGroups(patData)[i])
      suspectsDT[, colMS2] <- 0
      suspectsDT[, colCat] <- 0
      
      #Check for lowest categories, 4 only mass and 3 mass and rt
      temp <- suspectsDT[,c("name", "group", "mz", "ret", "d_rt", patRoon::replicateGroups(patData)[i]), drop = F]
      suspectsDT[, colCat] <- base::ifelse(base::as.numeric(temp[,patRoon::replicateGroups(patData)[i]]) > 0, 4, 0)
      suspectsDT[, colCat] <- base::ifelse(base::as.numeric(temp[,patRoon::replicateGroups(patData)[i]]) > 0 & !base::is.na(temp$d_rt), 3, suspectsDT[, colCat])
      
      #Extracts the MS2 data for each non 0 feature
      if(!base::is.null(listMS2[[patRoon::replicateGroups(patData)[i]]])){
        tempMS2 <- listMS2[[patRoon::replicateGroups(patData)[i]]]
      } else {
        tempMS2 <- temp[temp[,patRoon::replicateGroups(patData)[i]] > 0,]
        tempMS2 <-base::suppressWarnings(
          patRoon::generateMSPeakLists(patData[base::which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), tempMS2$group[drop = TRUE]],
                                                       "mzr", maxMSRtWindow = 5,
                                                       precursorMzWindow = 2,
                                                       avgFeatParams = MS2_avgPListParams,
                                                       avgFGroupParams = MS2_avgPListParams)
        )
        listMS2[[patRoon::replicateGroups(patData)[i]]] <- tempMS2
      }
      
      #Loop for all the rows in suspectsDT
      for (j in 1:base::nrow(temp)) {
        
        if (base::as.numeric(temp[j, patRoon::replicateGroups(patData)[i], drop = T]) > 0) {
          
          temp1 <- tempMS2[[temp$group[j]]]$MSMS
          
          #tentative to extract MS2 again if temp1 is NULL
          if(base::is.null(temp1)) {
            lastTemptMS2 <- base::suppressWarnings(
              patRoon::generateMSPeakLists(patData[which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), temp$group[j, drop = TRUE]],
                                  "mzr", maxMSRtWindow = 5,
                                  precursorMzWindow = 2,
                                  avgFeatParams = MS2_avgPListParams,
                                  avgFGroupParams = MS2_avgPListParams))
            temp1 <- lastTemptMS2[[temp$group[j]]]$MSMS
          }
          
          
          if (!base::is.null(temp1)) {
            
            temp1 <- dplyr::mutate(temp1, into_ind = intensity/base::max(temp1$intensity)) #normalize
            temp1 <- dplyr::filter(temp1, mz < temp$mz[j]+(5/1E6*temp$mz[j])) #remove mz higher than precursor, probably contamination
            
            # load MS2 from DB
            temp2 <- dplyr::filter(sDB, name == temp$name[j])
            if (base::is.na(temp2$hasMS2[drop = T])) {temp2$hasMS2 <- FALSE}
            
            if (temp2$hasMS2[drop = T]) {
              temp3 <- temp2$mzMS2[drop = T]
              temp3 <- base::as.data.frame(base::unlist(base::strsplit(temp3, split=";")))
              base::colnames(temp3) <- c("mz")
              temp4 <- temp2$intMS2[drop = T]
              temp4 <- base::as.data.frame(base::unlist(base::strsplit(temp4, split=";")))
              base::colnames(temp4) <- c("into_ind")
              temp5 <- base::cbind(temp3, temp4)
              temp5$mz <- base::as.numeric(base::as.character(temp5$mz))
              temp5$into_ind <- base::as.numeric(base::as.character(temp5$into_ind))
              #remove precurssor ion from fragments list
              temp5 <- dplyr::filter(temp5, mz < temp2$pos[drop=T] - (5/1E6*temp2$pos[drop=T]))
              
              # select top 5 in fragments from db, or top 10 if number of fragments is above 15
              if (base::nrow(temp5) < 15){
                top5_temp5 <- utils::head(dplyr::arrange(temp5, dplyr::desc(into_ind)), n = 5)
              } else { top5_temp5 <- utils::head(dplyr::arrange(temp5, dplyr::desc(into_ind)), n = 10) }
              
              # test match for MS2 correlation and excludes in diff > +/-5ppm
              temp6 <- fuzzyjoin::difference_inner_join(temp1, temp5, by = c("mz"), max_dist = 0.1, distance_col = "d_ppm")
              temp6$d_ppm <- base::abs(temp6$d_ppm)/temp6$mz.x*1E6
              temp6 <- dplyr::filter(temp6, d_ppm <= 10)
              
              # test match only in top 5
              top5_temp6 <- fuzzyjoin::difference_inner_join(temp1, top5_temp5, by = c("mz"), max_dist = 0.1, distance_col = "d_ppm")
              top5_temp6$d_ppm <- base::abs(top5_temp6$d_ppm)/top5_temp6$mz.x*1E6
              top5_temp6 <- dplyr::filter(top5_temp6, d_ppm <= 10)
              
              temp6 <- dplyr::distinct(temp6, mz.x, .keep_all= TRUE) # remove double entries for mz.x
              top5_temp6 <- dplyr::distinct(top5_temp6, mz.x, .keep_all= TRUE) # remove double entries for mz.x
              
              suspectsDT[j, colMS2] <- base::paste0(base::nrow(top5_temp6),"(",base::nrow(top5_temp5),")")
              
              if (base::nrow(top5_temp6) >= 2) {
                suspectsDT[j, colCat] <- 1
              } else {
                groupAndIsos <- dplyr::filter(df_patData, mz >= temp$mz[j, drop=T] & mz < 6+temp$mz[j, drop=T])
                groupAndIsos <- dplyr::filter(groupAndIsos, ret >= temp$ret[j, drop=T]-30 & ret <= temp$ret[j, drop=T]+30)
                tempMS2withIsotopes <- base::suppressWarnings(patRoon::generateMSPeakLists(patData[base::which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), groupAndIsos$group[drop = TRUE]], "mzr", maxMSRtWindow = 5, precursorMzWindow = 2, avgFeatParams = MS2_avgPListParams, avgFGroupParams = MS2_avgPListParams))
                #tentative to identify by insillico fragmentation
                tempFormulas <- patRoon::generateFormulasSIRIUS(patData[base::which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), groupAndIsos$group[drop = TRUE]],
                                                                tempMS2withIsotopes, relMzDev = 10,
                                                                adduct = "[M+H]+", elements = gsub("[^a-zA-Z]", "", temp2$formula[drop=T]), profile = "qtof",
                                                                database = NULL, noise = NULL, topMost = 20, extraOptsGeneral = NULL,
                                                                calculateFeatures = TRUE, featThreshold = 1)
                formulaResult <- tempFormulas[[temp$group[j, drop=T]]]
                formulaResult <- dplyr::filter(formulaResult, neutral_formula == temp2$formula[drop=T])
                
                suspectsDT[j, colMS2] <- base::paste0(base::nrow(formulaResult),"(",base::nrow(temp1),")")
                
                if (base::nrow(formulaResult) >= 2) {suspectsDT[j, colCat] <- 2}
              }
              
            } else {
              
              groupAndIsos <- dplyr::filter(df_patData, mz >= temp$mz[j, drop=T] & mz < 6+temp$mz[j, drop=T])
              groupAndIsos <- dplyr::filter(groupAndIsos, ret >= temp$ret[j, drop=T]-30 & ret <= temp$ret[j, drop=T]+30)
              tempMS2withIsotopes <-base::suppressWarnings(patRoon::generateMSPeakLists(patData[base::which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), groupAndIsos$group[drop = TRUE]], "mzr", maxMSRtWindow = 5, precursorMzWindow = 2, avgFeatParams = MS2_avgPListParams, avgFGroupParams = MS2_avgPListParams))
              temp2 <- dplyr::filter(sDB, name == temp$name[j])
              #tentative to identify by insillico fragmentation
              tempFormulas <- patRoon::generateFormulasSIRIUS(patData[base::which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), groupAndIsos$group[drop = TRUE]],
                                                              tempMS2withIsotopes, relMzDev = 10,
                                                              adduct = "[M+H]+", elements = base::gsub("[^a-zA-Z]", "", temp2$formula[drop=T]), profile = "qtof",
                                                              database = NULL, noise = NULL, topMost = 20, extraOptsGeneral = NULL,
                                                              calculateFeatures = TRUE, featThreshold = 1)
              
              # tempFormulas <- patRoon::generateFormulasGenForm(patData[which(patSampleDT$group == patRoon::replicateGroups(patData)[i]), groupAndIsos$group[drop = TRUE]],
              #                                                 tempMS2withIsotopes, relMzDev = 10, isolatePrec = TRUE,
              #                                                 adduct = "[M+H]+", elements = gsub("[^a-zA-Z]", "", temp2$formula[drop=T]),
              #                                                 topMost = 20, extraOpts = NULL,
              #                                                 calculateFeatures = TRUE, featThreshold = 1, timeout = 240, hetero = TRUE, oc = TRUE)
              
              formulaResult <- tempFormulas[[temp$group[j, drop=T]]]
              formulaResult <- dplyr::filter(formulaResult, neutral_formula == temp2$formula[drop=T])
              
              suspectsDT[j, colMS2] <- base::paste0(base::nrow(formulaResult),"(",base::nrow(temp1),")")
              
              if (base::nrow(formulaResult) >= 2) {suspectsDT[j, colCat] <- 2}
            }
            
          }
        }
      }
      rm(j, temp1, temp2, temp3)
    }
  }
  
  suspects <- base::list(patSuspects = patSuspects, suspectsDT = suspectsDT)
  
  if (!base::is.null(listMS2)) suspects[["MS2"]] <- listMS2
  
  return(suspects)
  
}
