


#' @title screenSuspectsFromCSV
#'
#' @param patData A \linkS4class{featureGroups} object with the feature groups of interest
#' @param polarity The polatiry of the MS files in \code{patData}. Possible values are 'positive' or 'negative'.
#' @param ppmWindow The allowed mass deviation, in ppm.
#' @param rtWindow The allowed retention time window, in seconds, when 'rt' is given in the screening list (\code{DB}).
#' @param screeningList A screeningList.
#' @param removeBlanks Logical, set to \code{TRUE} to ignore blank samples during suspect screening.
#' @param blankGroups A charcater vector with replicate group name of blanks to ignore.
#' @param withMS2 Logical, set to \code{TRUE} for using MS2 for confirming suspects.
#' @param listMS2 An MS2 list, if not given it will be calculated.
#' @param ppmMS2 Optional, sets a different mass deviation (in ppm) allowance for correlation of MS2 data.
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
                                  ppmWindow = 5, rtWindow = 10, screeningList = utils::read.csv(base::file.choose()),
                                  removeBlanks = TRUE, blankGroups = "Blank",
                                  withMS2 = TRUE, listMS2 = NULL, ppmMS2 = NULL) {
  
  #TODO Improve description of the function
  
  
  sDB <- screeningList
  sDB <- dplyr::filter(sDB, name != "")
  has_rt = FALSE
  if("rt" %in% base::colnames(sDB)) {has_rt <- TRUE}
  if(has_rt) {sDB$rt <- sDB$rt*60}
  
  if(base::length(polarity) > 1) stop("Multiple polarities are not supported! Select either pos for positive or neg for negative ionization.")
  if (polarity == "positive") adduct <- "[M+H]+"
  if (polarity == "negative") adduct <- "[M-H]-"
  
  if (base::is.null(ppmMS2)) ppmMS2 = ppmWindow
  
  rGroups <- patRoon::analysisInfo(patData)$group
  #if (removeBlanks) rGroups <- rGroups[rGroups != blankGroups]
  
  groupNames <- patRoon::replicateGroups(patData)
  if (removeBlanks) groupNames <- groupNames[!(groupNames %in% blankGroups)]
  
  df_patData <- base::as.data.frame(patRoon::as.data.table(patData, average = T))
  
  patSuspects <- patRoon::screenSuspects(patData, sDB, rtWindow = rtWindow, mzWindow = 0.03, adduct = adduct, onlyHits = TRUE)
  
  suspectsDT  <- dplyr::arrange(patRoon::screenInfo(patSuspects), group)
  suspectsDT  <- dplyr::select(suspectsDT, group, name, d_mz, d_rt)
  suspectsDT  <- dplyr::left_join(suspectsDT, df_patData, by = "group")
  suspectsDT <- dplyr::left_join(suspectsDT, sDB[,c("name", "formula", "comment", "int10")], by = "name")
  suspectsDT$d_ppm <- (base::abs(suspectsDT$d_mz)/suspectsDT$mz)*1E6
  suspectsDT <- dplyr::select(suspectsDT, group, name, formula, d_ppm, d_rt, mz, ret, dplyr::everything(), -d_mz)
  suspectsDT <- dplyr::filter(suspectsDT, d_ppm <= ppmWindow)
  suspectsDT <- base::as.data.frame(suspectsDT)
  
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
  
}
