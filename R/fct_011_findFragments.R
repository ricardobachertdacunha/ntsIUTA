

### ntsFragmentsData -----

#' @title ntsFragmentsData
#'
#' @slot targets A data frame with target values for finding fragments in MS2 data of a given feature list.
#' @slot param The parameters used for find fragments.
#' @slot data Extra data produced during screening for fragments.
#' @slot results A data.frame with summarized results per sample replicate group.
#'
#' @return An \linkS4class{ntsFragmentsData} object to be added to
#' the workflows slot of an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("ntsFragmentsData",
  slots = c(
    targets = "data.frame",
    param = "list",
    data = "list",
    results = "data.frame"
  ),
  prototype = list(
    targets = data.frame(),
    param = list(ppm = 10,
                 intMin = 10,
                 replicates = NULL,
                 ID = NULL),
    data = list(),
    results = data.frame()
  )
)

#' findFragments
#'
#' @param obj An \linkS4class{ntsData} object with features for finding fragment in MS2 data.
#' @param targets A data frame with target values for matching fragments.
#' See details for more information.
#' @param replicates A character vector with the name of the sample replicate groups to use for finding fragments.
#' @param ID A character vector with the ID of the features of interest.
#' @param title Optional, title to use for the entry in the workflows slot of the \linkS4class{ntsData} object.
#' @param ppm The mass deviation, in ppm, for matching fragments.
#' @param intMin The minimum intensity of fragments.
#' @param MS2param A list of parameters to use for extracting MS2 data.
#'
#' @return An \linkS4class{ntsData} object with an \linkS4class{ntsFragmentsData} object
#' add to the workflows slot.
#' 
#' @export
#'
findFragments <- function(obj, targets = NULL, replicates = NULL, ID = NULL,
                          title = NULL, ppm = 10, intMin = 10,
                          MS2param = NULL) {
  
  assertClass(obj, "ntsData")
  
  obj2 <- obj
  
  adduct <- ifelse(obj@polarity == "positive", "[M+H]+", ifelse(obj@polarity == "negative", "[M-H]-", NULL))
  if (is.null(adduct)) {
    warning("Polarity of ntsData not recognized!")
    return(obj)
  }
  
  
  
  ### Check Parameters -----
  
  data <- new("ntsFragmentsData")
  if (is.null(replicates)) replicates <- data@param$replicates
  if (is.null(ID)) ID <- data@param$ID
  if (is.null(ppm)) ppm <- data@param$ppm
  if (is.null(intMin)) intMin <- data@param$intMin
  if (is.null(targets)) targets <- data@targets
  if (is.null(MS2param)) MS2param <- obj@parameters@MS2 
  
  ### Collect Features -----
  
  #filter files for selected sampleGroups
  if (!is.null(replicates)) obj2 <- filterFileFaster(obj2, which(sampleGroups(obj2) %in% replicates))
  
  #select features of interest
  if (!is.null(ID)) obj2 <- obj2[, obj2@features$ID[obj2@features$ID %in% ID]]
  
  rg <- unique(sampleGroups(obj2))
  
  extra_data <- list()
  
  if (nrow(obj2@features) == 0) {
    warning("Features not found in the ntsData with the given conditions!")
    return(obj)
  }
  
  for (r in seq_len(length(rg))) {
  
    ft_org <- obj2@features
    
    ft <- dplyr::select(ft_org, ID, mz, rt, dppm, width, npeaks)
    ft[, c("int", "int_sd")] <- ft_org[, grep(rg[r], colnames(ft_org))]
    ft$group <- rg[r]
    ft <- dplyr::select(ft, group, ID, everything())
    ft <- ft[ft$int > 0, ]
    
    ### Load Targets -----
    
    targets_frag <- dplyr::filter(targets, process == "Fragment")
    targets_frag <- dplyr::mutate(targets_frag, MD = -MD)
    
    targets_cleav <- dplyr::filter(targets, process == "Cleavage")
    
    
    
    ### Load MS2 -----
    
    if (is.null(MS2param)) MS2param <- obj2@parameters@MS2
    
    MS2 <- extractMS2(obj2@patdata[which(sampleGroups(obj2) %in% rg[r]), ft$ID], param = MS2param)
    
    ft$hasFragments <- FALSE
    
    for (w in seq_len(nrow(ft))) if (!is.null(MS2[[ft$ID[w]]]$MSMS)) ft$hasFragments[w] <- TRUE
  
    
    
    ### Prepare Feature df -----
    
    ft <- ft[ft$hasFragments, ]
    
    ft <- dplyr::mutate(ft, hits = 0, ms2P = NA_character_, names = NA_character_, error_ppm = NA_character_)
    
    ft_info <- list()

    
    
    ### Find MS2 matches -----
    
    for (i in seq_len(nrow(ft))) {
      
      idf <- ft$ID[i]
      
      t <- MS2[[idf]]$MSMS
      
      # filter with for intensity thres
      t <- dplyr::filter(t, intensity > intMin)
      
      if (nrow(t) > 0) {
        
        if(nrow(t) > 1) { #when has more than one fragment
          
          #make full table for iteration
          t2 <- t
          t2 <- merge.data.frame(t,t2, by = NULL)
          t2 <- dplyr::mutate(t2, TP = "")
          
          # include precurssor if not present in MS2
          if (TRUE %in% t$precursor) {
            tmz <- t[1,]
            tmz$mz <- ft$mz[i]
            tmz$intensity <- "Fragmented"
            tmz$precursor <- TRUE
            tmz$ID <- max(t$ID) + 1
            tmz <- merge.data.frame(t, tmz, by = NULL)
            tmz <- dplyr::mutate(tmz, TP = "")
            t2 <- rbind(t2, tmz)
            ft$ms2P[i] <- "Added"
          } else {
            ft$ms2P[i] <- "Present"
          }
          
          # calculate differences between fragments
          t2 <- transform(t2, TP = as.numeric(TP))
          for (j in seq_len(nrow(t2))) t2$TP[j] <- t2$mz.y[j] - t2$mz.x[j]
          
          # include possible protonated induced cleavage
          t3 <- cbind(t, t)
          colnames(t3) <- c("ID.x", "mz.x", "intensity.x", "precursor.x", "ID.y", "mz.y", "intensity.y", "precursor.y")
          t3 <- dplyr::mutate(t3, TP = mz.x - 1.0073, precursor.x = "Direct", precursor.y = "Direct")
          t2 <- rbind(t2, t3)
          t2 <- dplyr::filter(t2, TP > 0)
          
          # search for MS2 neutral loss between fragments and possible protonated induced cleavage
          tTP <- fuzzyjoin::difference_inner_join(t2, targets_frag, by = c("TP" = "MD"), max_dist = 0.01, distance_col = "diff")
          tTP$diff <- (tTP$diff/abs(tTP$MD)) * 1E6
          tTP <- tTP[tTP$diff < ppm, ]
          
          # search for alpha cleavage fragments in MS2 list
          t4 <- cbind(t, t)
          colnames(t4) <- c("ID.x", "mz.x", "intensity.x", "precursor.x", "ID.y", "mz.y", "intensity.y", "precursor.y")
          t4 <- dplyr::mutate(t4, TP = mz.x, precursor.x = "Alpha", precursor.y = "Alpha")
          t4 <- dplyr::filter(t4, TP > 0)
          
          tTP2 <- fuzzyjoin::difference_inner_join(t4, targets_cleav, by = c("TP" = "MD"), max_dist = 0.01, distance_col = "diff")
          tTP2$diff <- (tTP2$diff/abs(tTP2$MD)) * 1E6
          tTP2 <- tTP2[tTP2$diff < ppm, ]
          
          # merge two resulting hit tables
          tTP <- rbind(tTP, tTP2)
          
        } else { #when only one fragment is found
          
          # Include the precussor if the MS2 is different
          if (TRUE %in% t$precursor) {
            tmz <- t[1,]
            tmz$mz <- ft$mz[i]
            tmz$intensity <- "Fragmented"
            tmz$precursor <- TRUE
            tmz$ID <- max(t$ID) + 1
            tmz <- merge.data.frame(t, tmz, by = NULL)
            tmz <- dplyr::mutate(tmz, TP = NA)
            tmz <- transform(tmz, TP = as.numeric(TP))
            for (j in 1:nrow(tmz)) tmz$TP[j] <- tmz$mz.y[j] - tmz$mz.x[j]
            
            # Include possible protonated induced cleavage
            t5 <- cbind(t, t)
            colnames(t5) <- c("ID.x", "mz.x", "intensity.x", "precursor.x", "ID.y", "mz.y", "intensity.y", "precursor.y")
            t5 <- dplyr::mutate(t5, TP = mz.x - 1.0073, precursor.x = "Direct", precursor.y = "Direct")
            tmz <- rbind(tmz,t5)
            tmz <- dplyr::filter(tmz, TP > 0)
            
            # Search for MS2 neutral loss between fragments and possible protonated induced cleavage
            tTP <- fuzzyjoin::difference_inner_join(tmz, targets_frag, by = c("TP" = "MD"), max_dist = 0.01, distance_col = "diff")
            tTP$diff <- (tTP$diff/abs(tTP$MD)) * 1E6
            tTP <- tTP[tTP$diff < ppm, ]
            
            ft$ms2P[i] <- "Single"
            
            
            # Search for alpha cleavage fragments in MS2 list
            t6 <- cbind(t, t)
            colnames(t6) <- c("ID.x", "mz.x", "intensity.x", "precursor.x", "ID.y", "mz.y", "intensity.y", "precursor.y")
            t6 <- dplyr::mutate(t6, TP = mz.x, precursor.x = "Alpha", precursor.y = "Alpha")
            t6 <- dplyr::filter(t6, TP > 0)
            
            tTP2 <- fuzzyjoin::difference_inner_join(t6, targets_cleav, by = c("TP" = "MD"), max_dist = 0.01, distance_col = "diff")
            tTP2$diff <- (tTP2$diff/abs(tTP2$MD)) * 1E6
            tTP2 <- tTP2[tTP2$diff < ppm, ]
            
            tTP <- rbind(tTP, tTP2)
            
          } else {
            
            ft$ms2P[i] <- "Only"
            
          }
        }
        
        # Writing of results in ft table
        if (exists("tTP")) {
          
          tTP <- dplyr::select(tTP, -process)
          if (nrow(tTP) > 0) {
            
            for (k in seq_len(nrow(tTP))) {
              if (tTP$mz.x[k] != tTP$mz.y[k]) { tTP$precursor.x[k] <- "Loss" }
            }
            
            ft$hits[i] <- nrow(tTP)
            ft$names[i] <- paste(as.vector(tTP$name), collapse = "; ")
            ft$error_ppm[i] <- paste(round(as.vector(tTP$diff), digits = 1), collapse = "; ")
            if (ft$hits[i] > 0) ft_info[[idf]] <- tTP
            
          }
        }
      }
      
      if (exists("t")) {rm(t)}
      if (exists("tmz")) {rm(tmz)}
      if (exists("tTP")) {rm(tTP)}
      if (exists("tTP2")) {rm(tTP2)}
      if (exists("t2")) {rm(t2)}
      if (exists("t3")) {rm(t3)}
      if (exists("t4")) {rm(t4)}
      if (exists("t5")) {rm(t5)}
      if (exists("t6")) {rm(t6)}
      if (exists("idf")) {rm(idf)}
      if (exists("j")) {rm(j)}
      if (exists("k")) {rm(k)}
      
    }
    
    if (exists("i")) {rm(i)}
    #end of for i loop
    
    
    
    ft_simp <- ft
    
    ft <- ft[ft$hits > 0, ]
    
    ft <- dplyr::select(ft, -ms2P)
    
    ft <- dplyr::mutate(ft, checked = 0, confirmed = NA_character_, process = NA_character_)
    
    ft_formulas <- list()
    
    for (i in seq_len(nrow(ft))) {
      
      idf <- ft$ID[i]
      
      t <- ft_info[[idf]]
      
      t$result <- FALSE
      
      tformula <- obj2@patdata[, idf]
      
      # tformula <- "holder"
      # TODO As is, there is a error for one of the Aopti samples. Task to find it and fix, not using the error catcher
      # TODO talk with Rick about the zeros, check if it can be changed in xdata
      # tryCatch( { tformula <- obj2@patdata[, idf] }
      #     , error = function(e) {
      #         #test headinf zero on feat ID
      #         test <- stringr::str_extract(idf, "(?<=[:punct:])[:digit:]+$")
      #         test <- stringr::str_replace(idf, "(?<=[:punct:])[:digit:]+$", paste0(0,test))
      #         tformula <<- obj2@patdata[, test]})
      
      tformula <- patRoon::generateFormulasSIRIUS(tformula, MS2,
                                                  relMzDev = ppm,
                                                  adduct = adduct,
                                                  elements = "CHNOPSClF",
                                                  profile = "qtof",
                                                  database = NULL,
                                                  noise = NULL,
                                                  topMost = 5,
                                                  absAlignMzDev = 0.0005,
                                                  extraOptsGeneral = NULL,
                                                  calculateFeatures = FALSE)
      
      if (length(tformula) > 0) {
        
        tformula2 <- tformula[[idf]]
        
        
        # transform the table
        for (l in seq_len(nrow(tformula2))) {
          frags <- tformula2$fragInfo[[1]]
          frags$neutral_formula_P <- tformula2$neutral_formula[l]
          frags$ion_formula_P <- tformula2$ion_formula[l]
          frags$mz_P <-  tformula2$ion_formula_mz[l]
          
          if (l == 1) {
            fragments <- frags
          } else {
            fragments <- rbind(fragments,  frags)
          }
          
        } #end of l loop
        
        tformula5 <- list()
        
        for (z in seq_len(nrow(t))) {
          
          t2 <- t[z, ]
          
          if (t2$precursor.x == "Loss") {
            
            # for Loss looks for both fragments annotation
            t4 <- dplyr::inner_join(fragments, t2, by = c("mz" = "mz.x"))
            t4 <- dplyr::select(t4, ion_formula, mz, mz_P, ion_formula_P)
            
            t5 <- dplyr::inner_join(fragments, t2, by = c("mz" = "mz.y"))
            t5 <- dplyr::select(t5, ion_formula, mz, mz_P, ion_formula_P)
            
            t6 <- dplyr::inner_join(t4, t5, by = "ion_formula_P")
            
            if (nrow(t6) > 0) {
            
              vx <- sapply(t6$ion_formula.x, function(x) CHNOSZ::makeup(x))
              
              vy <- sapply(t6$ion_formula.y, function(x) CHNOSZ::makeup(x))
              
              test <- vy[rownames(vy) %in% rownames(vx), , drop = FALSE] - vx[rownames(vx) %in% rownames(vy), , drop = FALSE]
              
              expect <- CHNOSZ::makeup(t2$FD)
              
              t6$loss <- apply(test, MARGIN = 2, function(x, expect) {
                x <- x[x > 0]
                x <- setequal(x, expect)
                return(x)
              }, expect = expect)
              
              t$result[z] <- TRUE %in% t6$loss
              
            }
            
            t6$name <- t2$name
            t6$precursor.x <- t2$precursor.x
            
            tformula5[[z]] <- t6
            
          } else {
            
            
            t3 <- dplyr::inner_join(fragments, t2, by = c("mz" = "mz.x"))
            t3 <- dplyr::mutate(t3, result = FALSE)
            
            if (nrow(t3) > 0) {
              
              for (j in seq_len(nrow(t3))) {
                
                insi <- CHNOSZ::makeup(t3$ion_formula[j])
                
                indb <- CHNOSZ::makeup(t3$FD[j])
                
                #Adds one more H for direct match of protonated fragments
                if (t3$precursor.x[j] == "Direct") indb[names(indb) == "H"] <-  indb[names(indb) == "H"] + 1
                
                if (TRUE == all.equal(insi, indb)) t3$result[j] <- TRUE
                
              } #end of j loop
              
              t$result[z] <- TRUE %in% t3$result
            
            }
            
            tformula5[[z]] <- t3
          }
          
          if (exists("t2")) {rm(t2)}
          if (exists("t3")) {rm(t3)}
          if (exists("t4")) {rm(t4)}
          if (exists("t5")) {rm(t5)}
          if (exists("t6")) {rm(t6)}
          if (exists("insi")) {rm(insi)}
          if (exists("indb")) {rm(indb)}
          if (exists("test")) {rm(test)}
          if (exists("vx")) {rm(vx)}
          if (exists("vy")) {rm(vy)}
          if (exists("expect")) {rm(expect)}
          if (exists("j")) {rm(j)}
          
        } # enf for z for loop
        
        
        ft_formulas[[idf]] <- tformula5
        
        checked <- t[t$result == TRUE, ]
        
        ft$checked[i] <-  nrow(checked)
        ft$confirmed[i] <- paste(as.vector(checked$name), collapse = "; ")
        ft$process[i] <-  paste(as.vector(checked$precursor.x), collapse = "; ")
        
      }
      if (exists("checked")) {rm(checked)}
    } #end of i loop
    if (exists("i")) {rm(i)}
    
    if (r == 1) {
      ft_final <- ft
    } else {
      ft_final <- rbind(ft_final, ft)
    }
  
    extra_data[[rg[r]]] <- list(c(ft_info, ft_formulas))
    
  } #end of r loop
  
  if (nrow(ft_final) == 0) {
    warning("MS2 fragments in targets not found in given features!")
    return(obj)
  }
  
  data@param$replicates <- replicates
  data@param$ID <- ID
  data@param$ppm <- ppm
  data@param$intMin <- intMin
  data@targets <- targets
  
  
  
  data@data <- extra_data
  
  data@targets <- targets
  
  data@results <- ft_final
  
  if (is.null(title)) title <- "MS2FragmentsScreening"

  obj@workflows[[title]] <- data
  
  return(obj)

}
