

findFragments <- function(obj, targets = NULL, replicates = NULL, ID = NULL,
                          MS2param = NULL) {
  
  assertClass(obj, "ntsData")
  
  obj2 <- obj
  
  adduct <- ifelse(obj@polarity == "positive", "[M+H]+", ifelse(obj@polarity == "negative", "[M-H]-", NULL))
  if (is.null(adduct)) {
    warning("Polarity of ntsData not recognized!")
    return(obj)
  }
  
  
  ### Collect Features -----
  
  #filter files for selected sampleGroups
  if (!is.null(replicates)) obj2 <- filterFileFaster(obj2, which(sampleGroups(obj2) %in% replicates))
  
  #select features of interest
  if (!is.null(ID)) obj2 <- obj2[, obj2@features$ID[obj2@features$ID %in% ID]]
  
  ft <- obj2@features
  
  
  
  ### Search parameters -----
  
  ppm <- 10
  
  intMin <- 10
 
  
  
  ### Load Targets -----
  
  targets_frag <- dplyr::filter(targets, process == "Fragment")
  targets_frag <- dplyr::mutate(targets_frag, MD = -MD)
  
  targets_cleav <- dplyr::filter(targets, process == "Cleavage")
  
  
  
  ### Load MS2 -----
  
  if (is.null(MS2param)) MS2param <- obj2@parameters@MS2
  
  MS2 <- extractMS2(obj2@patdata, param = MS2param)
  
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
          ft_info[[idf]] <- tTP
          
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
  
  ft <- dplyr::select(ft, -ms2P, -hasFragments)
  
  ft <- dplyr::mutate(ft, checked = 0, confirmed = NA_character_, process = NA_character_)
  
  ft_formulas <- list()
  
  for (i in seq_len(nrow(ft))) {
    
    idf <- ft$ID[i]
    
    t <- ft_info[[idf]]
    
    t$result <- FALSE
    
    tformula <- obj2@patdata[, idf]
    
    tformula <- patRoon::generateFormulasSIRIUS(tformula, MS2,
                                                relMzDev = 5,
                                                adduct = adduct,
                                                elements = "CHNOPSClF",
                                                profile = "qtof",
                                                database = NULL,
                                                noise = NULL,
                                                topMost = 5,
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
            
            test <- vy
            test[rownames(test) %in% rownames(vx), ] <- test[rownames(test) %in% rownames(vx), ] - vx
            
            expect <- CHNOSZ::makeup(t2$FD)
            
            t6$loss <- apply(test, MARGIN = 2, function(x, expect) {
              x <- x[x > 0]
              x <- all.equal(x, expect)
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
      
      checked <- lapply(tformula5, function(x) unique(dplyr::select(x, name, precursor.x)))
      checked <- bind_rows(checked, .id = "column_label")
      checked <- dplyr::select(checked, name, precursor.x)
      
      ft$checked[i] <-  nrow(t[t$result == TRUE, ])
      ft$process[i] <-  paste(as.vector(tformula6$precursor.x), collapse = "; ")
      Nitro_MS2$confirmed[i] <- paste(as.vector(tformula6$name), collapse = "; ")
      Nitro_MS2$MS2.x[i] <- paste(as.vector(tformula6$mz.x), collapse = "; ")
      Nitro_MS2$MS2.y[i] <- paste(as.vector(tformula6$mz.y), collapse = "; ")
    }
    
    
    
    
    
    
    
  } #end of i loop
  if (exists("i")) {rm(i)}
  
  
  
  
  Nitro_MS2_Groups <- dplyr::select(Nitro_MS2, group)
  
  nitro <- list()
  
  nitro[["nitro_details"]] <- Nitro_MS2_Details
  nitro[["nitro_formulas"]] <- Nitro_MS2_Formulas
  nitro[["nitro_list"]] <- Nitro_MS2

  if (do.monitoring == TRUE) {
    
    # Intensity Filter
    Monitor_Nitro_IntThres <- 3000
    
    Nitro_DB <- read.csv("Database_TP_Nitro.csv")
    Nitro_DB_WTP <- dplyr::filter(Nitro_DB, process == "Water TP")
    Nitro_DB_WTP <- dplyr::select(Nitro_DB_WTP, -process)
    
    Monitor_Nitro_T <- list()
    
    for ( i in 1:length(monitoring[["monitor_process"]])) {
      
      n1 <- dplyr::filter(monitoring[["monitor_process"]][[i]], cat == "N" | cat == "H")
      n1 <- dplyr::filter(n1, n1[, which(colnames(monitoring[["monitor_process"]][[i]]) == monitoring$monitor_count$sample[i+1])] > Monitor_Nitro_IntThres)
      #n1 <- dplyr::filter(n1, n1[, 4+i] > Monitor_Nitro_IntThres)
      n1 <- dplyr::select(n1, group, ret, mz, qtl, cat, CMP)
      n1 <- as.data.frame(n1)
      
      n2 <- dplyr::filter(monitoring[["monitor_process"]][[i]], cat == "R" | cat == "L")
      n2 <- dplyr::filter(n2, n2[, which(colnames(monitoring[["monitor_process"]][[i]]) == monitoring$monitor_count$sample[i])] > Monitor_Nitro_IntThres)
      #n2 <- dplyr::filter(n2, n2[, 3+i] > Monitor_Nitro_IntThres)
      n2 <- dplyr::select(n2, group, ret, mz, qtl, cat, CMP)
      n2 <- as.data.frame(n2)
      
      n2 <- merge(n2, Nitro_DB_WTP, by = NULL, all = TRUE)
      n2 <- dplyr::mutate(n2, mzTP = mz + MD)
      n1 <- fuzzyjoin::difference_inner_join(n1, n2, by = c("mz" = "mzTP"), max_dist = 0.005, distance_col = "mzdiff")
      Monitor_Nitro_T[[i]] <- n1
      rm(n1, n2)
    }
    
    
    Monitor_Nitro <- list()
    for ( i in 1:length(monitoring[["monitor_process"]])) {
      
      Monitor_Nitro[[i]] <-  Monitor_Nitro_T[[i]]
      Monitor_Nitro[[i]] <- inner_join(Monitor_Nitro[[i]],  nitro[["nitro_list"]], by = c("group.y" = "group"))
      Monitor_Nitro[[i]] <- dplyr::mutate(Monitor_Nitro[[i]], Nitro = NA, NitroCheck = NA)
      for (j in 1:nrow(Monitor_Nitro[[i]])) { #nrow(Monitor_Nitro[[i]])
        temp1 <- Monitor_Nitro[[i]]
        temp2 <- as.vector(temp1$name[j])
        temp3 <- temp1$names[j]
        temp3 <- unlist(strsplit(temp3, split=";"))
        temp4 <- temp1$confirmed[j]
        temp4 <- unlist(strsplit(temp4, split=";"))
        if (temp2 %in% temp3) {
          Monitor_Nitro[[i]]$Nitro[j] <- temp2
        } else {
          Monitor_Nitro[[i]]$Nitro[j] <- NA
        }
        if (temp2 %in% temp4) {
          Monitor_Nitro[[i]]$NitroCheck[j] <- TRUE
        } else {
          Monitor_Nitro[[i]]$NitroCheck[j] <- FALSE
        }
      }
    }
    if (exists("temp1")) {rm(temp1)}
    if (exists("temp2")) {rm(temp2)}
    if (exists("temp3")) {rm(temp3)}
    if (exists("temp4")) {rm(temp4)}
    if (exists("i")) {rm(i)}
    if (exists("j")) {rm(j)}
    
    nitro[["monitor_nitro"]] <- Monitor_Nitro
  }

return(nitro)
}


