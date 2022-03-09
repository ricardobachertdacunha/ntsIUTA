

#' peakAnnotation
#'
#' @description Annotation of isotopic and adduct peaks.
#'
#' @param object An \linkS4class{ntsData} object with peaks/features.
#' @param algorithm The algorithm for annotation.
#' @param settings The respective parameter settings.
#' @param excludeBlanks Logical, set to \code{TRUE} for excluding
#' blank replicates from annotation. Blank intensities are not considered.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @return An \linkS4class{ntsData} object with annotated peaks/features.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importMethodsFrom patRoon generateComponents
#' @importClassesFrom patRoon components
#'
peakAnnotation <- function(object = NULL,
                           algorithm = NULL,
                           settings = NULL,
                           save = FALSE) {

  checkmate::assertClass(object, "ntsData")

  pat <- object@pat

  if (is.null(algorithm)) algorithm <- annotationParameters(object)@algorithm

  if (is.null(settings)) settings <- annotationParameters(object)@settings

  if (is.na(algorithm)) {
    warning("Annotation algorihtm not defined!")
    return(object)
  }

  ag <- list(fGroups = pat, algorithm = algorithm)

  comp <- do.call(patRoon::generateComponents, c(ag, settings))

  object@comp <- comp

  if (polarity(object) %in% "positive") {
    prefAdduct <- "[M+H]+"
  } else {
    prefAdduct <- "[M-H]-"
  }

  object <- annotateFeaturesTable(object, prefAdduct)

  if (is.logical(save)) if (save) saveObject(object = object)

  if (is.character(save)) saveObject(object = object, filename = save)

  return(object)
}




#' @title annotateFeaturesTable
#'
#' @description Selects the ions from the components and amends the festure data table.
#'
#' @param object An \linkS4class{ntsData} object with peaks/features.
#' @param prefAdduct A character string with the preferred adduct.
#' Possible values are \code{"[M+H]+"} (the default) and \code{"[M-H]-"} for
#' positive and negative ionization, respectively.
#'
#' @importMethodsFrom patRoon componentTable
#' @importFrom data.table rbindlist setnames fread
#'
annotateFeaturesTable <- function(object, prefAdduct = "[M+H]+") {

  comp <- object@comp
  feats <- object@features

  cat("Annotating features based on components... \n")

  pb <- txtProgressBar(
    min = 0,
    max = nrow(aft),
    style = 3,
    width = 50,
    char = "="
  )

  #change to ntsIUTA package folder
  if ("[M+H]+" %in% prefAdduct) {
    db_adducts <- fread(system.file("rules/primary_adducts_pos.csv", package = "CAMERA"), header = TRUE)
  } else {
    db_adducts <- fread(system.file("rules/primary_adducts_neg.csv", package = "CAMERA"), header = TRUE)
  }

  comp_df <- componentTable(comp)
  comp_df <- rbindlist(comp_df, idcol = "component")

  #In cliqueMS colnames are: neutralMass (changes to M_adduct) isonr charge adduct_ion  intensity intensity_rel
  #In ramclustr colnames are: intensity intensity_rel group isogroup isonr charge adduct_ion ppm
  setnames(comp_df, "neutralMass", "M_adduct", skip_absent = TRUE)

  if (!"M_adduct" %in% colnames(comp_df)) {
    comp_df[, M_adduct := as.numeric(NA)]
  }

  if (!"isogroup" %in% colnames(comp_df)) {
    comp_df[, isogroup := as.numeric(stringr::str_extract(comp_df$component, "[:digit:]"))]
  }
  comp_df$isonr <- as.numeric(comp_df$isonr)

  aft <- feats[, .(id, mz, rt)]
  aft <- aft[, `:=`(
    "neutralMass" = round(mz - db_adducts[name %in% prefAdduct, massdiff], digits = 4),
    "isonr" = 0,
    "monoiso" = id,
    "isogroup" = NA_character_,
    "charge" = 1,
    "adduct_ion" = prefAdduct,
    "intensity" = NA,
    "rel_intensity" = NA,
    "component" = NA_character_,
    "annotated" = NA_integer_
  )]

  aft$rel_intensity <- as.numeric(aft$rel_intensity)
  aft$intensity <- as.numeric(aft$intensity)

  #priority for [M+H]+ or [M-H]-
  #smallest rt difference from correspondent prefAdduct
  #smallest mz differnece from correspondence prefAdduct
  for (f in seq_len(nrow(aft))) {

    x <- aft$id[f]

    temp <- comp_df[group %in% x, ]

    aft[id %in% x, annotated := nrow(temp)]

    #filter preferential adduct, isotopes are not likely to be duplicated
    if (nrow(temp) > 1) {
      if (prefAdduct %in% temp$adduct_ion) {
        temp <- temp[adduct_ion %in% prefAdduct, ]
      } else {
        temp[, `:=`(rt_d = 0, mz_d = 0)]
        for (i in seq_len(nrow(temp))) {
          temp2 <- comp_df[component %in% temp$component[i] & M_adduct == temp$M_adduct[i], ]

          if (prefAdduct %in% temp2$adduct_ion) {
            temp$rt_d[i] <-  abs(temp$ret[i] - temp2[adduct_ion %in% prefAdduct, ret])
            temp$mz_d[i] <-  abs(
              temp$mz[i] -
              temp2[adduct_ion %in% prefAdduct, mz] -
              db_adducts[name %in% temp$adduct_ion[i], massdiff] +
              db_adducts[1, massdiff]
            )
          } else if (length(unique(temp2$adduct_ion)) == 1) {

            temp[1, colnames(temp)[stringr::str_detect(colnames(temp), "ad")] := NA]
            temp <- temp[1, ]
            break

          } else {
            temp$rt_d[i] <-  abs(temp$ret[i] - temp2[M_adduct %in% temp$M_adduct[i] & !adduct_ion %in% temp$adduct_ion[i], ret])
            temp$mz_d[i] <-  abs(
              temp$mz[i] -
              temp2[M_adduct %in% temp$M_adduct[i] & !adduct_ion %in% temp$adduct_ion[i], mz] -
              db_adducts[name %in% temp$adduct_ion[i], massdiff] +
              db_adducts[name %in% temp2[M_adduct %in% temp$M_adduct[i] & !adduct_ion %in% temp$adduct_ion[i], adduct_ion], massdiff]
            )
          }
        }
        temp <- temp[rt_d == min(temp$rt_d), ] #lowest rt diff
        temp <- temp[mz_d == min(temp$mz_d), ] # lowest mz diff
      }
    }

    if (nrow(temp) == 1) {

      ntm <- temp$M_adduct

      #amend neutralMass for multiple charged isotopes
      if (TRUE %in% (temp$charge > 1) & TRUE %in% (temp$isonr == 0)) {
        ntm <- (aft[id %in% x, mz] - db_adducts[name %in% prefAdduct, massdiff]) * temp$charge
      }

      #amend direct ions to the protonted ion
      if (TRUE %in% grepl("[M]", temp$adduct_ion, fixed = TRUE)) {
        #retrieve the mass of the respective adduct
        mono_id <- comp_df[component %in% temp$component & isogroup %in% temp$isogroup & isonr %in% 0, ]

        if (nrow(mono_id) > 1) {
          mono_id[, dist := abs(temp$mz - mz)]
          mono_id <- mono_id[dist == min(dist), ]
        }

        mono_id <- mono_id[, group]

        if (length(mono_id) > 0) {
          ntm <- aft[id %in% mono_id, neutralMass]
          aft[id %in% x, isonr := temp$isonr]
          aft[id %in% x, monoiso := mono_id]
          aft[id %in% x, rel_intensity := as.numeric(temp$intensity) / comp_df[group %in% mono_id, intensity]]
        }
        mono_id <- character()
      }

      #amend isotopic peaks to related monoisotope
      if (TRUE %in% (temp$isonr > 0)) {

        #retrieve the mass of the respective adduct
        mono_id <- comp_df[component %in% temp$component & isogroup %in% temp$isogroup & isonr %in% 0, ]

        if (nrow(mono_id) > 1) {
          mono_id[, dist := abs(temp$mz - mz)]
          mono_id <- mono_id[dist == min(dist), ]
        }

        mono_id <- unique(mono_id[, group])

        if (length(mono_id) > 0) {
          ntm <- aft[id %in% mono_id, neutralMass]
          aft[id %in% x, isonr := temp$isonr]
          aft[id %in% x, monoiso := mono_id]
          temp$adduct_ion <- paste0("[M+", temp$isonr, "]")
          aft[id %in% x, rel_intensity := as.numeric(temp$intensity) / aft[id %in% mono_id, intensity]]
        }
        mono_id <- character()
      }

      if (!is.na(ntm)) {
        aft[id %in% x, neutralMass := ntm]
        aft[id %in% x, adduct_ion := temp$adduct_ion]
      }
      ntm <- NA

      aft <- aft[id %in% x, `:=`(
        charge = temp$charge,
        intensity = as.numeric(temp$intensity),
        component = temp$component,
        isogroup = temp$isogroup
      )]
    }

    setTxtProgressBar(pb, f)
  }

  close(pb)
  cat("Done!")

  aft[is.na(rel_intensity), rel_intensity := 1]
  aft$rel_intensity <- round(aft$rel_intensity, digits = 4)
  aft[, c("mz", "rt") := NULL]

  colN <- colnames(aft)
  colN <- colN[!colN %in% "id"]
  keepColN <- colnames(feats)[!colnames(feats) %in% colN]

  feats <- dplyr::left_join(feats[, keepColN, with = FALSE], aft, by = "id")

  object@features <- feats

  return(object)
}

#length(feats$id) - length(unique(comp_df$group))

#plotFeatures(object, targets = c("M100_R585_1", "M105_R585_4"), interactive = TRUE)

  # View(object@pat@features@features[[1]])

  # x <- "M239_R936_638"
  # x <- "M233_R646_572"
  # x <- "M125_R587_36"
  # x <- "M435_R584_2383"
  # x <- "M516_R712_2924"
  # x <- "M244_R936_684"

  # test2 <- comp_df[duplicated(comp_df$group), ]

  # tar <- aft[id %in% c("M239_R936_638"), ]
  # tar <- aft[neutralMass == tar$neutralMass, ]
  # tar

#  View(tar)

  # feats[id %in% "M109_R612_12"]

  # test3 <- comp_df[component %in% "CMP27", ]
  # test3 <- comp_df[isogroup %in% "65", ]

  # comp_df[component == unique(tar$component), ]

  # head(comp_df)

  # length(feats$id) == length(unique(comp_df$group))

  # head(aft)

  # nrow(comp_df)
  # nrow(feats)

  # test <- selectIons(fGroups = pat, components = comp, prefAdduct = "[M+H]+", onlyMonoIso = FALSE, chargeMismatch = "isotope")
  # test@annotations[group %in% c("M239_R936_638", "M240_R936_647"), ]

# test <- comp_df[charge > 1 & !is.na(adduct_ion), ]

# test2 <- comp_df[isogroup == 470, ]

# test3 <- comp_df[component %in% "CMP2", ]

# #prepare data frame
#     for (r in seq_len(length(xAL))) {
#       df_temp <- CAMERA::getPeaklist(xAL[[r]], intval = "maxo")
#       df_temp$group <- names(xAL)[r]
      
#       #for one sample, peaks might be missing
#       sp <- samples(obj)[sampleGroups(obj) == names(xAL)[r]]
      
#       if (length(sp) > 1) {
        
#         df_temp$ID <- obj@features$ID
#         df_temp$intensity <- apply(df_temp[, gsub("[-.]", "", colnames(df_temp)) %in% gsub("[-.]", "", samples(obj))], MARGIN = 1, function(x) mean(x, na.rm = TRUE))
#         df_temp$intensity_sd <- apply(df_temp[, gsub("[-.]", "", colnames(df_temp)) %in% gsub("[-.]", "", samples(obj))], MARGIN = 1, function(x) sd(x, na.rm = TRUE))
#         df_temp$intensity[is.na(df_temp$intensity)] <- 0
#         df_temp$intensity_sd[is.na(df_temp$intensity_sd)] <- 0
#         df_temp <- df_temp[, !gsub("[-.]", "", colnames(df_temp)) %in% gsub("[-.]", "", samples(obj))]
#         df_temp <- df_temp[, !colnames(df_temp) %in% rg]
        
#       } else {
        
#         df_temp$ID <- obj@features$ID[obj@features[, names(xAL)[r], drop = TRUE] > 0]
#         df_temp <- rename(df_temp, intensity = maxo)
#         df_temp$intensity_sd <- 0
#         df_temp$intensity[is.na(df_temp$intensity)] <- 0
#         df_temp <- select(df_temp, -into, -sample)
#         df_temp$npeaks <- 1
#         df_temp <- select(df_temp, mz, mzmin, mzmax, rt, rtmin, rtmax, npeaks, isotopes, adduct, pcgroup, group,  ID, intensity, intensity_sd)

#       }
      
#       if (r == 1 | length(xAL) == 1) {
#         df <- df_temp
#       } else {
#         df <- rbind(df, df_temp)
#       }

#     }
#     df$pcgroup <- paste0(df$pcgroup, "_", df$group)
#     df <- rename(df, comp = pcgroup)


  
#     # TODO Addapt to other algorithms
#     ag <- list(fGroups = obj@patdata, algorithm = algorithm)
#     pat <- do.call(generateComponents, c(ag, param))



#   ## Isotopic information for each feature per sample replicate group
#   isotopes <- df$isotopes
#   isotopes <- strsplit(isotopes, split = "\\]\\[")

#   isogroup <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][1])
#   isogroup <- str_extract(isogroup, pattern = "([0-9]+)")
#   noIsogroup <- is.na(isogroup)
#   isogroup <- paste0(isogroup, "_", df$group)
#   isogroup[noIsogroup] <- NA

#   isoclass <- lapply(X = seq_len(length(isotopes)), function(x) isotopes[[x]][2])
#   isoclass <- as.vector(ifelse(!is.na(isoclass), paste0("[",isoclass), NA))

#   isonr <- gsub("\\[M\\]\\+|\\[M\\]\\-", 0, isoclass)
#   isonr <- as.numeric(str_extract(isonr, "[0-9]"))

#   df$isogroup <- isogroup
#   df$isoclass <- isoclass
#   df$isonr <- isonr


#   ## Make Mions
#   # When not annotated, charge is 1 by default and respective Mion is calculated
#   Mion <- df$mz
#   index <- seq_len(length(Mion))
#   isMion <- grepl(pattern = "[M]", isoclass, fixed = TRUE)

#   charge <- sapply(index, FUN =  function(x) {
#    ifelse(isMion[x], as.numeric(str_extract(isoclass[x], pattern = "([0-9]+)")), NA)
#   })

#   charge[is.na(charge)] <- 1
#   if (obj@polarity == "positive") {
#     baseadduct <- "[M+H]+"
#     Mion <- (Mion - 1.007276) * charge
#   } else {
#     baseadduct <- "[M-H]-"
#     Mion <- (Mion + 1.007276) * charge
#   }
#   Mion <- round(Mion, digits = 3)


#   ## Collect possible adducts and add conflicts information
#   adducts <- strsplit(df$adduct, split = " ")
#   conflicts <- unlist(lapply(X = seq_len(length(adducts)), function(x) ifelse(length(adducts[[x]]) > 2, TRUE, FALSE)))

#   adductclass <- lapply(X = seq_len(length(adducts)), function(x) grep("M", adducts[[x]], value = TRUE))
#   adductMion <- lapply(X = seq_len(length(adducts)), function(x) as.numeric(grep("M", adducts[[x]], value = TRUE, invert = TRUE)))

#   possibleAdducts <- unlist(lapply(adductMion, length))
#   nAdducts <- possibleAdducts
#   isAdduct <- possibleAdducts != 0


#   ## Updates Mion with non-conflicting adducts
#   Mion <- unlist(lapply(index, FUN = function(x) ifelse(isAdduct[x] & !conflicts[x], as.numeric(adductMion[[x]]), Mion[x])))

#   df$conflicts <- conflicts
#   df$nAdducts <- nAdducts

#   adductclass <- lapply(adductclass, function(x) {
#     y <- x
#     if (length(x) == 0) y <- NA
#     return(y)
#   })

#   adductMion <- lapply(adductMion, function(x) {
#     y <- x
#     if (length(x) == 0) y <- NA
#     return(y)
#   })

#   df$adductMion <- I(adductMion)

#   df$adductclass <- I(adductclass)


#   ## solve conflicts
#   # Rules
#   # - Give priority to [M+H]+
#   # - When ties between adducts return the one with lowest mass error

#   consolidate <- lapply(index, FUN =  function(x, dt, adductclass) {

#     hasH <- numeric()

#     hasP <- numeric()

#     y <- adductclass[[x]]

#     if (length(y) > 1) {

#       hold <- dt[x, ]

#       #priority
#       hasH <- grep("H]", y, fixed = TRUE)

#       if (length(hasH) > 0) {

#         df$adductMion[x] <<- unlist(hold$adductMion)[hasH]
#         df$adductclass[x] <<- y[hasH]
#         Mion[x] <<- unlist(hold$adductMion)[hasH]

#       } else {

#         #when priority adduct not present, calculate mass error and return the lowest
#         massrules <- rules[rules$name %in% y, ]
#         massrules <- massrules[match(y, massrules$name), ]
#         errors <- data.frame(mions = unlist(hold$adductMion),
#                              adu = y,
#                              mz = rep(hold$mz, length(y)),
#                              massdiff = massrules$massdiff)
#         errors$expect <- errors$mz - errors$massdiff
#         errors$error <- abs(errors$mions - errors$expect)

#         hasP <- min(errors$error, na.rm = TRUE)
#         hasP <- which(errors$error == hasP)

#         if (length(hasP) < 2) {
#           df$adductMion[x] <<- unlist(hold$adductMion)[hasP]
#           df$adductclass[x] <<- y[hasP]
#           Mion[x] <<- unlist(hold$adductMion)[hasP]
#         } else {
#           #tied case for same mass error, leave it NA with default Mion
#           df$adductMion[x] <<- NA
#           df$adductclass[x] <<- NA
#         }
#       }
#     }
#   }, adductclass = adductclass, dt = df)


#   ## Update Mion with isotopologues after settling adducts priority when multiple hits
#   isopolog <- !isMion
#   isopolog[is.na(isoclass)] <- FALSE

#   names(Mion) <- isogroup

#   Mion <- unlist(lapply(index, function(x) ifelse(isopolog[x], Mion[names(Mion) %in% isogroup[x]], Mion[x])))

#   df$Mion <- as.numeric(Mion)


#   ## finishing data.frame and adding info to obj
#   df <- select(df, ID, group, comp, Mion, isonr, isoclass, isogroup, adductclass, adductMion, conflicts, nAdducts, everything())

#   obj@annotation$comp <- df

#   ## TODO make componentsFeatures to integrate in patRoon workflow

#   obj@annotation$raw <- xAL

#   obj <- annotationParameters(obj, algorithm = algorithm, param = param)
  
#   if (is.character(save)) {
#     saveObject(obj = obj, filename = save)
#   } else {
#     if (save) saveObject(obj = obj)  
#   }
