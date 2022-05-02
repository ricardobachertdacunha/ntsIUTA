

### ntsSuspectData -----

#' @title ntsSuspectData
#'
#' @slot suspectList The \linkS4class{suspectList} object used for suspect screening.
#' @slot param The parameters used for suspect screening.
#' @slot data The raw \linkS4class{featureGroupsScreening}
#' for each sample replicate group.
#' @slot inSilico When present, the results from the matched in silico
#' fragmentation experiment.
#' @slot results A data.frame with summarized results, including
#' the identification level according to the guideline.
#'
#' @return An \linkS4class{ntsSuspectData} object to be added to
#' the workflows slot of an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("ntsSuspectData",
  slots = c(
    suspectList = "suspectList",
    param = "list",
    data = "list",
    inSilico = "list",
    results = "data.frame"
  ),
  prototype = list(
    suspectList = new("suspectList"),
    param = list(ppm = 5,
                 rtWindow = 30,
                 adduct = NULL,
                 excludeBlanks = TRUE,
                 withMS2 = TRUE,
                 MS2param = new("MS2param")),
    data = list(),
    inSilico = list(),
    results = data.frame()
  )
)




#' @title suspectScreening
#'
#' @description Method to perform suspect screening from a given
#' list of candidates in a csv file. The method uses the S4 method
#' \code{screenSuspects} from \pkg{patRoon}.
#'
#' @param obj An \linkS4class{ntsData} object with features
#' to preform suspect screening.
#' @param samples A character vector with the name of the
#' samples to perform the suspect screening.
#' @param ID A character vector with the ID of the features to perform suspect screening.
#' @param title An optional string (without spaces) to be used as identifier for the
#' suspect screening entry in the \code{workflows} slot of the \linkS4class{ntsData} object.
#' @param suspectList A \linkS4class{suspectList} object containing suspect compounds.
#' See \link{getSuspectListTemplate} for more details about the information required.
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
#' @param param An \linkS4class{MS2param} object with parameters for MS2 extraction.
#'
#' @return An \linkS4class{ntsData} object with screening results
#' added to the workflows slot.
#'
#' @references \insertRef{Helmus2021}{ntsIUTA}
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom dplyr filter arrange left_join select everything mutate distinct desc
#' @importFrom utils read.csv head
#' @importFrom data.table rbindlist
#' @importMethodsFrom patRoon screenSuspects as.data.frame screenInfo annotateSuspects generateFormulas
#'
suspectScreening <- function(obj,
                             samples = NULL,
                             ID = NULL,
                             title = NULL,
                             suspectList = NULL,
                             ppm = 5,
                             rtWindow = 30,
                             adduct = NULL,
                             excludeBlanks = TRUE,
                             withMS2 = TRUE,
                             param = MS2param()) {

  assertClass(obj, "ntsData")

  assertClass(suspectList, "suspectList")

  if (suspectList@length == 0) {
    warning("The suspectList is empty!")
    return(obj)
  }

  susdf <- suspectList@data

  #selects top 5 or 10 fragment if MS2 data is present for suspects
  if ("hasFragments" %in% colnames(susdf)) {
    for (n in seq_len(nrow(susdf))) {
      if (susdf$hasFragments[n]) {
        top <- unlist(susdf$fragments_mz[n])
        top <- as.data.frame(unlist(strsplit(top, split = ";")))
        if (nrow(top) > 5) {
          colnames(top) <- c("mz")
          top$mz <- as.numeric(top$mz)
          top$int <- as.numeric(unlist(strsplit(susdf$fragments_int[n], split = ";")))
          #remove precursor ion from fragments list
          top <- dplyr::filter(top, mz < susdf$mz[n] - (5 / 1E6 * susdf$mz[n]))
          if (nrow(top) < 15) ntop <- 5 else ntop <- 10
          top <- utils::head(arrange(top, desc(int)), n = ntop)
          susdf$fragments_mz[n] <- paste(top$mz, collapse = ";")
          susdf$fragments_int[n] <- paste(top$int, collapse = ";")
          # TODO add top for formulas but might not be necessary
        }
      }
    }
  }

  if (is.null(adduct)) {
    if (!("adduct" %in% colnames(susdf))) {
      adduct <- ifelse(obj@polarity == "positive", "[M+H]+",
                       ifelse(obj@polarity == "negative", "[M-H]-",
                              stop("polarity argument must be 'positive' or 'negative'")))
    }
  }

  obj2 <- obj


  if (!is.null(samples)) obj2 <- filterFileFaster(obj2, samples)

  rg <- unique(sampleGroups(obj2))
  
  if (excludeBlanks) { # just take the blanks from the replicates for MS2 check
    #if (TRUE %in% (sampleGroups(obj2) %in% blanks(obj2))) {
    #  obj2 <- filterFileFaster(obj2, samples(obj2)[!(sampleGroups(obj2) %in% blanks(obj2))])
    #}
    rg <- rg[!rg %in% blanks(obj2)]
  }
  
  susdf$rt <- as.numeric(susdf$rt)
  
  screen <- screenSuspects(
    fGroups = obj2@patdata, 
    suspects = select(susdf, -mz),
    rtWindow = rtWindow, 
    mzWindow = 0.03,
    adduct = adduct, 
    onlyHits = FALSE
  )
  
  df <- patRoon::as.data.frame(x = screen, average = FALSE, onlyHits = TRUE)
  df <- arrange(df, group)
  df <- rename(df, name = susp_name)
  df <- left_join(df, susdf[, colnames(susdf) %in% c("name", "formula", "adduct", "hasFragments", "intControl")], by = "name")
  df <- left_join(df, select(arrange(patRoon::screenInfo(obj = screen), group), group, d_mz, d_rt), by = "group")
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- dplyr::rename(df, rt = ret, ID = group)
  df <- select(df, name, formula, adduct, ID, mz, rt, everything(), -d_mz, d_ppm, d_rt)

  df <- df[df$ID %in% obj2@features$ID, ] #removed features that were removed during filtering

  df <- dplyr::filter(df, d_ppm <= ppm)
  df$d_ppm <- round(df$d_ppm, digits = 1)
  df$d_rt <- round(df$d_rt, digits = 1)
  
  #TODO Check why Formulas occasionally become NA, breaking Genform with NACHO elements
  #SOLVED The problem is that names appear in duplicate in df. I attempt to fix by spliting when formula is NA
  #TODO reimplement Code !!!!
  df <- df[!(is.na(df$formula)), ]
  # dup <- df[is.na(df$formula),]
  # if (nrow(dup) > 0) {
  #   df <- df[!(is.na(df$formula)), ]
  #   for (i in seq_len(nrow(dup))) {
  #     nDups <- unlist(stringr::str_split(dup[i, 1, drop = TRUE], pattern = ","))
  #     for (nDup in nDups) {
  #       newrow <- dup[i, ]
  #       newrow[1, 1] <- nDup
  #       newrow[1, 2:3] <- susdf[susdf$name %in% nDup, c("formula", "adduct"), drop = TRUE]
  #       df <- rbind(df, newrow)
  #     }
  #   }
  #   df <- dplyr::arrange(df, mz)
  # }
  
  df <- dplyr::distinct(df)
  
  screen <- screen[, df$ID]

  elements <- gsub("[^a-zA-Z]", "", df$formula)
  elements <- paste0(elements, collapse = "")
  elements <- gsub("([[:upper:]])", " \\1", elements)
  elements <- unique(strsplit(elements, " ")[[1]])
  elements <- elements[!elements %in% ""]
  elements <- paste0(elements, collapse = "")

  data <- list()
  inSilico <- list()
  results <- list()
  
  if (withMS2) {

    for (g in seq_len(length(rg))) {

      temp <- screen[which(obj2@samples$group == rg[g])]

      MS2 <- extractMS2(temp, param = param)

      # TODO Adduct is [M+H]+ by default but should take the value from screening list
      formulas <- patRoon::generateFormulasGenForm(temp, MS2,
                                                   relMzDev = ppm,
                                                   isolatePrec = TRUE,
                                                   absAlignMzDev = 0.005,
                                                   adduct = "[M+H]+", # TODO make NULL
                                                   elements = elements,
                                                   topMost = 50, extraOpts = NULL,
                                                   calculateFeatures = FALSE,
                                                   featThreshold = 1,
                                                   timeout = 240, hetero = TRUE, oc = TRUE)

      temp <- patRoon::annotateSuspects(temp, MSPeakLists = MS2,
                                        formulas = formulas,
                                        absMzDev = 0.008,
                                        relMinMSMSIntensity = 0.05,
                                        simMSMSMethod = "cosine",
                                        checkFragments = "mz",
                                        formulasNormalizeScores = "max",
                                        compoundsNormalizeScores = "max")

      tempScreen <- temp@screenInfo[match(df$ID, temp@screenInfo$group)]
      
      # TODO Verify if this solution works in all instances
      if (all(is.na(tempScreen$rt))) {
        tempScreen$rt <- ifelse(df$ID == tempScreen$group, df$rt, NA)
        
      }
      tempScreen <- tempScreen[!is.na(tempScreen$group)]
      

      
      tempScreen$d_ppm <- (abs(tempScreen$d_mz) / tempScreen$mz) * 1E6

      isoScore <- tempScreen$group
      for (iso in seq_len(nrow(tempScreen))) {
        tempForm <- formulas@groupAnnotations[[isoScore[iso]]]
        if (!is.null(tempForm)) {
          tempForm <- tempForm[tempForm$neutral_formula %in% tempScreen$formula, ]
          tempForm$ID <- isoScore[iso]
          isoScore[iso] <- round(unique(tempForm$isoScore)[1], digits = 2)
          if (iso == 1 | !exists("Form")) {
            Form <- tempForm
          } else {
            Form <- rbind(Form, tempForm, fill = TRUE)
          }
        } else {
          isoScore[iso] <- NA
        }
      }

      tempScreen$isoScore <- isoScore

      tempScreen$hasExpFrag <- FALSE
      
      for (i in seq_len(nrow(tempScreen))) {
        fragments <- MS2[[tempScreen$group[i]]]$MSMS
        if (!is.null(fragments)) tempScreen$hasExpFrag[i] <- TRUE
      }

      tempScreen$FragMatch <- paste0(tempScreen$maxFragMatches, "(", tempScreen$maxFrags, ")")

      tempScreen <- dplyr::rename(tempScreen, ID = group, IdLevel = estIDLevel)

      tempScreen$hasFrag <- df$hasFragments[df$ID %in% tempScreen$ID]

      tempScreen <- select(tempScreen, name, formula, adduct, ID, mz, rt, IdLevel, d_rt, d_ppm, isoScore, FragMatch, hasExpFrag, hasFrag, everything())

      tempScreen <- dplyr::left_join(tempScreen, df[, colnames(df) %in% c("ID", obj2@samples$sample)], by = "ID")

      tempScreen <- dplyr::distinct(tempScreen)
      
      data[[rg[g]]] <- temp

      inSilico[[rg[g]]] <- Form

      results[[rg[g]]] <- tempScreen

    }

    results <- data.table::rbindlist(results, idcol = "group", fill = TRUE)

    results <- base::as.data.frame(results)

  } else {

    data <- list(screen)

    results <- df

  }
  
  sS <- new("ntsSuspectData")

  sS@suspectList <- suspectList
  
  sS@param$ppm <- ppm
  sS@param$rtWindow <- rtWindow
  sS@param$adduct <- adduct
  sS@param$excludeBlanks <- excludeBlanks
  sS@param$withMS2 <- withMS2
  sS@param$MS2param <- param
  
  sS@data <- data
  sS@inSilico <- inSilico
  sS@results <- results
  
  if (is.null(title)) title <- "SuspectScreening"

  obj@workflows[[title]] <- sS

  return(obj)

}
