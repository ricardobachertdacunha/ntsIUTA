

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
#' @param MS2param An \linkS4class{MS2param} object with parameters for MS2 extraction.
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
                             suspectList = utils::read.csv(base::file.choose()),
                             ppm = 5,
                             rtWindow = 30,
                             adduct = NULL,
                             excludeBlanks = TRUE,
                             withMS2 = TRUE,
                             MS2param = MS2param()) {

  assertClass(obj, "ntsData")

  assertClass(suspectList, "suspectList")

  if (suspectList@rtUnit != "sec") {
    warning("The rt should be in seconds!")
    return(obj)
  }

  if (suspectList@length == 0) {
    warning("The suspectList is empty!")
    return(obj)
  }

  suspects <- suspectList@data

  #selects top 5 or 10 fragment if MS2 data is present for suspects
  if ("hasFragments" %in% colnames(suspects)) {
    for (n in seq_len(nrow(suspects))) {
      if (suspects$hasFragments[n]) {
        top <- unlist(suspects$fragments_mz[n])
        top <- as.data.frame(unlist(strsplit(top, split = ";")))
        if (nrow(top) > 5) {
          colnames(top) <- c("mz")
          top$mz <- as.numeric(top$mz)
          top$int <- as.numeric(unlist(strsplit(suspects$fragments_int[n], split = ";")))
          #remove precursor ion from fragments list
          top <- dplyr::filter(top, mz < suspects$mz[n] - (5 / 1E6 * suspects$mz[n]))
          if (nrow(top) < 15) ntop <- 5 else ntop <- 10
          top <- utils::head(arrange(top, desc(int)), n = ntop)
          suspects$fragments_mz[n] <- paste(top$mz, collapse = ";")
          suspects$fragments_int[n] <- paste(top$int, collapse = ";")
          # TODO add top for formulas but might not be necessary
        }
      }
    }
  }

  if (is.null(adduct)) {
    if (!("adduct" %in% colnames(suspects))) {
      adduct <- ifelse(obj@polarity == "positive", "[M+H]+",
                       ifelse(obj@polarity == "negative", "[M-H]-",
                              stop("polarity argument must be 'positive' or 'negative'")))
    }
  }

  obj2 <- obj


  if (!is.null(samples)) obj2 <- filterFileFaster(obj2, samples)

  if (excludeBlanks) {
    if (TRUE %in% (sampleGroups(obj2) %in% blanks(obj2))) {
      obj2 <- filterFileFaster(obj2, samples(obj2)[!(sampleGroups(obj2) %in% blanks(obj2))])
    }
  }

  rg <- unique(sampleGroups(obj2))

  if (!is.null(ID)) {
    obj2@patdata <- obj2@patdata[, ID]
    obj2@features <- obj2@features[obj2@features$ID %in% ID, ]
  }

  screen <- screenSuspects(obj2@patdata, select(suspects, -mz),
                           rtWindow = rtWindow, mzWindow = 0.03,
                           adduct = adduct, onlyHits = TRUE)

  df <- arrange(patRoon::as.data.frame(screen, average = FALSE), group)
  df <- rename(df, name = susp_name)
  df <- left_join(df, suspects[, colnames(suspects) %in% c("name", "formula", "adduct", "hasFragments", "intControl")], by = "name")
  df <- left_join(df, select(arrange(screenInfo(screen), group), group, d_mz, d_rt), by = "group")
  df$d_ppm <- (abs(df$d_mz) / df$mz) * 1E6
  df <- dplyr::rename(df, rt = ret, ID = group)
  df <- select(df, name, formula, adduct, ID, mz, rt, everything(), -d_mz, d_ppm, d_rt)

  df <- dplyr::filter(df, d_ppm <= ppm)
  df$d_ppm <- round(df$d_ppm, digits = 1)
  df$d_rt <- round(df$d_rt, digits = 1)

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

      MS2 <- extractMS2(temp, param = MS2param)

      # TODO Adduct is [M+H]+ by default but should take the value from screening list
      formulas <- patRoon::generateFormulasGenForm(temp, MS2,
                                                   relMzDev = ppm,
                                                   isolatePrec = TRUE,
                                                   adduct = "[M+H]+", # TODO make NULL
                                                   elements = elements,
                                                   topMost = 25, extraOpts = NULL,
                                                   calculateFeatures = FALSE,
                                                   featThreshold = 1,
                                                   timeout = 240, hetero = TRUE, oc = TRUE)

      temp <- patRoon::annotateSuspects(temp, MSPeakLists = MS2,
                                        formulas = formulas,
                                        compounds = NULL)

      tempScreen <- temp@screenInfo[match(df$ID, temp@screenInfo$group)]

      tempScreen$d_ppm <- (abs(tempScreen$d_mz) / tempScreen$mz) * 1E6

      isoScore <- tempScreen$group
      for (iso in seq_len(nrow(tempScreen))) {
        tempForm <- formulas@groupAnnotations[[isoScore[iso]]]
        tempForm <- tempForm[tempForm$neutral_formula %in% tempScreen$formula, ]
        tempForm$ID <- isoScore[iso]
        isoScore[iso] <- round(unique(tempForm$isoScore)[1], digits = 2)
        if (iso == 1) {
          Form <- tempForm
        } else {
          Form <- rbind(Form, tempForm, fill = TRUE)
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

      tempScreen$hasFrag <- df$hasFragments

      tempScreen <- select(tempScreen, name, formula, adduct, ID, mz, rt, IdLevel, d_rt, d_ppm, isoScore, FragMatch, hasExpFrag, hasFrag, everything())

      tempScreen <- cbind(tempScreen, df[, colnames(df) %in% obj2@samples$sample])

      data[[rg[g]]] <- temp

      inSilico[[rg[g]]] <- Form

      results[[rg[g]]] <- tempScreen

    }

    results <- data.table::rbindlist(results, idcol = "group")

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
  sS@param$MS2param <- MS2param

  sS@data <- data
  sS@inSilico <- inSilico
  sS@results <- results

  if (is.null(title)) title <- "SuspectScreening"

  obj@workflows[[title]] <- sS

  return(obj)

}
