

### ntsSuspectData -----

#' @title ntsSuspectData
#'
#' @slot suspects The data.frame listing all the suspects of interest.
#' @slot param The parameters used for suspect screening.
#' @slot data The raw \linkS4class{featureGroupsScreening}
#' for each sample replicate group.
#' @slot inSilico When present, the results from the matched In Silico
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
    suspects = "data.frame",
    param = "list",
    data = "list",
    inSilico = "list",
    results = "data.frame"
  ),
  prototype = list(
    suspects = data.frame(),
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



#' @title screenSuspectsFromCSV
#'
#' @description Method to perform suspect screening from a given
#' list of candidates in a csv file. The method uses the S4 method
#' \code{screenSuspects} from \pkg{patRoon}.
#'
#' @param obj An \linkS4class{ntsData} object with features
#' to preform suspect screening.
#' @param name An optional string (without spaces) to be used as identifier for the
#' suspect screening entry in the \code{workflows} slot of the \linkS4class{ntsData} object.
#' @param suspects A \code{data.frame} or a location for a csv
#' with details for each suspect compound.
#' See details for more information about the required data.frame structure.
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
#' @return A data.frame with the suspect screening results.
#'
#' @details The \code{suspects} data.frame should follow the template requires as
#' obtained via \code{getScreeningListTemplate()}.
#' Add other details of the template.
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
#'
#' @examples
#'
screenSuspectsFromCSV <- function(obj,
                                  name = NULL,
                                  suspects = utils::read.csv(base::file.choose()),
                                  ppm = 5,
                                  rtWindow = 30,
                                  adduct = NULL,
                                  excludeBlanks = TRUE,
                                  withMS2 = TRUE,
                                  MS2param = MS2param()) {

  # TODO A ID and sampleGroups argument for running suspect screening only for a set of features.

  assertClass(obj, "ntsData")

  if (!is.data.frame(suspects)) suspects <- utils::read.csv(suspects)

  if (max(suspects$rt) < 120) suspects$rt <- suspects$rt * 60

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

  rg <- unique(sampleGroups(obj))

  if (excludeBlanks) rg <- rg[!(rg %in% blanks(obj))]

  screen <- screenSuspects(obj@patdata, select(suspects, -mz),
                           rtWindow = rtWindow, mzWindow = 0.03,
                           adduct = adduct, onlyHits = TRUE)

  df <- arrange(patRoon::as.data.frame(screen, average = FALSE), group)
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

      temp <- filterFeatureGroups(screen, which(obj@samples$group == rg[g]))

      MS2 <- extractMS2(temp, param = MS2param)

      # TODO Adduct is [M+H]+ by default but should take the value from screening
      formulas <- patRoon::generateFormulasGenForm(temp, MS2,
                                                   relMzDev = ppm,
                                                   isolatePrec = TRUE,
                                                   adduct = "[M+H]+",
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
        tempForm <- formulas@formulas[[isoScore[iso]]]
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

      tempScreen <- cbind(tempScreen, df[, colnames(df) %in% obj@samples$sample])

      data[[rg[g]]] <- temp

      inSilico[[rg[g]]] <- Form

      results[[rg[g]]] <- tempScreen

    }

    results <- data.table::rbindlist(results, idcol = "group")

    results <- base::as.data.frame(results)

  } else {

    data <- screen

    results <- df

  }

  sS <- new("ntsSuspectData")

  sS@suspects <- suspects

  sS@param$ppm <- ppm
  sS@param$rtWindow <- rtWindow
  sS@param$adduct <- adduct
  sS@param$excludeBlanks <- excludeBlanks
  sS@param$withMS2 <- withMS2
  sS@param$MS2param <- MS2param

  sS@data <- data
  sS@inSilico <- inSilico
  sS@results <- results

  if (is.null(name)) name <- "SuspectScreening"

  obj@workflows[[name]] <- sS

  return(obj)

}
