

#' getSuspectList
#'
#' @param path The complete path to load the csv file of the list of suspects.
#' Note that the \emph{suspectListTemplate} from \pkg{ntsIUTA}
#' should be used to build the suspect list without changing the columns.
#' The template can be obtained with \link{getSuspectListTemplate}.
#' @param comment Optional comment for the suspect list.
#'
#' @return An \linkS4class{suspectList} object containing information of
#' suspect compounds used for screening workflows within \pkg{ntsIUTA}.
#'
#' @importFrom utils read.csv
#'
#' @export
#'
getSuspectList <- function(path = NULL, comment = NA_character_) {

  sl <- read.csv(path)

  if (FALSE %in% (colnames(sl) %in%
                  c("name",
                    "formula",
                    "neutralMass",
                    "SMILES",
                    "InChIKey",
                    "adduct",
                    "mz",
                    "rt",
                    "intControl",
                    "logP",
                    "comment",
                    "hasFragments",
                    "fragments_mz",
                    "fragments_int",
                    "fragments_formula",
                    "fragments_pre"))) {

    stop("Columns do not match the suspectListTemplate!
         See ?getSuspectListTemplate for more information.")

  }

  sl[is.na(sl)] <- ""

  if (max(sl$rt) != 0 & max(sl$rt) != "") {
    if (max(sl$rt) < 60) {
      warning("It seems that the rt is in minutes.
              It is recommended to change for seconds!")
    }
  }

  sl2 <- new("suspectList",
             path = path,
             data = sl,
             length = nrow(sl),
             comment = comment)

  return(sl2)

}




#' getSuspectListTemplate
#'
#' @param saveTo The path to save the \emph{suspectListTemplate.csv} file.
#'
#' @return A csv template for building a \linkS4class{suspectList} object.
#'
#' @importFrom utils read.csv write.csv
#'
#' @export
#'
getSuspectListTemplate <- function(saveTo = NULL) {

  template <- data.frame(name = "Carbamazepine",
                         formula = "C15H12N2O",
                         neutralMass = 236.094955,
                         SMILES = "c1ccc2c(c1)C=Cc3ccccc3N2C(=N)O",
                         InChIKey = "FFGPTBGBLSHEPO-UHFFFAOYSA-N",
                         adduct = "[M+H]+",
                         mz = 237.102231,
                         rt = "",
                         intControl = "",
                         logP = 0.79,
                         comment = "",
                         hasFragments = FALSE,
                         fragments_mz = "",
                         fragments_int = "",
                         fragments_formula = "",
                         fragments_pre = "")

  write.csv(template, paste0(saveTo, "/suspectListTemplate.csv"),
            row.names = FALSE)

  #TODO add details for the suspectList in the description as details

}




#' addMS2Info
#'
#' @param wfobj An \linkS4class{ntsSuspectData} object
#' @param sampleGroup A character object with the name of a
#' representative sample replicate group to extract the MS2 information.
#' @param suspectList The \linkS4class{suspectList} object to add the MS2 information.
#' @param MS2param The parameters to extract the MS2 from the \linkS4class{ntsSuspectData} object.
#' @param updateRT Logical, set to \code{TRUE} to update or add the retention time
#' to the \linkS4class{suspectList}.
#' @param updateIntControl Logical, set to \code{TRUE} to add the intensity
#' of each compound to the \emph{intControl} column in the \linkS4class{suspectList}.
#' @param save Logical, set to \code{TRUE} to update the csv file in the disk.
#' Optionally, a character object with the new name for the csv file.
#' Note, do not add the file extension to the name.
#'
#' @return A \linkS4class{suspectList} object with updated MS2 information
#' for each compound identified in the \linkS4class{ntsSuspectData} object.
#'
#' @export
#'
addMS2Info <- function(wfobj = NULL,
                       sampleGroup = NULL,
                       suspectList = NULL,
                       MS2param = MS2param(),
                       updateRT = TRUE,
                       updateIntControl = TRUE,
                       save = TRUE) {

  sl <- suspectList@data

  df <- wfobj@results

  if (length(unique(df$name)) != nrow(df)) {
    warning("Duplicate compound names found. Filter workflow object results,
            keeping with desired correspondence between compuond name and feature ID.")
    return()
  }

  if (length(wfobj@data) > 1) {
    if (!is.null(sampleGroup)) {
      pat <- wfobj@data[[sampleGroup]]
    } else {
      warning("More than one sample replicate group in the workflow object.
              Select a representative sample replicate group
              to extract MS2 using the argument sampleGroup.")
      return()
    }
  } else {
    pat <- wfobj@data[[1]]
    if (!is.null(sampleGroup)) {
      pat <- pat[which(patRoon::analysisInfo(pat)$group == sampleGroup), ]
    }
  }

  MS2 <- extractMS2(pat, param = MS2param)

  for (i in seq_len(nrow(df))) {

    xgroup <- df$ID[i]

    xfrag <- MS2[[xgroup]]$MSMS

    if (!base::is.null(xfrag)) {

      sl$hasFragments[sl$name %in% df$name[i]] <- TRUE
      sl$fragments_mz[sl$name %in% df$name[i]] <- base::paste(xfrag$mz, collapse = ";")
      sl$fragments_int[sl$name %in% df$name[i]] <- base::paste(xfrag$intensity, collapse = ";")
      sl$fragments_pre[sl$name %in% df$name[i]] <- base::paste(xfrag$precursor, collapse = ";")

    } else {

       sl$hasFragments[sl$name %in% df$name[i]] <- FALSE

    }

    if (updateRT) sl$rt[sl$name %in% df$name[i]] <- df$rt[i]

    if (updateIntControl) {
      samples <- patRoon::analysisInfo(pat)$analysis
      sl$intControl[sl$name %in% df$name[i]] <- apply(df[i, colnames(df) %in% samples], MARGIN = 1, mean)
    }

  }

  suspectList@data <- sl

  if (TRUE == save) utils::write.csv(sl, file = suspectList@path)

  if (is.character(save)) {
    utils::write.csv(sl, file = base::paste0(dirname(suspectList@path), "/", save, ".csv"), row.names = FALSE)
    suspectList@path <- base::paste0(dirname(suspectList@path), "/", save, ".csv")
  }

return(suspectList)

}
