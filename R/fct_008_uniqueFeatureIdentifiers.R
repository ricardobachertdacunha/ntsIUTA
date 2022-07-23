

#' @title makeUFI
#'
#' @description Function to collapse istopes to the respective feature,
#' building a unique feature identifier (UFI) based on the calculated
#' neutral monoisotopic mass of the (de)protonated ion.
#'
#' @param object A \linkS4class{ntsData} object containing annotated features.
#' @param save Logical, set to \code{TRUE} to save updated
#' \linkS4class{ntsData} object in the \strong{rdata} folder.
#' Note that \code{TRUE} overwrites the existing \linkS4class{ntsData} object.
#' Optionally, a character string can be given instead of \code{TRUE}
#' to be used as file name, avoiding overwriting.
#'
#' @export
#'
#' @importFrom checkmate assertClass
#' @importFrom dplyr select left_join
#' @importFrom pbapply pblapply
#' @importFrom data.table copy
#'
makeUFI <- function(object, save = FALSE) {

  checkmate::assertClass(object, "ntsData")

  feats <- features(object)

  if ("ufi" %in% colnames(feats)) {
    feats[, ufi := NULL]
    feats[filterTag %in% "isotope", filterTag := NA_character_]
  }

  dir <- system.file(package = "ntsIUTA", dir = "extdata")
  if (polarity(object) %in% "negative") {
    db_adducts <- fread(paste0(dir, "/adducts_neg.csv"), header = TRUE)
    prefAdduct <- db_adducts$name
  } else {
    db_adducts <- fread(paste0(dir, "/adducts_pos.csv"), header = TRUE)
    prefAdduct <- db_adducts$name
  }

  if (!"neutralMass" %in% colnames(feats)) {
    warning("Annotation not found in the given ntsData object.")
    return(object)
  }

  cat("Merging isotopic peaks... \n")

  unify <- feats[, .(id, component, neutralMass, d_ppm, rt, d_sec, charge, adduct_ion)]
  unify <- unify[adduct_ion %in% prefAdduct, ]
  unify[is.na(charge), charge := 1]
  setnames(unify,
    c("neutralMass", "d_ppm", "rt", "d_sec", "charge", "adduct_ion"),
    c("m", "d", "r", "t", "z", "a")
  )
  unify[, `:=`(i1 = 0, i2 = 0, i3 = 0, i4 = 0, i5 = 0)]
  unify[, `:=`(i1_sd = 0, i2_sd = 0, i3_sd = 0, i4_sd = 0, i5_sd = 0)]
  unify[, `:=`(f1 = 0, f2 = 0, f3 = 0)]
  unify[, `:=`(features = list())]

  unifyList <- pbapply::pblapply(unify$id, function(x, unify, feats) {

    temp <- unify[id %in% x, ]
    mainft <- feats[id %in% x, ]

    isotopes <- feats[monoiso %in% mainft$monoiso & isonr > 0, ]

    if (nrow(isotopes) > 0) {
      mainpks <- peaks(object, samples = samples(object)[!blanks(object) == replicates(object)],  targets = x)
      for (iso in isotopes$id) {
        isoNumber_sd <- paste0("i", isotopes[id %in% iso, isonr], "_sd")
        isoNumber <- paste0("i", isotopes[id %in% iso, isonr])
        pks <- peaks(object, samples = samples(object)[!blanks(object) == replicates(object)],  targets = iso)
        pks <- pks[intensity > 0, ]
        if (nrow(pks) > 0) {
          temp[, features := list(c(unlist(unlist(temp$features)), iso))]
          pks <- dplyr::left_join(pks[, .(sample, intensity)], mainpks[, .(sample, intensity)], by = "sample")
          pks_mean <- round(mean(pks$intensity.x / pks$intensity.y, na.rm = TRUE) * 100, digits = 1)
          pks_sd <- round(sd(pks$intensity.x / pks$intensity.y, na.rm = TRUE) * 100, digits = 1)
          temp[, (isoNumber) := pks_mean]
          temp[, (isoNumber_sd) := pks_sd]
        }
      }
    }

    return(temp)
  }, unify = unify, feats = feats)

  unify <- rbindlist(unifyList)

  ufi <- paste(
    "m", round(unify$m, digits = 4),
    "_d", round(unify$d, digits = 0),
    "_r", round(unify$r, digits = 0),
    "_t", round(unify$t, digits = 0),
    "_z", round(unify$z, digits = 0),
    "_a", unify$a,
    "_i", round(unify$i1, digits = 2),
    "/", round(unify$i2, digits = 2),
    "/", round(unify$i3, digits = 2),
    "/", round(unify$i4, digits = 2),
    "/", round(unify$i5, digits = 2),
    "_f", round(unify$f1, digits = 3),
    "/", round(unify$f2, digits = 3),
    "/", round(unify$f3, digits = 3),
    sep = ""
  )

  unify[, ufi := ufi]

  unify <- select(unify, ufi, everything())

  object@unified <- unify

  org_feats <- features(object)

  org_feats <- copy(left_join(org_feats, unify[, .(ufi, id)], by = "id"))

  if (!"filterTag" %in% colnames(org_feats)) org_feats[, filterTag := NA_character_]

  org_feats[is.na(ufi), filterTag := "isotope"]

  org_feats <- select(org_feats, ufi, everything())

  object@features <- org_feats

  pks <- peaks(object)

  if ("ufi" %in% colnames(pks)) pks[, ufi := NULL]

  pks <- copy(dplyr::left_join(pks, org_feats[, .(ufi, id)], by = c("feature" = "id")))

  object@peaks <- pks

  cat("Done! \n")

  if (is.logical(save)) if (save) saveObject(object = object)

  if (is.character(save)) saveObject(object = object, filename = save)

  return(object)
}
