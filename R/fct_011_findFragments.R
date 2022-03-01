

### ntsFragments -----

#' @title ntsFragments
#'
#' @slot database A \link[data.table]{data.table} with target fragments
#' for searching in MS2 data of a given feature list.
#' @slot settings The parameter settings used for find the fragments.
#' @slot data Extra data produced during screening for fragments.
#' @slot results A \link[data.table]{data.table}
#' with summarized results per sample replicate group.
#'
#' @return An \linkS4class{ntsFragments} object to be added to
#' the workflows slot of an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("ntsFragments",
  slots = c(
    database = "data.table",
    settings = "list",
    data = "list",
    results = "data.table"
  ),
  prototype = list(
    database = data.table::data.table(),
    settings = list(ppm = numeric(),
                    minFeatureIntensity = numeric(),
                    replicates = character(),
                    targets = NULL),
    data = list(),
    results = data.table::data.table()
  )
)




#' @title findFragments
#'
#' @param object An \linkS4class{ntsData} object with features for finding fragments in MS2 data.
#' @param database A \link[data.table]{data.table} with target fragments for matching in the MS2 data.
#' @param replicates A character vector with the name of the sample replicate groups to use for finding fragments.
#' @param targets A character vector with the ID of the features of interest.
#' @param title Optional, title to use for the entry in the workflows slot of the \linkS4class{ntsData} object.
#' @param ppm The mass deviation, in ppm, for matching fragments.
#' @param minFeatureIntensity The minimum intensity of the features for searching for fragments.
#' @param inSilico A character vector with the name of the in silico software
#' for confirmation of fragments formulae. Possible values are "sirius" and "fenform".
#' @param MS2settings A list of parameter settings to use for extracting MS2 data.
#'
#' @return An \linkS4class{ntsData} object with an \linkS4class{ntsFragments} object
#' add to the workflows slot.
#'
#' @export
#'
#' @importFrom data.table merge.data.table as.data.table setnames
#'
findFragments <- function(object,
                          database = NULL,
                          replicates = NULL,
                          targets = NULL,
                          title = NULL,
                          ppm = 15,
                          ppmLoss = 40,
                          minFeatureIntensity = 5000,
                          topMost = 5,
                          inSilico = "sirius",
                          MS2settings = NULL) {

  # object <- dtxcms
  # targets <- c("M207_R932_124", "M233_R941_217", "M748_R883_1287", "M242_R884_264")
  # database <- data.table::fread(paste0(path(object), "//tp_nitro.csv"))
  # object <- fragmentsParameters(object, algorithm = "default")

  checkmate::assertClass(object, "ntsData")

  info <- samplesTable(object)[, .(replicate, blank, polarity)]
  info <- info[!duplicated(info)]
  info$adduct <- sapply(info$polarity, function(x) ifelse(x == "negative", "[M-H]-", "[M+H]+"))

  feats <- features(object)
  if (!is.null(targets)) feats <- feats[id %in% targets, ]

  database <- as.data.table(database)

  if (nrow(database) == 0 | nrow(feats) == 0) {
    warning("Features or database data not found! Check function inputs.")
    return(object)
  }

  dbt <- database
  dbt[type == "induced", type := "loss"]
  setnames(dbt, c("mz", "id"), c("mz_db", "id_hits"), skip_absent = TRUE)

  rpl <- info$replicate
  rpl <- rpl[rpl %in% replicates]

  extra <- list()
  ft_final <- list()

  for (r in rpl) {

    cat(paste0("Screening the ", r, " replicate... \n"))

    ### setup -----

    logVec <- feats[[r]] > minFeatureIntensity
    ft <- feats[logVec, ]
    ft <- setnames(ft, c(r, paste0(r, "_sd")), c("intensity", "intensity_sd"))
    ft$replicate <- r
    ft <- ft[, .(replicate, id, mz, rt, d_ppm, d_sec, intensity, intensity_sd)]

    obj_s <- object[which(replicates(object) %in% r), ft$id]

    pat_s <- obj_s@pat

    if (is.null(MS2settings)) MS2settings <- fragmentsParameters(object)

    MS2 <- generateMS2(
      obj_s,
      algorithm = MS2settings@algorithm,
      settings = MS2settings@settings
    )

    MS2 <- patRoon::filter(
      MS2,
      isolatePrec = list(
        maxIsotopes = 6,
        mzDefectRange = c(-0.005, 0.005),
        intRange = c(0.0005, 2),
        z = 1,
        maxGap = 2
      ),
      retainPrecursorMSMS = TRUE
    )

    ft$hasMS2 <- sapply(ft$id, function(x)  !is.null(MS2[[x]]$MSMS))

    ft <- ft[ft$hasMS2, ]

    ft <- dplyr::mutate(ft,
      hits = 0,
      names = NA_character_,
      id_hits = NA_character_,
      error_ppm = NA_character_
    )

    ft_info <- list()

    #database table, create ionized/deprotonated fragments
    dbt <- dbt[type != "ionized", ]
    addSourceIons <- dbt[type == "loss", ]
    polarity <- info[replicate == r, polarity]
    if (polarity == "positive") addSourceIons$mz_db <- addSourceIons$mz_db + 1.0073
    if (polarity == "negative") addSourceIons$mz_db <- addSourceIons$mz_db - 1.0073
    addSourceIons$type <- "ionized"
    dbt <- rbind(dbt, addSourceIons)


    ### find MS2 -----

    for (i in seq_len(nrow(ft))) {

      idf <- ft$id[i]

      t_x <- as.data.table(MS2[[idf]]$MSMS)
      if ("CE" %in% colnames(t_x)) t_x[, CE := NULL]
      if ("preMZ" %in% colnames(t_x)) t_x[, preMZ := NULL]
      t_x[, ID := seq_len(nrow(t_x))]
      setcolorder(t_x, c("ID", "mz", "intensity", "precursor"))

      if (nrow(t_x) > 0) {

        cat(paste0("Looking at MS2 from feature ", idf, " (", i, "/", nrow(ft), ")... \n"))

        #make full table for iteration
        t_y <- t_x
        colnames(t_y) <- paste0(colnames(t_y), "_y")
        t_y <- t_x[, as.list(t_y), by = t_x]
        t_y <- t_y[, neutralDiff := 0]

        # include precurssor if not present in MS2
        if (!TRUE %in% t_x$precursor) {
          t_p <- t_x[1, ]
          t_p$mz <- ft$mz[i]
          t_p$intensity <- "Fragmented"
          t_p$precursor <- TRUE
          t_p$ID <- max(t_x$ID) + 1
          t_x_temp <- t_x
          colnames(t_x_temp) <- paste0(colnames(t_x_temp), "_y")
          t_p <- t_p[, as.list(t_x_temp), by = t_p]
          t_p <- t_p[, neutralDiff := 0]
          t_y <- rbind(t_y, t_p)
          ft$addedPrecursorIon[i] <- TRUE
        }

        # calculate differences between fragments
        t_y$neutralDiff <- t_y$mz - t_y$mz_y

        #remove negative differences as they repeat
        t_y <- t_y[neutralDiff >= 0, ]

        #convert zeros into fragment ions and remove duplicated neutralDiff
        t_y <- t_y[neutralDiff == 0, neutralDiff := mz]
        t_y <- t_y[!duplicated(t_y$neutralDiff), ]

        #search for MS2 neutral loss between fragments and possible protonated induced cleavage
        t_hits <- fuzzyjoin::difference_inner_join(
          t_y, dbt,
          by = c("neutralDiff" = "mz_db"),
          max_dist = 0.05,
          distance_col = "diff"
        )
        t_hits <- as.data.table(t_hits)
        t_hits$diff <- (t_hits$diff / abs(t_hits$mz_db)) * 1E6

        #ending when t_hits
        if (nrow(t_hits) > 0) {

          #change neutral loss to fragment, when found as a mass trace
          logChangeLoss <- sapply(as.numeric(rownames(t_hits)), function(x, t_hits) {
            t_hits$mz[x] == t_hits$mz_y[x] & t_hits$type[x] == "loss"
          }, t_hits = t_hits)

          t_hits[logChangeLoss, type := "fragment"]

          #remove direct traces by ppm
          t_hits <- t_hits[!(type != "loss" & diff > ppm), ]

          #remove neutral losses by pppmLoss
          t_hits <- t_hits[!(type == "loss" & diff > ppmLoss), ]

          ft$hits[i] <- nrow(t_hits)
          ft$id_hits[i] <- paste(t_hits$id_hits, collapse = "; ")
          ft$names[i] <- paste(t_hits$name, collapse = "; ")
          ft$error_ppm[i] <- paste(round(t_hits$diff, digits = 1), collapse = "; ")
          if (ft$hits[i] > 0) ft_info[[idf]] <- t_hits
        }
      }

      if (exists("t_x")) rm(t_x)
      if (exists("t_y")) rm(t_y)
      if (exists("t_p")) rm(t_p)
      if (exists("t_x_temp")) rm(t_x_temp)
      if (exists("t_hits")) rm(t_hits)
      if (exists("idf")) rm(idf)
    }
    if (exists("i")) rm(i)

    ### in silico -----

    ft_simp <- copy(ft)

    ft <- ft[ft$hits > 0, ]

    ft <- dplyr::mutate(ft,
      checked = 0,
      confirmed = NA_character_,
      confirmed_type = NA_character_
    )

    ft_formulas <- list()

    for (i in seq_len(nrow(ft))) {

      idf <- ft$id[i]

      cat(paste0("Confirming hits for feature ", idf, " (", i, "/", nrow(ft), ")... \n"))

      t <- ft_info[[idf]]

      t$result <- FALSE

      pat_f <- pat_s[, idf]

      if (inSilico == "genform") {
        tformula <- patRoon::generateFormulasGenForm(
          pat_f[1],
          MS2,
          relMzDev = ppm,
          adduct = info[replicate == r, adduct],
          elements = "CHNOPSClF", #CHNOPSClF CHBrClFINOPSSi
          hetero = TRUE,
          oc = TRUE,
          extraOpts = paste0("acc=", ppm, " rej=,", ppm * 2),
          calculateFeatures = TRUE,
          featThreshold = 0,
          featThresholdAnn = 0,
          absAlignMzDev = 0.01,
          MSMode = "both",
          isolatePrec = TRUE,
          timeout = 250,
          topMost = topMost,
          batchSize = 4
        )
      } else {
        tformula <- patRoon::generateFormulasSIRIUS(
          pat_f, MS2,
          relMzDev = ppm,
          adduct = info[replicate == r, adduct],
          elements = "CHBrClFINOPSSi", #CHNOPSClF CHBrClFINOPSSi
          profile = "qtof",
          database = NULL,
          noise = NULL,
          topMost = topMost,
          calculateFeatures = TRUE,
          featThreshold = 0,
          featThresholdAnn = 0,
          absAlignMzDev = 0.01,
          extraOptsGeneral = NULL
        )
      }

      if (length(tformula) > 0) {

        tformula2 <- tformula[[idf]]

        # transform the table
        for (l in seq_len(nrow(tformula2))) {
          frags <- tformula2$fragInfo[[l]]
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

          if (t2$type == "loss") {

            # for Loss looks for both fragments annotation
            t4 <- dplyr::inner_join(fragments, t2, by = c("mz"))
            t4 <- dplyr::select(t4, ion_formula, mz, mz_P, ion_formula_P)

            t5 <- dplyr::inner_join(fragments, t2, by = c("mz" = "mz_y"))
            t5 <- dplyr::select(t5, ion_formula, mz, mz_P, ion_formula_P)

            t6 <- dplyr::inner_join(t4, t5, by = "ion_formula_P")

            if (nrow(t6) > 0) {

              vx <- lapply(t6$ion_formula.x, function(x) CHNOSZ::makeup(x))

              vy <- lapply(t6$ion_formula.y, function(x) CHNOSZ::makeup(x))

              #vy <- lapply(vy, function(x) -x)

              test <- lapply(seq_len(length(vx)), function(x, vx, vy) {
                temp <- vx[[x]]

                for (el in names(vy[[x]])) {
                  temp[names(temp) %in% el] <- temp[names(temp) %in% el] - vy[[x]][names(vy[[x]]) %in% el]
                }
                return(temp)
              }, vx = vx, vy = vy)

              expect <- CHNOSZ::makeup(t2$formula)

              t6$loss <- lapply(test, function(x, expect) {
                x <- x[x > 0]
                x <- setequal(x, expect)
                return(x)
              }, expect = expect)

              t$result[z] <- TRUE %in% t6$loss
            }

            t6$name <- t2$name
            t6$id_check <- t2$id_hits
            t6$type <- t2$type

            tformula5[[z]] <- t6

          } else {

            t3 <- dplyr::inner_join(fragments, t2, by = c("mz"))
            t3 <- dplyr::mutate(t3, result = FALSE)

            if (nrow(t3) > 0) {

              for (j in seq_len(nrow(t3))) {

                insi <- CHNOSZ::makeup(t3$ion_formula[j])

                indb <- CHNOSZ::makeup(t3$formula[j])

                #Adds one more H for direct match of protonated fragments
                if (t3$type[j] == "ionized") indb[names(indb) == "H"] <-  indb[names(indb) == "H"] + 1

                if (TRUE == all.equal(insi, indb)) t3$result[j] <- TRUE

              } #end of j loop

              t$result[z] <- TRUE %in% t3$result
            }

            tformula5[[z]] <- t3
          }

          if (exists("t2")) rm(t2)
          if (exists("t3")) rm(t3)
          if (exists("t4")) rm(t4)
          if (exists("t5")) rm(t5)
          if (exists("t6")) rm(t6)
          if (exists("insi")) rm(insi)
          if (exists("indb")) rm(indb)
          if (exists("vx")) rm(vx)
          if (exists("vy")) rm(vy)
          if (exists("expect")) rm(expect)
          if (exists("j")) rm(j)

        } # enf for z for loop
        if (exists("z")) rm(z)

        ft_formulas[[idf]] <- tformula5
      }

      checked <- t[t$result == TRUE, ]
      ft$checked[i] <-  nrow(checked)
      if (nrow(checked) > 0) {
        ft$confirmed[i] <- paste(as.vector(checked$name), collapse = "; ")
        ft$confirmed_type[i] <-  paste(as.vector(checked$type), collapse = "; ")
      }
      if (exists("checked")) rm(checked)

    } #end of i loop
    if (exists("i")) rm(i)

    ft_final[[r]] <- ft
    extra[[r]] <- list(c(ft_info, ft_formulas))

  } #end of r loop

  ft_final <- rbindlist(ft_final)

  if (nrow(ft_final) == 0) {
    warning("MS2 fragments in targets not found in given features!")
    return(object)
  }

  data <- new("ntsFragments")
  data@database <- database
  data@settings$replicates <- replicates
  data@settings$ppm <- ppm
  data@settings$minFeatureIntensity <- minFeatureIntensity
  data@settings$targets <- targets
  data@data <- extra
  data@results <- ft_final

  if (is.null(title)) title <- "MS2FragmentsScreening"

  object@workflows[[title]] <- data

  return(object)
}
