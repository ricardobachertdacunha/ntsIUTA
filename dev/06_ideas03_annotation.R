


# View(df[df$comp == "1_Blank", ])
#   
# View(components(obj,
#                 samples = NULL,
#                 ID = "M239_R936_1314",
#                 mz = NULL, ppm = NULL,
#                 rt = NULL, rtWindow = 1, rtUnit = "min",
#                 compNumber = NULL,
#                 entireComponents = TRUE,
#                 onlyAnnotated = TRUE,
#                 onlyRelated = TRUE))
# 
# plotFeaturePeaks(obj = dtxcms2,
#                  samples = NULL,
#                  ID = "M123_R657_72",
#                  mz = NULL,
#                  rt = NULL,
#                  rtUnit = "min",
#                  ppm = NULL,
#                  rtWindow = NULL,
#                  interactive = TRUE)




makeFeatureComponents <- function(obj) {
  
  replicateList <- split(obj@annotation$comp, obj@annotation$comp$group)
  
  #prepare
  repliComponents <- sapply(replicateList, function(x) {
    
    tbl <- dplyr::select(x, ID, rt, mz, comp, isogroup, isoclass, isonr, adductclass, Mion)
    
    tbl$group <- tbl$ID
    
    tbl$ID <- seq_len(nrow(tbl))
    
    tbl$adduct_ion <- unlist(tbl$adductclass)
    
    tbl$comp <- as.numeric(stringr::str_extract(tbl$comp, "[^_]+"))
    
    tbl$isogroup <- as.numeric(stringr::str_extract(tbl$isogroup, "[^_]+"))
    
    tbl$charge <- sub('.*(?=.{2}$)', '', tbl$isoclass, perl = T)
    
    tbl$charge <- gsub("\\]\\+|\\]\\-", 1, tbl$charge)
    
    tbl$charge <- as.numeric(stringr::str_extract(tbl$charge, "[0-9]"))
    
    tbl$score <- 100
    
    tbl <- dplyr::rename(tbl, ret = rt,
                         cliqueGroup = comp,
                         neutralMass = Mion)
    
    tbl <- tbl[, c("ID", "ret", "mz", "cliqueGroup",
                   "isogroup", "isonr", "charge",
                   "adduct_ion",
                   "score", "neutralMass", "group")]
    
    ret <- unname(split(data.table::data.table(tbl), by = c("cliqueGroup", "neutralMass"), keep.by = TRUE))
    
    return(ret)
    
  })
  
  ## When needed per sample not sample replicate group
  # sp <- samples(obj)[sampleGroups(obj) %in% names(replicateList)]
  # rg <- sampleGroups(obj)[samples(obj) %in% sp]
  # featComponents <- list()
  # for (i in seq_len(length(sp))) featComponents[[sp[i]]] <- repliComponents[[rg[i]]]
  
  featureComponents <- repliComponents
  
  ## get info from patdata
  fGroups <- obj@patdata
  
  ftindex <- groupFeatIndex(fGroups)
  gNames <- names(fGroups)
  gInfo <- groupInfo(fGroups)
  anas <- analyses(fGroups)
  gTable <- groupTable(fGroups)
  
  
  
  
  ### utils functions -----
  
  pruneList <- function(l, checkEmptyElements = FALSE, checkZeroRows = FALSE)
  {
    ret <- l[!sapply(l, is.null)]
    if (checkEmptyElements)
      ret <- ret[lengths(ret) > 0]
    if (checkZeroRows)
      ret <- ret[sapply(ret, nrow) > 0]
    return(ret)
  }
  
  numGTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x > y
  
  numEQ <- function(x, y, tol = sqrt(.Machine$double.eps)) abs(x - y) <= tol
  
  numLTE <- function(x, y, tol = sqrt(.Machine$double.eps)) numEQ(x, y, tol) | x < y
  
  allSame <- function(l, func = identical)
  {
    if (length(l) > 1)
    {
      if (all(is.na(l)))
        return(TRUE)
      if (any(is.na(l)))
        return(FALSE)
      
      return(all(sapply(l[-1], func, l[[1]])))
    }
    
    return(TRUE)
  }
  
  calculateComponentIntensities <- function(comps, fGroups)
  {
    getGroupInt <- function(grp)
    {
      ints <- fGroups[[grp]]
      return(mean(ints[ints != 0]))
    }
    return(lapply(comps, function(cmp)
    {
      cmp <- copy(cmp)
      cmp[, intensity := sapply(group, getGroupInt)]
      cmp[, intensity_rel := intensity / max(intensity)]
      return(cmp[])
    }))
  }
  
  
  
  ### from patRoon -----
  
  
  
  # NO need as groups are included already
  
  # featureComponents <- Map(featureComponents, split(ftindex, seq_len(nrow(ftindex))), f = function(fCmpL, fti)
  # {
  #   # assign group names and prune features without groups
  #   fti <- unlist(fti)
  #   fCmpL <- lapply(fCmpL, function(cmp)
  #   {
  #     set(cmp, j = "group", value = gNames[match(cmp$ID, fti)])
  #     return(cmp[!is.na(group)])
  #   })
  #   fCmpL <- pruneList(fCmpL, checkZeroRows = TRUE)
  #   return(fCmpL)
  # })
  
  
  
  
  # fCMP: unique feature component ID within an analysis
  ## Change for each sample replicate group
  cmpTab <- rbindlist(lapply(featureComponents, rbindlist, idcol = "fCMP"), idcol = "analysis")
  
  
  # fCMPID: unique identifier throughout all analyses (sample replicate groups)
  
  #cmpTab[, fCMPID := paste0(match(analysis, analyses(fGroups)), "-", fCMP)]
  cmpTab[, fCMPID := paste0(match(analysis, unique(sampleGroups(obj))), "-", fCMP)]
  
  # NOTE: abundance only takes assigned features into account, as unassigned won't be present
  cmpTab[!is.na(adduct_ion), abundance := sapply(adduct_ion, function(a) sum(a == adduct_ion)) / .N, by = "group"]
  
  
  # this shouldn't happen for now!!!!
  # Filter adducts not abundantly assigned to same feature group
  #cmpTab <- cmpTab[is.na(abundance) | numGTE(abundance, relMinAdductAbundance)]
  
  
  if (adductConflictsUsePref && length(prefAdducts) > 0)
  {
    # for fGroups with features that have a preferential adduct: remove all others or ones that are lower ranked
    cmpTab[!is.na(adduct_ion), prefInd := match(adduct_ion, prefAdducts, nomatch = length(prefAdducts) + 1),
           by = "group"]
    # NOTE: below leaves features untouched if none of the adducts are preferential, since prefInd will be the same
    # for all and thus all are equal to min(prefInd)
    cmpTab[, keep := is.na(adduct_ion) | prefInd == min(prefInd), by = "group"]
    cmpTab <- cmpTab[is.na(adduct_ion) | keep == TRUE][, keep := NULL]
  }
  
  
  
  
  # Only the most abundantly assigned adduct for each feature group. NOTE: if preferential adducts were selected above
  # then these are now always the most abundant.
  # TODO UNDONE: handle ties?
  #solve ties by mass error
  # cmpTab$massdiff <- rules$massdiff[match(cmpTab$adduct_ion, rules$name)]
  # 
  # cmpTab$masserror <- abs((cmpTab$mz - cmpTab$massdiff) - cmpTab$neutralMass)
  
  cmpTab[!is.na(adduct_ion), keep :=
           uniqueN(adduct_ion) == 1 | adduct_ion == adduct_ion[which.max(abundance)], by = "group"]
  
  # cmpTab[!is.na(adduct_ion), keep :=
  #          uniqueN(abundance) == 1 | adduct_ion == adduct_ion[which.min(masserror)], by = "group"]
  
  
  cmpTab <- cmpTab[is.na(adduct_ion) | keep == TRUE][, keep := NULL]
  
  
  
  
  # Start making group components; for each feature group:
  # - find all feature components that this "parent group" is in
  # - assume that this group and all groups in the feature components are related
  # - mark all these groups so they won't be used as parent groups next iterations
  # NOTE: we can do the latter since groups present in multiple components will be removed afterwards, hence nothing
  # will be 'missed'.
  usedGroups <- setNames(rep(FALSE, length(gNames)), gNames)
  
  linkedFGs <- rbindlist(sapply(gNames, function(gn) {
    
    if (usedGroups[gn]) return(NULL)
    
    IDs <- unique(cmpTab[group == gn]$fCMPID)
    ct <- cmpTab[fCMPID %chin% IDs]
    
    ##### Solved by given default neutral mass to features.
    # special case: if parent has NA neutralMass, ensure all others also have and vice versa. This ensures a split
    # between those with NAs and those without. If the parentGroup has mixed NAs, the ones with NA will be ignored.
    ct <- if (any(!is.na(ct[group == gn]$neutralMass))) ct[!is.na(neutralMass)] else ct[is.na(neutralMass)]
    
    usedGroups[ct$group] <<- TRUE
    
    if (!allSame(ct$neutralMass, function(x1, x2) numLTE(abs(x1 - x2), absMzDev))) {
      # conflict in neutral masses --> only retain those with most abundant neutral mass
      # UNDONE: handle ties?
      
      nrowPrior <- nrow(ct)
      colsPrior <- copy(names(ct)) # NOTE: need a copy: https://stackoverflow.com/a/15913648
      
      # since we must work with mass tolerances, first cluster presumably masses together
      hc <- fastcluster::hclust(dist(ct$neutralMass))
      ct[, clust := cutree(hc, h = absMzDev)]
      ct[, clust_size := .N, by = "clust"]
      
      for (cnf in NMConflicts)
      {
        bestCL <- integer()
        if (cnf == "preferential" && length(prefAdducts) > 0)
        {
          ct[, prefAdductMatch := match(adduct_ion, prefAdducts, nomatch = length(prefAdducts) + 1),
             by = "clust"]
          ct[, topRankedMatch := min(prefAdductMatch), by = "clust"]
          
          if (uniqueN(ct$topRankedMatch) == max(ct$clust)) # all clusters ranked differently?
            bestCL <- ct[which.min(topRankedMatch)]$clust
        }
        else if (cnf == "mostAbundant")
          bestCL <- ct[which.max(clust_size)]$clust
        else # mostIntense
        {
          ct[, intensity := mapply(analysis, group, FUN = function(a, g) gTable[[g]][match(a, anas)])]
          ct[, maxClustInt := max(intensity) / clust_size, by = "clust"]
          bestCL <- ct[which.max(intensity)]$clust
        }
        
        if (length(bestCL) == 1)
        {
          ct <- ct[clust == bestCL]
          break
        }
      }
      if (nrow(ct) == nrowPrior)
      {
        # Could not resolve neutral mass conflict --> just default to first cluster...
        ct <- ct[clust == 1]
      }
      
      ct <- ct[, colsPrior, with = FALSE] # remove temporary work columns
    }
    
    return(ct)
    
  }, simplify = FALSE), idcol = "parentGroup")
  
  # collapse features: only retain one row per feature group
  linkedFGs <- unique(linkedFGs, by = c("parentGroup", "group"))
  
  dups <- function(v) names(which(table(v) > 1))
  
  # remove feature groups that occur in multiple to be components
  linkedFGs <- linkedFGs[!group %chin% dups(group)]
  
  # prepare for components
  cols <- intersect(c("parentGroup", "group", "neutralMass", "isonr", "charge", "adduct_ion"),
                    names(linkedFGs))
  
  linkedFGs <- linkedFGs[, cols, with = FALSE]
  
  linkedFGs[, c("ret", "mz") := gInfo[group, c("rts", "mzs")]]
  
  setcolorder(linkedFGs, c("group", "ret", "mz"))
  
  comps <- split(linkedFGs, by = "parentGroup", keep.by = FALSE)
  
  # Remove any fGroups from components with equal adducts (unless assigned to different isotope)
  # if (!is.null(linkedFGs[["isonr"]])) {
  #   comps <- lapply(comps, function(ct) ct[is.na(adduct_ion) | !paste0(adduct_ion, isonr) %chin% dups(paste0(adduct_ion, isonr))])
  # } else {
  #   comps <- lapply(comps, function(ct) ct[is.na(adduct_ion) | !adduct_ion %chin% dups(adduct_ion)])
  # }
  
  # NOTE: minSize should be >= 1 to filter out empty components
  comps <- comps[sapply(comps, nrow) >= 1] #minSize
  
  calculateComponentIntensities <- function(comps, fGroups)
  {
    getGroupInt <- function(grp)
    {
      ints <- fGroups[[grp]]
      return(mean(ints[ints != 0]))
    }
    return(lapply(comps, function(cmp)
    {
      cmp <- copy(cmp)
      cmp[, intensity := sapply(group, getGroupInt)]
      cmp[, intensity_rel := intensity / max(intensity)]
      return(cmp[])
    }))
  }
  
  
  if (length(comps) > 0) {
    comps <- calculateComponentIntensities(comps, fGroups)
    names(comps) <- paste0("CMP", seq_along(comps))
  }
  
  cInfo <- data.table(name = names(comps), cmp_ret = sapply(comps, function(cmp) mean(cmp$ret)),
                      cmp_retsd = sapply(comps, function(cmp) sd(cmp$ret)),
                      neutral_mass = sapply(comps, function(cmp) mean(cmp$neutralMass)),
                      size = sapply(comps, nrow))
  
  
  return(callNextMethod(.Object, featureComponents = featureComponents, components = comps,
                        componentInfo = cInfo, ...))
  
  
  
  try <- new("componentsCliqueMS", fGroups = new("featureGroupsXCMS3", features = new("featuresXCMS3")),
             absMzDev = 0.008,
             minSize = 2,
             relMinAdductAbundance = 1,
             adductConflictsUsePref = TRUE,
             NMConflicts = c("mostIntense"),
             prefAdducts = c("[M+H]+"),
             featureComponents = featComponents)
  
  try <- new("componentsCliqueMS", fGroups = new("featureGroupsXCMS3", features = new("featuresXCMS3")))
  try@cliques <- list()
  try@featureComponents <- featureComponents
  try@components <- comps
  try@componentInfo <- cInfo
  
  
  
  
  
  
  
}




consolidateAnnotation <- function(obj,
                                  rule = c("preferable",
                                           "mostIntense",
                                           "mostAbundant"),
                                  adduct = ("[M+H]+")) {
  
  rule <- "mostIntense"
  
  ft <- obj@features[, c("ID", "mz", "rt")]
  
  ft$compN <- 0
  
  ft$Mion <- NA
  
  ft$adduct <- NA_character_
  
  ft$multiMion <- FALSE
  
  comp <- obj@annotation$comp
  
  rg <- unique(comp$group)
  
  makecomp <- ft[, "ID", drop = FALSE]
  for (r in seq_len(length(rg))) {
    comp_g <- comp[comp$group %in% rg[r], c("ID", "comp"), drop = FALSE]
    colnames(comp_g)[2]    <- rg[r]
    makecomp <- dplyr::left_join(makecomp, comp_g, by = "ID")
  }
  makecomp$compN <- 0
  
  
  for (i in seq_len(nrow(ft))) {
    
    temp <- comp[comp$ID %in% ft$ID[i], ]
    
    compN <- unique(temp$comp)
    
    mion <- unique(temp$Mion)
    
    iso <- unique(temp$isoclass)
    
    adu <- unlist(unique(temp$adductclass))
    
    
    
    
    
    
    
    
    #if (length(mion) > 1) ft$multiMion[i] <- TRUE
    
    if (length(mion) > 1) {
      ft$multiMion[i] <- TRUE
      if (rule == "mostIntense") {
        mion <- temp$Mion[temp$intensity == max(temp$intensity)] 
        adu <- unlist(temp$adductclass[temp$intensity == max(temp$intensity)])
      }
    }
    
    ft$Mion[i] <- mion
    
    ft$adduct[i] <- unlist(adu)
    
  }
}
