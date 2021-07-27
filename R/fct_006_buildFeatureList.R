
# TODO improve buildFeatureList function

#' buildFeatureList
#'
#' @param obj An \linkS4class{ntsData} object.
#'
#' @return
#'
#' @export
#'
#' @importMethodsFrom patRoon as.data.frame groupFeatIndex
#' @importFrom stats na.omit
#'
#' @examples
#'
buildFeatureList <- function(obj) {

  obj2 <- obj@patdata

  feat <- patRoon::as.data.frame(obj2, average = TRUE)

  feat <- rename(feat, rt = ret, ID = group)

  rgs <- obj@samples$group

  rg <- unique(rgs)

  sp <- obj@samples$sample

  names(rgs) <- sp

  for (i in seq_len(length(rg))) {
    feat[, paste0(rg[i], "_sd%")] <- apply(
      X = patRoon::as.data.table(obj2, average = FALSE)[, .SD, .SDcols = sp[rgs == rg[i]]],
      MARGIN = 1, function(x) {
        round(ifelse(sd(x) != 0, (sd(x) / mean(x)) * 100, NA), digits = 0)})
  }

  pl <- obj@peaks

  fidx <- patRoon::groupFeatIndex(obj2)
  fidx[fidx == 0] <- NA
  total <- nrow(pl[pl$sample %in% sp[1], ])
  for (i in 2:length(sp)) {
    newtotal <- nrow(pl[pl$sample %in% sp[i], ])
    fidx[i, ] <- fidx[i, ] + total
    total <- total + newtotal
  }
  fidx <-  as.list(fidx)

  evalZeros <- feat[, rg] == 0

  if (TRUE %in% evalZeros) {
    feat <- feat[!apply(evalZeros, MARGIN = 1, function(h) sum(h) == length(rg)), ]
  }

  fidx <- fidx[names(fidx) %in% feat$ID]

  feat$mzmin <- unlist(lapply(fidx, function(h) {min(pl[na.omit(h), "mzmin", drop = TRUE])}))
  feat$mzmax <- unlist(lapply(fidx, function(h) {max(pl[na.omit(h), "mzmax", drop = TRUE])}))
  feat$rtmin <- unlist(lapply(fidx, function(h) {min(pl[na.omit(h), "rtmin", drop = TRUE])}))
  feat$rtmax <- unlist(lapply(fidx, function(h) {max(pl[na.omit(h), "rtmax", drop = TRUE])}))

  feat$dppm <- round(((feat$mzmax - feat$mzmin) / feat$mz) * 1E6, digits = 1)

  feat$width <- round(feat$rtmax - feat$rtmin, digits = 0)

  feat$pIdx <- I(fidx)

  fidx2 <- fidx
  fidx2 <- do.call(rbind, fidx2)
  colnames(fidx2) <- rgs
  fidx2[is.na(fidx2)] <- 0
  fidx2[fidx2 > 0] <- 1
  fidx2 <- sapply(rg, function(h) rowSums(fidx2[, grepl(h, rgs), drop = FALSE]))
  feat$npeaks <- I(lapply(rownames(fidx2), function(x) {
    fidx2[x,, drop = TRUE]
  }))

  feat$hasFilled <- unlist(lapply(fidx, function(x) {
    1 %in% pl[na.omit(x), "is_filled"]
  }))

  feat <- select(feat, ID, mz, rt, dppm, width, everything())

  obj@features <- feat

  return(obj)

}
