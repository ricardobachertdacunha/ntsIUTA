#' @title getSnRatio
#' @description Function to calculate the s/n ratios for the current feature set. The features \code{data.frame} is amended with columns of the respective noise, noise sd and signal to noise ratio. 
#' The noise is calculated from the mean value of the intensities of the surrounding peaks within the same mass window.
#' @param obj A \linkS4class{XCMSnExp} object.
#' @param xrt The size of the time interval before and after the feature for determination of the surrounding noise. Value is in seconds.
#'
#' @return A \linkS4class{XCMSnExp} object.
#' 
#' @export
#' 
#' @importFrom stats sd
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' 
#'
#' @examples
#' 
#' 
#' 
getSnRatio <-  function(obj, xrt = 30) {
  

    feat_sn <- obj@features
    
    if (!("noise" %in% colnames(feat_sn))) {
      feat_sn$noise <- NA
    }
    
    if (!("noise_sd" %in% colnames(feat_sn))) {
      feat_sn$noise_sd <- NA
    }
    if (!("sn" %in% colnames(feat_sn))) {
      feat_sn$sn <- NA
    }
    
    
    feats_sn2  <- feats_sn[feats_sn$isFiltered == TRUE,]
    feats_sn <- feats_sn[feats_sn$isFiltered == FALSE,]
    
    pb <- txtProgressBar(min = 0, max = nrow(feat_sn), style = 3)
    for (jj in  seq_len(nrow(feat_sn))) {
       featEIC <- extractEIC(obj = obj,
                            fileIndex = NULL,
                           mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                            ppm = NULL,
                            rt = NULL,
                            rtWindow = c((feat_sn$rtmin[jj] - xrt),(feat_sn$rtmax[jj] +xrt)),
                            rtUnit = "sec")
       
       featEIC <- featEIC[!(featEIC$rt>feat_sn$rtmin[jj] & featEIC$rt<feat_sn$rtmax[jj]),]

       feat_sn$noise[jj] <- round(mean(featEIC$i, trim = 0.05), digits = 0)
       feat_sn$noise_sd[jj] <- round(sd(featEIC$i), digits = 0)
       feat_sn$sn[jj] <- round(feat_sn$IN[jj]/feat_sn$noise[jj], digits = 0)

      setTxtProgressBar(pb, jj)
    }
    
    feat_sn <- rbind(feat_sn, feat_sn2)
    feat_sn <- feats_sn[order(feat_sn$mz),]
    
      obj@features <- feat_sn

      return(obj)
}

