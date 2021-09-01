library(ntsIUTA)
rm(list = ls())

feat_sn <- dtxcms2@features

##############################################################

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
#' @importFrom stats mean sd
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
      obj@features <- feat_sn

      return(obj)
}


tic("getsnratio 1030")
obj2 <- getSnRatio(obj, 30)
toc()



total <- 20
# create progress bar
pb <- txtProgressBar(min = 0, max = total, style = 3)
for(i in 1:total){
   Sys.sleep(0.1)
   # update progress bar
   setTxtProgressBar(pb, i)
}
close(pb)
####################################################################

    for (jj in 1:100) {
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

      obj@features <- feat_sn

      }



tic("lapply 100")
base::lapply(jj = 1:100, feat_sn = feat_sn,function(jj, feat_sn) {     

       featEIC <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[jj] - xrt),(feat_sn$rtmax[jj] +xrt)),
                          rtUnit = "sec")
       
       featEIC <- featEIC[!(featEIC$rt>feat_sn$rtmin[jj] & featEIC$rt<feat_sn$rtmax[jj]),]

       #if ((is.null( length(which(featEIC$mz %in% feat_sn$mz))))) print(paste0("index: ", jj))
       #featEIC <- featEIC[ !(featEIC$mz %in% feat_sn$mz),]
       
      #peaks <- peaks_sn[(peaks_sn$rt>feat_sn$rtmin[jj] & peaks_sn$rt<feat_sn$rtmax[jj]),]
      #peaks <- peaks[(peaks_sn$mz>feat_sn$mzmin[jj] & peaks_sn$mz<feat_sn$mzmax[jj]),]

      #peaks <- peaks_sn[(peaks_sn$mz>feat_sn$mzmin[jj] & peaks_sn$mz<feat_sn$mzmax[jj]),]
      #peaks <- peaks[(peaks$rt>(feat_sn$rtmin[jj] - xrt ) & peaks$rt<(feat_sn$rtmax[jj] + xrt)),]






       feat_sn$noise[jj] <- mean(featEIC$i, trim = 0.05)
       #noise_med <- median(featEIC3$i)
       feat_sn$noise_sd[jj] <- sd(featEIC$i)
       feat_sn$sn[jj] <- feat_sn$IN[jj]/feat_sn$noise[jj]


  } )

toc()



  tic("test for loop 100")
    for (jj in 1:100) {

       featEIC <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[jj] - xrt),(feat_sn$rtmax[jj])),
                          rtUnit = "sec")

      featEIC2 <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[jj]),(feat_sn$rtmax[jj] +xrt)),
                          rtUnit = "sec")

      featEIC <- rbind(featEIC, featEIC2)
       
       #featEIC <- featEIC[!(featEIC$rt>feat_sn$rtmin[jj] & featEIC$rt<feat_sn$rtmax[jj]),]

       feat_sn$noise[jj] <- mean(featEIC$i, trim = 0.05)
       feat_sn$noise_sd[jj] <- sd(featEIC$i)
       feat_sn$sn[jj] <- feat_sn$IN[jj]/feat_sn$noise[jj]

  }
 toc() 



tic("lapply 100")
base::lapply(jj = 1:100, feat_sn = feat_sn,function(jj, feat_sn) {

       featEIC <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[jj] - xrt),(feat_sn$rtmin[jj])),
                          rtUnit = "sec")

      featEIC2 <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmax[jj]),(feat_sn$rtmax[jj] +xrt)),
                          rtUnit = "sec")

      featEIC <- rbind(featEIC, featEIC2)
       
       featEIC <- featEIC[!(featEIC$rt>feat_sn$rtmin[jj] & featEIC$rt<feat_sn$rtmax[jj]),]






       feat_sn$noise[jj] <- mean(featEIC$i, trim = 0.05)
       #noise_med <- median(featEIC3$i)
       feat_sn$noise_sd[jj] <- sd(featEIC$i)
       feat_sn$sn[jj] <- feat_sn$IN[jj]/feat_sn$noise[jj]


  } )

toc()

tic("lapply EIC extraction to list 100")
EICextract <- base::lapply(X = 1:100, feat_sn = feat_sn[1:100,], obj = obj, xrt = xrt, function(X, feat_sn, obj, xrt){ 
              extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[X],feat_sn$mzmax[X]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[X] - xrt),(feat_sn$rtmax[X] +xrt)),
                          rtUnit = "sec")
})

toc()



tic("lapply EIC extraction to list 100")
feat_sn_amend <- base::lapply(X = 1:100, feat_sn = feat_sn[1:100,], obj = obj, xrt = xrt, function(X, feat_sn, obj, xrt){ 
              mean( extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[X],feat_sn$mzmax[X]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[X] - xrt),(feat_sn$rtmax[X] +xrt)),
                          rtUnit = "sec") )



                })

toc()




  maxMultiProcess = TRUE
  if (maxMultiProcess) {
    snow <- registered("SnowParam")
    if (snow$workers < detectCores()) {
      snow <- SnowParam(workers = detectCores() - 1,
                        type = "SOCK",
                        exportglobals = FALSE,
                        progressbar = TRUE)
      register(snow, default = TRUE)
    }
  }

tic("bplapply EIC extraction to list 100")

EICextract2 <- BiocParallel::bplapply( X = 1:100, feat_sn = feat_sn[1:100,], obj = obj, xrt = xrt, function(X, feat_sn, obj, xrt){ 

  ntsIUTA::extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[X],feat_sn$mzmax[X]), 
                          ppm = NULL,
                          rt = NULL,
                          rtWindow = c((feat_sn$rtmin[X] - xrt),(feat_sn$rtmax[X] +xrt)),
                          rtUnit = "sec")

}, BPPARAM = BiocParallel::bpparam("SnowParam"))

toc()

#' @title getSnRatio
#' @description Function to extract to s/n ratio of a certain feature.
#' 
#' @param x A \linkS4class{XCMSnExp} object.
#' @param xA A list of \linkS4class{xsAnnotate} objects per replicate group.
#' @param xPat A \linkS4class{featureGroups} object converted by \code{\link{getPatData}}.
#' @param snWindow Time in seconds to expand the peak width for calculation of the signal-to-noise ratio.
#' @param save Logical, set to \code{TRUE} to save the generated \code{fl} object in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#'
#' @return A \code{data.frame} with detailed information for each feature in the given objects.
#' 
#' @export
#' 
#' @import magrittr
#' @importFrom BiocParallel registered SnowParam register bpparam bplapply
#' @importFrom parallel detectCores
#' @importFrom xcms filterMsLevel chromPeaks featureDefinitions
#' @importFrom methods as
#' @importMethodsFrom MSnbase fileNames
#' @importFrom patRoon importFeatureGroupsXCMS3 as.data.table
#' @importFrom dplyr select mutate group_by count filter between all_of everything
#' @importFrom stats quantile sd na.omit
#' @importFrom CAMERA getPeaklist
#' @importFrom stringr str_extract
#' 
#' 
#'
#' @examples
#' 
#' 
#' 
getSnRatio <-  function(obj, xrt = 30)
  
  
  
  
  
  
  feat_sn <- obj@features
  
  #feat_sn$sn <- 1
  
  feat_sn[1,]
  
  xrt <- 30

  peaks_sn <- obj@peaks

  peaks_sn[1,]
  
  
  #get centroids
  exEIC <- extractEIC(dtcent, fileIndex = NULL,
                      mz = 242.1434, ppm = 20,
                      rt = 14.8, rtWindow = 0.5,
                      rtUnit = "min")
  # #1
  exEIC <- extractEIC(dtcent, fileIndex = NULL,
                      mz = 100.0756, ppm = 18.1 ,
                      rt = 585.885, rtWindow =  17,
                      rtUnit = "sec")
  # #2
  exEIC <- extractEIC(dtcent, fileIndex = NULL,
                      mz = 100.1119, ppm = 10.6 ,
                      rt = 585.7506, rtWindow =  20,
                      rtUnit = "sec")
  
  
  feat[1,]
  
  peaks <- obj@peaks
  
  peaks[,peaks$mz ==  100.0756]
  
  dplyr::filter(peaks, mz ==  100.0756, drop = TRUE)
  
  extractEIC()
  xcms::getEIC()
  
  jj <- 27
 #for (jj in 1000:seq_len(nrow(feat))) {
 # for (jj in 1:100) {
 tic("test") 
 rtWin <- ((feat$rtmax[jj] - feat$rtmin[jj])/(2))+xrt
# print(rtWin)
# print(feat$mz[jj])
# print(feat$ret[jj])
  tic("EIC")
  featEIC <- extractEIC(obj = obj,
             fileIndex = NULL,
             mz = feat$mz[jj], 
             ppm = feat$dppm[jj],
             rt = feat$ret[jj],
             #rtWindow = c(feat$rtmax[jj],feat$rtmin[jj]),
             #rtWindow = ((feat$rtmax[jj] - feat$rtmin[jj])/(2))+xrt,
             rtWindow = rtWin,
             rtUnit = "sec"
             )
  toc()
  #featEIC2 <- featEIC[featEIC$rt<feat$rtmin[jj] & featEIC$rt>feat$rtmax[jj]]
  
  featEIC2 <- featEIC[!(featEIC$rt>feat$rtmin[jj] & featEIC$rt<feat$rtmax[jj]),]
  
  featEIC3 <- featEIC2[ !(featEIC2$mz %in% feat$mz),]

  #if (!is.null(length(featEIC3))) print("match")
  
 #}
  
  noise <- mean(featEIC3$i, trim = 0.05)
  noise_med <- median(featEIC3$i)
  noise_sd <- sd(featEIX$i)
  
  sn 
  
  #stats::mad(featEIC3$i)

  #stats::sd(featEIC3$i)
  #stats::var(featEIC3$i)
  toc()
  
  
  test <- feat[1,]
  
  #BiocGenerics::sd(featEIC3$i)
  #BiocGenerics::var(featEIC3$i)
  
  
  #for (jj in 1:seq_len(nrow(feat))) {
  for (jj in 1:1000) {
    
    if (!feat$hasFilled[jj]) {
    
       featEIC <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = feat$mz[jj], 
                          ppm = feat$dppm[jj],
                          rt = feat$ret[jj],
                          #rtWindow = c(feat$rtmax[jj],feat$rtmin[jj]),
                          #rtWindow = ((feat$rtmax[jj] - feat$rtmin[jj])/(2))+xrt,
                          rtWindow = ((feat$rtmax[jj] - feat$rtmin[jj])/(2))+xrt,
                          rtUnit = "sec")
       
       featEIC <- featEIC[!(featEIC$rt>feat$rtmin[jj] & featEIC$rt<feat$rtmax[jj]),]
       featEIC <- featEIC[ !(featEIC2$mz %in% feat$mz),]
       
       feat_sn$noise[jj] <- mean(featEIC$i, trim = 0.05)
       #noise_med <- median(featEIC3$i)
       feat_sn$noise_sd[jj] <- sd(featEIC$i)
       feat_sn$sn[jj] <- feat_sn$IN[jj]/feat_sn$noise[jj]
       
    }
  
    
    
  }
  # 
  # feat$mz %in% featEIC2$mz
  # featEIC2$mz %in% feat$mz
  # 
  # mz2 = feat$mz[jj]
  # ppm2 = feat$dppm[jj]
  # rt2 = feat$ret[jj]
  # #rtWindow2 = c(feat$rtmin[1],feat$rtmax[1])
  # rtWindow2 = (feat$rtmax[jj] - feat$rtmin[jj])/2
  #  
  # 
  # extractEIC(obj = obj,
  #            fileIndex = NULL,
  #            mz = mz2, 
  #            ppm = ppm2,
  #            rt = rt2,
  #            rtWindow = rtWindow2[1],
  #            #rtWindow = ((feat$rtmax[jj] - feat$rtmin[jj])/(2*60)),
  #            rtUnit = "sec")
  # 
  

    for (jj in 1:1000) {
    
    if (!feat$hasFilled[jj]) {
    
    rtrange = c(feat_sn$rtmin[jj],feat_sn$rtmax[jj])
       featEIC <- extractEIC(obj = obj,
                          fileIndex = NULL,
                          mz = c(feat_sn$mzmin[jj],feat_sn$mzmax[jj]), 
                          ppm = NULL,
                          #rt = c(feat_sn$rtmin[jj],feat_sn$rtmax[jj]),
                          rt = NULL,
                          #rt = c(feat_sn$rtmax[jj],feat_sn$rtmin[jj]),
                          rtWindow = c(feat$rtmin[jj],feat$rtmax[jj]),
                          #rtWindow = NULL,
                          #rtWindow = ((feat$rtmax[jj] - feat$rtmin[jj])/(2))+xrt,
                          #rtWindow = ((feat$rtmax[jj] - feat$rtmin[jj])/(2))+xrt,
                          #rtWindow = rtrange,
                          rtUnit = "sec")
       
       featEIC <- featEIC[!(featEIC$rt>feat$rtmin[jj] & featEIC$rt<feat$rtmax[jj]),]
       featEIC <- featEIC[ !(featEIC2$mz %in% feat$mz),]
       
       feat_sn$noise[jj] <- mean(featEIC$i, trim = 0.05)
       #noise_med <- median(featEIC3$i)
       feat_sn$noise_sd[jj] <- sd(featEIC$i)
       feat_sn$sn[jj] <- feat_sn$IN[jj]/feat_sn$noise[jj]
       
    }
  
    
    
  }

