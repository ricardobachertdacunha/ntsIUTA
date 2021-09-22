


#' @title FindIsotopesWithValidationAltered
#' 
#' @description Altered version of the \code{findIsotopesWithValidation} function from \code{CAMERA} package.
#'
#' @param xA An \linkS4class{xsAnnotate} object.
#' @param obj An \linkS4class{ntsData} object containing features.
#' @param sampleidxs The index of samples for the replicate groups of interest.
#' @param maxcharge Maximum number of possible charges.
#' @param ppm The mass deviation, in ppm, allowed for searching isotopes within each component.
#' @param mzabs The mass deviation, in \emph{m/z}, allowed for searching isotopes within each component.
#' @param noise The overall noise level. When given the sn is calculated based on the given noise level.
#' @param intval The intensity value to check isotopes. Default and recommended is "maxo", which is the peak height.
#' @param validateIsotopePatterns Logical, set to \code{TRUE} to validate isotopes using the kegg library.
#'
#' @return An \linkS4class{xsAnnotate} object with annotated isotopes.
#' 
#' @references
#' \insertRef{CAMERA}{ntsIUTA}
#'
#' @importFrom HyperbolicDist is.wholenumber
#' @importClassesFrom xcms xcmsSet
#' @importClassesFrom CAMERA xsAnnotate
#' @importMethodsFrom xcms groupval groups filterMsLevel filterFile
#' @importFrom methods as
#' @importFrom stats median quantile
#' @importFrom BiocParallel registered SnowParam register bpparam bplapply
#' @importFrom parallel detectCores
#' @importFrom dplyr filter between
#' 
FindIsotopesWithValidationAltered <- function(xA = xA,
                                              obj = obj,
                                              sampleidxs = sampleidxs,
                                              maxcharge = 3,
                                              ppm = 40,
                                              mzabs = 0.005,
                                              noise = NULL,
                                              intval = "maxo",
                                              validateIsotopePatterns = TRUE) {
  
  
  
  # object <- xA
  # obj <- obj
  # sampleidxs = sampleidxs
  # maxcharge = 3
  # ppm = 40
  # mzabs = 0.01
  # noise = 350
  # intval = "maxo"
  # validateIsotopePatterns = TRUE
  
  
### checks ------------------------------------------------------------------------------------------------
  
  ## test maxcharge
  if (!HyperbolicDist::is.wholenumber(maxcharge) || maxcharge < 1)
    stop("Invalid argument 'maxcharge'. Must be integer and > 0.\n")
  ## test ppm
  if (!is.numeric(ppm) || ppm < 0)
    stop("Invalid argument 'ppm'. Must be numeric and not negative.\n")
  ## test mzabs
  if (!is.numeric(mzabs) || mzabs < 0)
    stop("Invalid argument 'mzabs'. Must be numeric and not negative.\n")
 
  
### init --------------------------------------------------------------------------------------------------

  object <- xA
  numberOfPS <- length(object@pspectra)
  
  ## scaling
  devppm <- ppm / 1000000
  
  ## number of peaks in pseudospectrum
  numberOfPeaks <- sum(sapply(object@pspectra, length))
  
  if (numberOfPeaks != nrow(obj@features))
    warning("Features in spectra are not the same number as in the ntsData.")
  
  cat("Generating peak matrix!\n")
  
  ftID <- obj@features[, "ID", drop = FALSE]
  
  peaks <- obj@peaks[obj@peaks$sample %in% ntsIUTA::samples(obj)[sampleidxs], ]
  
  peaks <- peaks[, c("mz","rt", "feature", "intensity", "sn"), drop = FALSE]
  
  peaks <- peaks %>%
    dplyr::group_by(feature) %>%
    summarize(mz = mean(mz),
              rt = mean(rt),
              int = mean(intensity),
              sn = mean(sn, na.rm = TRUE))
  
  peaks$sn[is.nan(peaks$sn)] <- 0
  
  ftID <- dplyr::left_join(ftID, peaks, by = c("ID" = "feature"))
  
  for (i in seq_len(nrow(ftID))) {
    if (is.na(ftID$mz[i])) {
      ftID[i, c("ID", "mz", "rt")] <- obj@features[obj@features$ID %in% ftID$ID[i], c("ID", "mz", "rt"), drop = TRUE]
      ftID[i, c("int", "sn")] <- 0
    }
  }
  
  mzValues  <- ftID[, "mz", drop = FALSE]
  
  rtValues <- ftID[, "rt", drop = FALSE]
  
  
  ## isotope matrix
  isoMatrix <- matrix(ncol = 4, nrow = 0)
  colnames(isoMatrix) <- c("mpeak", "isopeak", "iso", "charge")
  
  
### find isotopes -----------------------------------------------------------------------------------------
  
  cat("Run isotope peak annotation\n % finished: ")
  
  df_temp <- CAMERA::getPeaklist(object, intval = "maxo")
  df_temp$ID <- ftID$ID
  
  ## look for isotopes in every pseudospectrum
  numberOfPS <- length(object@pspectra)
  isotopeClusterCounter <- 0
  
  for (psIdx in 1:numberOfPS) {
    
    ## progress
    if (numberOfPS >= 10) if ((psIdx %% round(numberOfPS / 10)) == 0)   cat((psIdx / round(numberOfPS / 10) * 10), ' ')
    
    ## get peak indizes for psIdx-th pseudospectrum
    peakIndeces <- object@pspectra[[psIdx]]
    if (length(peakIndeces) <= 1)
      ## Pseudospectrum has only one peak
      next
    
### calculate isotopes ----------------------------------------------------------------------------------

    isoMatrixForPS.list <- findIsotopesForPSAltered(
      peakIndeces,
      mzValues = ftID$mz[peakIndeces],
      intValues = ftID$int[peakIndeces],
      snValues = ftID$sn[peakIndeces],
      maxcharge = maxcharge,
      devppm = devppm,
      mzabs = mzabs,
      validateIsotopePatterns = validateIsotopePatterns
    )
    
    if (length(isoMatrixForPS.list) == 0)
      ## no isotope cluster found
      next
    
    ## add isotope cluster to isoMatrix
    for (isotopeCluster in 1:length(isoMatrixForPS.list)) {
      isoClusterMatrix <- isoMatrixForPS.list[[isotopeCluster]]
      numberOfClusterPeaks <- nrow(isoClusterMatrix)
      isotopeClusterCounter <- isotopeClusterCounter + 1
      
      ## assign cluster index and box
      isoCluster2 <- cbind(
        isoClusterMatrix[, "peak index"], 
        rep(x = isotopeClusterCounter, times = numberOfClusterPeaks), 
        isoClusterMatrix[, "isotope"], 
        isoClusterMatrix[, "charge"]
      )
      isoMatrix <- rbind(isoMatrix, isoCluster2)
    }
  }


### finishing ---------------------------------------------------------------------------------------------
  
  ## create isotope matrix for peak list
  isotopeMatrix <- vector(mode = "list", length = numberOfPeaks)
  
  for (peakIdx in 1:numberOfPeaks) {
    row <- match(peakIdx, isoMatrix[, "mpeak"])
    if (is.na(row))
      ## no isotope found
      next
    
    ## register isotope data
    isotopeMatrix[[peakIdx]] <- list(
      "y"       = isoMatrix[row[[1]], "isopeak"],
      "iso"     = isoMatrix[row[[1]], "iso"] - 1,
      "charge"  = isoMatrix[row[[1]], "charge"],
      "val"     = 0
    )
  }
  
  ## assign isotope matrix to object
  object@isoID <- isoMatrix
  object@isotopes <- isotopeMatrix
  
  ## out
  cat("\nNumber of isotopes:", nrow(object@isoID), "\n")
  return(object)
}




#' @title findIsotopesForPSAltered
#' @description Function used within \code{FindIsotopesWithValidationAltered}.
#'
#' @param peakIndeces Indices of peaks.
#' @param mzValues \emph{m/z} mzvalues of peaks.
#' @param intValues Intensity values of peaks.
#' @param snValues Signal-to-noise ration of the peaks.
#' @param maxcharge Maximum charge allowed.
#' @param devppm Maximum mass deviation, in ppm.
#' @param mzabs Maximum mass deviation, in \emph{m/z}.
#' @param validateIsotopePatterns Logical, when \code{TRUE} applies validation using the kegg library.
#'
#' @return A list of isotopic annotation for the given peaks.
#'
#' @references
#' \insertRef{CAMERA}{ntsIUTA}
#'
#' @export
#'
#' @importFrom CAMERA compoundQuantiles getIsotopeProportion
#' 
findIsotopesForPSAltered <- function(peakIndeces = peakIndeces,
                                     mzValues = mzValues,
                                     intValues = intValues,
                                     snValues = snValues,
                                     maxcharge = maxcharge,
                                     devppm = devppm,
                                     mzabs = mzabs,
                                     validateIsotopePatterns = TRUE) {
  
  
  
  # mzValues = mzValues[peakIndeces]
  # intValues = intValues[peakIndeces]
  # snValues = snValues[peakIndeces]
  
  ## peakIndeces  - peak indeces
  ## mzValues     - m/z vector, contains all m/z values from specific pseudospectrum
  ## intValues    - int vector, see above
  ## snValues     - signal-to-noise vector, see above
  ## maxcharge    - maximum allowed charge
  ## devppm       - scaled ppm error
  ## mzabs        - absolut error in m/z
  
  ###################################################################################################
  ## algorithm
  ## 
  ## 1) compute triangular matrix of peak differences regarding m/z
  ## 2) compute for all charges possible isotope chains M, M+1, M+2, ...
  ##    (in case of multiple possible successors take the successor with minimum m/z distance from expected m/z)
  ## 3) select the longest chain and proceed with the remaining peaks
  ## 4) check isotope cluster validity and split chains to chain segments
  ## 5) box results
  ## 
  
  ###################################################################################################
  
  ## init
  noiseEstimate <- intValues / snValues
  noiseEstimate[is.nan(noiseEstimate)] <- 0
  intensityMin  <- intValues - noiseEstimate
  intensityMax  <- intValues + noiseEstimate
  
  spectrum <- cbind(peakIndeces, mzValues, intValues, intensityMin, intensityMax)
  colnames(spectrum) <- c("peak index", "mz", "intensity", "intensityMin", "intensityMax")
  
  ## order peaks by m/z
  spectrum <- spectrum[order(spectrum[, "mz"]), ]
  numberOfPeaksHere <- nrow(spectrum)
  
  if (numberOfPeaksHere <= 1) return(list())
  
  ## calculate allowed m/z errors for the given masses; at least mzabs
  mzErrors <- devppm * spectrum[, "mz"]
  #mzErrors[mzErrors < mzabs] <- mzabs
  
  ###################################################################################################
  ## compute m/z difference and possible isotopic connection for every peak pair in pseudospectrum
  isotopeDifference <- 1.0033548378
  expectedDistances <- isotopeDifference / 1:maxcharge
  hitsForCharge_c_p1_p2 <- array(
    dim = c(maxcharge, numberOfPeaksHere - 1, numberOfPeaksHere - 1),
    dimnames = list(c("charge", "peak1", "peak2"))
  )
  
  for (peakIdx in 1:(numberOfPeaksHere - 1)) {
    
    ## create distance matrix
    mzDiff <- spectrum[(peakIdx + 1):numberOfPeaksHere, "mz"] - spectrum[peakIdx, "mz"]
    
    ## compare distances to expected distances
    mzDiffExp <- outer(-expectedDistances, mzDiff, FUN = "+")
    
    ## create hit matrix
    hits <- abs(mzDiffExp) <= mzErrors[[peakIdx]]
    
    ## add hit matrix to big picture
    for (chargeIdx in 1:maxcharge) {
      hitsForCharge_c_p1_p2[chargeIdx, peakIdx, 1:(numberOfPeaksHere - peakIdx)] <- hits[chargeIdx, ]
    }
    
  }#end for peakIdx
  
  ###################################################################################################
  
  ## find and select isotope chains
  goOnSearching <- TRUE
  peakIsProcessed <- rep(x = FALSE, times = numberOfPeaksHere)
  resultChains <- list()
  chainCharges <- list()
  
  while (goOnSearching) {
    candidateChains <- list()
    candidateCharge <- list()
    
    ## check all charges and potential starting peaks to built all possible isotope chains
    for (chargeIdx in 1:maxcharge) {
      
      ## get potential monoisotopic peaks to built all possible isotope chains
      potentialStartPeaks <- which(!peakIsProcessed[1:(numberOfPeaksHere - 1)])
      
      for (peakIdx in potentialStartPeaks) {
        
        ## start new chain with this peak
        peakIdx1 <- peakIdx
        candidateChain <- c(peakIdx1)
        
        ## assemble matching m/z distances to a chain of peaks
        while (
          ## peak is not in any prior chain
          (!peakIsProcessed[peakIdx1]) && 
          
          ## peak is not the last peak in m/z dimension in the current pseudo spectrum
          (peakIdx1 <= (numberOfPeaksHere - 1)) && 
          
          ## there are peaks in the right m/z distance
          any(hitsForCharge_c_p1_p2[chargeIdx, peakIdx1, ], na.rm = TRUE)) {
          
            ## get matching peak index
            matchingIdx <- which(hitsForCharge_c_p1_p2[chargeIdx, peakIdx1, ])
            
            peakIdx2 <- matchingIdx + peakIdx1
            
            if (length(peakIdx2) > 1) {
              ## more than one peak is candidate as successor: take the best matching one regarding m/z
              
              ## get mass of current peak and successor peaks
              mass1   <- spectrum[peakIdx1, "mz"]
              masses2 <- spectrum[peakIdx2, "mz"]
              
              ## take peak with minimum deviation from the expected value
              masses2 <- masses2 - mass1 - expectedDistances[chargeIdx]
              masses2 <- abs(masses2)
              tempIdx <- which.min(masses2)
              peakIdx2 <- peakIdx2[[tempIdx]]
            }
            
            if (peakIsProcessed[peakIdx2])
              ## next peak is already part of a isotope chain --> do not elongate chain
              break
            
            ## elongate chain
            candidateChain <- c(candidateChain, peakIdx2)
            
            ## set peak for next iteration
            peakIdx1 <- peakIdx2
            
          }#end of while loop for current peak chain
        
        if (length(candidateChain) == 1)
          ## a single peak
          next
        
        ## add candidate chain
        candidateChains[[length(candidateChains) + 1]] <- candidateChain
        candidateCharge[[length(candidateCharge) + 1]] <- chargeIdx
      
      }#end for peakIdx
    }#end for chargeIdx
    
    ## select the longest chain
    if (length(candidateChains) == 0) {
      ## no more chains left -> stop searching
      goOnSearching <- FALSE
    } else {
      ## select the longest chain
      maxChainIdx <- which.max(sapply(candidateChains, function(x) length(x)))
      maxChain    <- candidateChains[[maxChainIdx]]
      maxCharge   <- candidateCharge[[maxChainIdx]]
      ## add chain
      resultChains[[length(resultChains) + 1]] <- maxChain
      chainCharges[[length(chainCharges) + 1]] <- maxCharge
      ## mark comprised peaks
      peakIsProcessed[maxChain] <- TRUE
    }
  }#end of while loop for peak chains
  
  if (length(resultChains) == 0)
    ## no isotope cluster found
    return(list())
  
  ###################################################################################################
  ## validate chains
  validatedResultChains <- list()
  validatedChainCharges <- list()
  ## for validation
  cpObj <- compoundQuantiles(compoundLibrary = "kegg")
  maximumIsotopeNumber <- max(cpObj@isotopeSet)
  
  ## available quantiles:
  ## 0.000005 0.999995 0.000010 0.999990 0.000050 0.999950 0.000100 0.999900 0.000500 0.999500 0.001000 0.999000
  ## 0.005000 0.995000 0.010000 0.990000 0.025000 0.975000 0.050000 0.950000 0.100000 0.900000 0.500000
  #quantileLow  <- 0.00500
  #quantileHigh <- 0.99500
  quantileLow  <- 0.00005
  quantileHigh <- 0.99995
  
  for (chainIdx in seq_len(length.out = length(resultChains))) {
    
    ## get data
    charge  <- chainCharges[[chainIdx]]
    chain   <- resultChains[[chainIdx]]
    numberOfIsotopes <- length(chain)
    
    mass <- spectrum[[chain[[1]], "mz"]] * charge
    compoundMassInRange <- mass < cpObj@maxCompoundMass #the maximum mass allowed, default is 1000
    
    if (validateIsotopePatterns & compoundMassInRange) {
      
      ## validate and decompose chain to segments
      monoisotopicPeakIdx <- 1
      isotopePatternMembership <- integer(length = length(chain))
      isotopePatternLabel <- 1
      isotopePatternMembership[[1]] <- isotopePatternLabel
      
      for (peakIdx in 2:length(chain)) {
        
        proportionObservedMin <- spectrum[[chain[[monoisotopicPeakIdx]], "intensityMin"]] / spectrum[[chain[[peakIdx]], "intensityMax"]]
        proportionObservedMax <- spectrum[[chain[[monoisotopicPeakIdx]], "intensityMax"]] / spectrum[[chain[[peakIdx]], "intensityMin"]]
        
        isotopeNumber <- peakIdx - monoisotopicPeakIdx
        
        if (isotopeNumber > maximumIsotopeNumber) {
          ## pattern is too long to be checkable
          proportionExpectedMin <- 0
          proportionExpectedMax <- Inf
        } else {
          ## fetch expected proportion interval
          proportionExpectedMin <- getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = isotopeNumber, mass = mass, quantile = quantileLow)
          proportionExpectedMax <- getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = isotopeNumber, mass = mass, quantile = quantileHigh)
        }
        
        centerObserved = (proportionObservedMin + proportionObservedMax) / 2;
        centerObserved = ifelse(is.nan(centerObserved),0,centerObserved);
        centerExpected  = (proportionExpectedMin + proportionExpectedMax) / 2;
        radiusObserved  = (proportionObservedMax - proportionObservedMin) / 2;
        radiusObserved  = ifelse(is.nan(radiusObserved),0,radiusObserved)
        radiusExpected  = (proportionExpectedMax - proportionExpectedMin) / 2;
        isotopeProportionFits <- abs(centerObserved - centerExpected) <= (radiusObserved + radiusExpected)
        
        if (isotopeProportionFits != TRUE & isotopeProportionFits != FALSE | is.na(isotopeProportionFits)) isotopeProportionFits <- FALSE
        
        if (!isotopeProportionFits) {
          ## isotope proportion does not fit --> start new pattern
          isotopePatternLabel <- isotopePatternLabel + 1
          monoisotopicPeakIdx <- peakIdx
        }
        
        ## assign pattern membership
        isotopePatternMembership[[peakIdx]] <- isotopePatternLabel
        
      }## end for peakIdx
      
      ## extract chain segments
      chainSegments <- list()
      for (chainSegmentLabel in 1:isotopePatternLabel) {
        chainSegment <- chain[which(isotopePatternMembership == chainSegmentLabel)]
        if (length(chainSegment) <= 1)
          next
        
        chainSegments[[chainSegmentLabel]] <- chainSegment
      }## end for chainSegment
    } else {
      chainSegments <- list()
      chainSegments[[1]] <- chain
    }
    
    ## box chain segments
    for (chainSegment in chainSegments) {
      validatedResultChains[[length(validatedResultChains) + 1]] <- chainSegment
      validatedChainCharges[[length(validatedChainCharges) + 1]] <- charge
    }## end for chainSegment
  }## end for chain
  
  ###################################################################################################
  ## assemble list of isotope clusters as result
  isoMatrixForPS.list <- list()
  for (chainIdx in seq_len(length.out = length(validatedResultChains))) {
    ## get data
    chain  <- validatedResultChains[[chainIdx]]
    charge <- validatedChainCharges[[chainIdx]]
    
    ## create and fill isotope cluster
    numberOfIsotopes <- length(chain)
    
    isoClusterMatrix <- matrix(nrow = numberOfIsotopes, ncol = 3)
    colnames(isoClusterMatrix) <- c("peak index", "isotope", "charge")
    for (isoIdx in 1:numberOfIsotopes) {
      ## translate local peak index to global peak index
      peakIdx <- spectrum[chain[[isoIdx]], "peak index"]
      ## add isotope entry
      isoClusterMatrix[isoIdx, ] <- c(peakIdx, isoIdx, charge)
    }
    ## add isotope cluster to result list
    isoMatrixForPS.list[[length(isoMatrixForPS.list) + 1]] <- isoClusterMatrix
  }## end for chainSegment
  
  return(isoMatrixForPS.list)
  
}

### OLD CODE -----

# mzValues  <- object@groupInfo[, "mz", drop = FALSE]
# 
# rtValues  <- object@groupInfo[, "rt", drop = FALSE]



# #Samples in replicate group
# if (nrow(object@xcmsSet@groups) > 0) {
#   
#   intValues <- apply(X = xcms::groupval(object@xcmsSet, value = intval)[, sampleidxs, drop = FALSE],
#                      MARGIN = 1, FUN = function(x) median(x = x, na.rm = TRUE))
#   
#   #Altered the way sn is extracted
#   snValues  <- unlist(lapply(X = 1:nrow(object@groupInfo),
#                              FUN = function(x){
#    filterSamples <- object@xcmsSet@peaks[object@xcmsSet@groupidx[[x]], "sample"] %in% sampleidxs  
#    median(object@xcmsSet@peaks[object@xcmsSet@groupidx[[x]], "sn"][filterSamples], na.rm = TRUE)
#   }))
#   
#   names(snValues) <- names(intValues)
#   
#   #give 0 to sn with int also NA and give intensity 0 to NA intValues
#   intValues[is.na(intValues)] <- 0
#   snValues[is.na(intValues)] <- 0
#   
#   #set threshold values for intensities
#   if (is.null(noise)) noise <- 0
#   intValues[intValues < noise] <- 0
#   snValues[intValues < noise] <- 0
#   
#   # TODO S/N when NA gives error. Fix by calculating the S/N after recursive integration or to all the peaks for other algorithms than xcms.
#   # #collect centroids
#   # if (TRUE %in% is.na(snValues)) {
#   #   
#   #   set <- base::as.data.frame(xcms::groups(object@xcmsSet[,sampleidxs]))
#   #   
#   #   cent <- xcms::filterMsLevel(xcms::filterFile(featData, file = sampleidxs), msLevel. = 1)
#   #   cent <- base::as.data.frame(methods::as(cent, "data.frame"))
#   #   cent$rt <- base::as.numeric(cent$rt)
#   #   cent$mz <- base::as.numeric(cent$mz)
#   #   cent$i <- base::as.numeric(cent$i)
#   #   base::colnames(cent) <- base::c("file", "rt", "mz", "into")
#   #   
#   #   #Enable full parallel processing
#   #   snow <- BiocParallel::registered("SnowParam")
#   #   if (snow$workers < parallel::detectCores())
#   #   {
#   #     snow <- BiocParallel::SnowParam(workers = parallel::detectCores(), type = "SOCK", exportglobals = FALSE)
#   #     BiocParallel::register(snow, default = TRUE)
#   #   }
#   #   
#   #   snValues <- base::unlist(BiocParallel::bplapply(X = 1:base::length(snValues), cent=cent, set=set, intValues=intValues, snValues=snValues,
#   #                                      function(sn, cent, set, intValues, snValues){
#   #     if (base::is.na(snValues[sn])) {
#   #       temp_cent <- cent[cent$mz >= set$mzmed[sn]-(5/1E6*set$mzmed[sn]) &
#   #                         cent$mz <= set$mzmed[sn]+(5/1E6*set$mzmed[sn]) &
#   #                         cent$rt >= set$rtmin[sn]-240 & cent$rt <= set$rtmax[sn]+240, ]
#   #       temp_cent <- dplyr::filter(temp_cent, !dplyr::between(temp_cent$rt, set$rtmin[sn], set$rtmax[sn]))
#   #       temp_sn <-   intValues[sn]/base::round(stats::quantile(temp_cent$into, probs = base::seq(0,1,0.25))[2], digits = 0)
#   #       base::names(temp_sn) <- base::names(intValues[sn])
#   #       sn <- temp_sn
#   #     } else {
#   #       sn <- snValues[sn]
#   #     }
#   #   }, BPPARAM = BiocParallel::bpparam("SnowParam")))
#   # }
#   
#   snValues <- round(snValues, digits = 1)
#   
#   snValues[snValues < 3] <- 0
#   intValues[snValues < 3] <- 0
#   
#   #last check for NA, if present remove both int and sn values
#   snValues[is.na(snValues)] <- 0
#   intValues[is.na(snValues)] <- 0
#   
# } else {
#   
#   ## one sample case
#   intValues <- object@groupInfo[, intval, drop = FALSE]
#   snValues  <- object@groupInfo[, "sn", drop = FALSE]
#   snValues[is.na(snValues)] <- 3
#   snValues[is.na(intValues)] <- 0
#   intValues[is.na(intValues)] <- 0
#   
# }
