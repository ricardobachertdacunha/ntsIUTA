


#' @title FindIsotopesWithValidationAltered
#' 
#' @description Altered version of the \code{findIsotopesWithValidation} function from \code{CAMERA} package.
#'
#' @param object An \linkS4class{xsAnnotate} object.
#' @param featData An \linkS4class{XCMSnExp} object containing features.
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
#' @export
#'
#' @importFrom HyperbolicDist is.wholenumber
#' @importFrom xcms groupval groups filterMsLevel filterFile
#' @importFrom methods as
#' @importFrom stats median quantile
#' @importFrom BiocParallel registered SnowParam register bpparam bplapply
#' @importFrom parallel detectCores
#' @importFrom dplyr filter between
#' 
#'
#' @examples
#' 
#' 
#' 
FindIsotopesWithValidationAltered <- function(object, featData = featData, sampleidxs = sampleidxs, maxcharge = 3,
                                              ppm = 50, mzabs = 0.01, noise = NULL,
                                              intval = "maxo", validateIsotopePatterns = TRUE) {
  
  
  
  # object <- xA
  # featData <- featData
  # sampleidxs = sampleidxs
  # maxcharge = 3
  # ppm = 40
  # mzabs = 0.01
  # noise = 350
  # intval = "maxo"
  # validateIsotopePatterns = TRUE
  
  #searches in every pseudospectrum mass differences, which match isotope distances
  
  ################################################################################################
  ## sanity checks
  
  ## test maxcharge  
  if(!HyperbolicDist::is.wholenumber(maxcharge) || maxcharge < 1)
    stop("Invalid argument 'maxcharge'. Must be integer and > 0.\n")
  ## test ppm
  if(!base::is.numeric(ppm) || ppm < 0)
    stop("Invalid argument 'ppm'. Must be numeric and not negative.\n")
  ## test mzabs
  if(!base::is.numeric(mzabs) || mzabs < 0)
    stop("Invalid argument 'mzabs'. Must be numeric and not negative.\n")
  ## test intval
  #intval <- match.arg(intval)
  
  ################################################################################################
  ## init  
  numberOfPS <- base::length(object@pspectra)
  
  ## scaling
  devppm <- ppm / 1000000
  

  ## Check if object have been preprocessed with groupFWHM
  if(numberOfPS < 1) {
    base::cat("xsAnnotate contains no pseudospectrum. Regroup all peaks into one!\n")
    numberOfPS <- 1
    object@pspectra[[1]] <- base::seq(1:base::nrow(object@groupInfo))
    object@psSamples  <- 1
  }
  
  ## number of peaks in pseudospectrum
  numberOfPeaks <- base::sum(base::sapply(object@pspectra, length))
  
  ## get mz, rt, and intensity values from peaktable
  
  ## "mz"     "mzmin"  "mzmax"  "rt"     "rtmin"  "rtmax"  "into"   "intb"   "maxo"   "sn"
  ## "egauss" "mu"     "sigma"  "h"      "f"      "dppm"   "scale"  "scpos"  "scmin"  "scmax"  "lmin"   "lmax"   "sample"
  
  base::cat("Generating peak matrix!\n")
  mzValues  <- object@groupInfo[, "mz", drop=FALSE]
  rtValues  <- object@groupInfo[, "rt", drop=FALSE]
  
  
  
  #Edit, selects the samples given by sampleidxs
  if(base::nrow(object@xcmsSet@groups) > 0){
    
    ## multiple sample or grouped single sample
    if(base::is.na(object@sample[1])){
      index <- 1:base::length(object@xcmsSet@filepaths) #not used as sample is always given by the for loop
    }else{
      index <- object@sample
    }
    
    #intValues <- groupval(object@xcmsSet,value=intval)[,index,drop=FALSE]
    intValues <- base::apply(X = xcms::groupval(object@xcmsSet,value=intval)[,sampleidxs,drop=FALSE], MARGIN = 1,
                       FUN = function(x){stats::median(x = x, na.rm = TRUE)})
    
    #Altered the way sn is extracted
    snValues  <- base::unlist(base::lapply(X = 1:base::nrow(object@groupInfo),
                               FUN = function(x){
                                                 filterSamples <- object@xcmsSet@peaks[object@xcmsSet@groupidx[[x]], "sample"] %in% sampleidxs  
                                                 stats::median(object@xcmsSet@peaks[object@xcmsSet@groupidx[[x]], "sn"][filterSamples], na.rm = TRUE)
                                                 }))
    
    base::names(snValues) <- base::names(intValues)
    
    #give 0 to sn with int also NA and give intensity 0 to NA intValues
    intValues[base::is.na(intValues)] <- 0
    snValues[base::is.na(intValues)] <- 0
    
    #set threshold values for intensities
    intValues[intValues < noise] <- 0
    snValues[intValues < noise] <- 0
    
    #TODO recalculate the sn based on neigbouring centroids
    #replace NA values in snValues vector by the minimum from peak picking which is 5
    #snValues[is.na(snValues)] <- 3
    
   
    
    #collect centroids
    if (TRUE %in% base::is.na(snValues)) {
      
      set <- base::as.data.frame(xcms::groups(object@xcmsSet[,sampleidxs]))
      
      cent <- xcms::filterMsLevel(xcms::filterFile(featData, file = sampleidxs), msLevel. = 1)
      cent <- base::as.data.frame(methods::as(cent, "data.frame"))
      cent$rt <- base::as.numeric(cent$rt)
      cent$mz <- base::as.numeric(cent$mz)
      cent$i <- base::as.numeric(cent$i)
      base::colnames(cent) <- base::c("file", "rt", "mz", "into")
      
      #Enable full parallel processing
      snow <- BiocParallel::registered("SnowParam")
      if (snow$workers < parallel::detectCores())
      {
        snow <- BiocParallel::SnowParam(workers = parallel::detectCores(), type = "SOCK", exportglobals = FALSE)
        BiocParallel::register(snow, default = TRUE)
      }
      
      snValues <- base::unlist(BiocParallel::bplapply(X = 1:base::length(snValues), cent=cent, set=set, intValues=intValues, snValues=snValues,
                                         function(sn, cent, set, intValues, snValues){
        if (base::is.na(snValues[sn])) {
          temp_cent <- cent[cent$mz >= set$mzmed[sn]-(5/1E6*set$mzmed[sn]) &
                            cent$mz <= set$mzmed[sn]+(5/1E6*set$mzmed[sn]) &
                            cent$rt >= set$rtmin[sn]-240 & cent$rt <= set$rtmax[sn]+240, ]
          temp_cent <- dplyr::filter(temp_cent, !dplyr::between(temp_cent$rt, set$rtmin[sn], set$rtmax[sn]))
          temp_sn <-   intValues[sn]/base::round(stats::quantile(temp_cent$into, probs = base::seq(0,1,0.25))[2], digits = 0)
          base::names(temp_sn) <- base::names(intValues[sn])
          sn <- temp_sn
        } else {
          sn <- snValues[sn]
        }
      }, BPPARAM = BiocParallel::bpparam("SnowParam")))
    }
    
    snValues <- base::round(snValues, digits = 1)
    
    snValues[snValues < 3] <- 0
    intValues[snValues < 3] <- 0
    
    #last check for NA, if present remove both int and sn values
    snValues[base::is.na(snValues)] <- 0
    intValues[base::is.na(snValues)] <- 0
    
  }else{
    ## one sample case
    intValues <- object@groupInfo[, intval, drop=FALSE]
    snValues  <- object@groupInfo[, "sn", drop=FALSE]
    snValues[base::is.na(snValues)] <- 3
    snValues[base::is.na(intValues)] <- 0
    intValues[base::is.na(intValues)] <- 0
  }
  
  ## isotope matrix
  isoMatrix <- base::matrix(ncol=4, nrow=0)
  base::colnames(isoMatrix) <- base::c("mpeak", "isopeak", "iso", "charge")
  
  ################################################################################################
  ## find isotopes
  base::cat("Run isotope peak annotation\n % finished: ")
  
  ## look for isotopes in every pseudospectrum
  numberOfPS <- base::length(object@pspectra)
  isotopeClusterCounter <- 0
  
  #testing code
  # test <- CAMERA::getPeaklist(xA)
  # psIdx <- 8
  
  for(psIdx in 1:numberOfPS){ #numberOfPS
    ## progress
    if(numberOfPS >= 10)      if((psIdx %% base::round(numberOfPS / 10)) == 0)   base::cat((psIdx / base::round(numberOfPS / 10) * 10), ' ')
    
    ## get peak indizes for psIdx-th pseudospectrum
    peakIndeces <- object@pspectra[[psIdx]]
    if(base::length(peakIndeces) <= 1)
      ## Pseudospectrum has only one peak
      next
    
    ## calculate isotopes
    isoMatrixForPS.list <- ntsIUTA::findIsotopesForPSAltered(
      peakIndeces, mzValues[peakIndeces], intValues[peakIndeces], snValues[peakIndeces],
      maxcharge=maxcharge, devppm=devppm, mzabs=mzabs, validateIsotopePatterns=validateIsotopePatterns
    )
    
    if(base::length(isoMatrixForPS.list) == 0)
      ## no isotope cluster found
      next
    
    ## add isotope cluster to isoMatrix
    for(isotopeCluster in 1:base::length(isoMatrixForPS.list)){
      isoClusterMatrix <- isoMatrixForPS.list[[isotopeCluster]]
      numberOfClusterPeaks <- base::nrow(isoClusterMatrix)
      isotopeClusterCounter <- isotopeClusterCounter + 1
      
      ## assign cluster index and box
      isoCluster2 <- base::cbind(
        isoClusterMatrix[, "peak index"], 
        base::rep(x = isotopeClusterCounter, times = numberOfClusterPeaks), 
        isoClusterMatrix[, "isotope"], 
        isoClusterMatrix[, "charge"]
      )
      isoMatrix <- base::rbind(isoMatrix, isoCluster2)
    }
  }
  
  ################################################################################################
  ## create isotope matrix for peak list
  isotopeMatrix <- base::vector(mode="list", length=numberOfPeaks)
  
  for(peakIdx in 1:numberOfPeaks){
    row <- base::match(peakIdx, isoMatrix[, "mpeak"])
    if(base::is.na(row))
      ## no isotope found
      next
    
    ## register isotope data
    isotopeMatrix[[peakIdx]] <- base::list(
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
  base::cat("\nNumber of isotopes:", base::nrow(object@isoID), "\n")
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
#' @examples
#' 
#' 
#' 
findIsotopesForPSAltered <- function(peakIndeces = peakIndeces, mzValues = mzValues,
                                     intValues = intValues,
                                     snValues = snValues,
                                     maxcharge = maxcharge,
                                     devppm = devppm, mzabs = mzabs,
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
  
  ## matrix with all important informationen
  #snValues[is.na(snValues)] <- 1
  #snValues[snValues==0] <- 1
  noiseEstimate <- intValues / snValues
  noiseEstimate[base::is.nan(noiseEstimate)] <- 0
  intensityMin  <- intValues - noiseEstimate
  intensityMax  <- intValues + noiseEstimate
  
  spectrum <- base::cbind(peakIndeces, mzValues, intValues, intensityMin, intensityMax)
  base::colnames(spectrum) <- c("peak index", "mz", "intensity", "intensityMin", "intensityMax")
  
  ## order peaks by m/z
  spectrum <- spectrum[base::order(spectrum[, "mz"]), ]
  numberOfPeaksHere <- base::nrow(spectrum)
  
  if(numberOfPeaksHere <= 1)
    return(base::list())
  
  ## calculate allowed m/z errors for the given masses; at least mzabs
  mzErrors <- devppm * spectrum[, "mz"]
  mzErrors[mzErrors < mzabs] <- mzabs
  
  ###################################################################################################
  ## compute m/z difference and possible isotopic connection for every peak pair in pseudospectrum
  isotopeDifference <- 1.0033548378
  expectedDistances <- isotopeDifference / 1:maxcharge
  hitsForCharge_c_p1_p2 <- base::array(
    dim = base::c(maxcharge, numberOfPeaksHere - 1, numberOfPeaksHere - 1),
    dimnames = base::list(base::c("charge", "peak1", "peak2"))
  )
  
  for(peakIdx in 1:(numberOfPeaksHere - 1)){
    ## create distance matrix
    mzDiff <- spectrum[(peakIdx + 1):numberOfPeaksHere, "mz"] - spectrum[peakIdx, "mz"]
    ## compare distances to expected distances
    mzDiffExp <- base::outer(-expectedDistances, mzDiff, FUN="+")
    ## create hit matrix
    hits <- base::abs(mzDiffExp) <= mzErrors[[peakIdx]]
    ## add hit matrix to big picture
    for(chargeIdx in 1:maxcharge)
      hitsForCharge_c_p1_p2[chargeIdx, peakIdx, 1:(numberOfPeaksHere - peakIdx)] <- hits[chargeIdx, ]
  }#end for peakIdx
  
  ###################################################################################################
  ## find and select isotope chains
  goOnSearching <- TRUE
  peakIsProcessed <- base::rep(x = FALSE, times = numberOfPeaksHere)
  resultChains <- base::list()
  chainCharges <- base::list()
  
  while(goOnSearching){
    candidateChains <- base::list()
    candidateCharge <- base::list()
    
    ## check all charges and potential starting peaks to built all possible isotope chains
    for(chargeIdx in 1:maxcharge){
      ## get potential monoisotopic peaks to built all possible isotope chains
      potentialStartPeaks <- base::which(!peakIsProcessed[1:(numberOfPeaksHere - 1)])
      for(peakIdx in potentialStartPeaks){
        ## start new chain with this peak
        peakIdx1 <- peakIdx
        candidateChain <- base::c(peakIdx1)
        
        ## assemble matching m/z distances to a chain of peaks
        while(
          ## peak is not in any prior chain
          (!peakIsProcessed[peakIdx1]) && 
          ## peak is not the last peak in m/z dimension in the current pseudo spectrum
          (peakIdx1 <= (numberOfPeaksHere - 1)) && 
          ## there are peaks in the right m/z distance
          base::any(hitsForCharge_c_p1_p2[chargeIdx, peakIdx1, ], na.rm = TRUE)){
          ## get matching peak index
          matchingIdx <- base::which(hitsForCharge_c_p1_p2[chargeIdx, peakIdx1, ])
          
          peakIdx2 <- matchingIdx + peakIdx1
          
          if(base::length(peakIdx2) > 1){
            ## more than one peak is candidate as successor: take the best matching one regarding m/z
            
            ## get mass of current peak and successor peaks
            mass1   <- spectrum[peakIdx1, "mz"]
            masses2 <- spectrum[peakIdx2, "mz"]
            
            ## take peak with minimum deviation from the expected value
            masses2 <- masses2 - mass1 - expectedDistances[chargeIdx]
            masses2 <- base::abs(masses2)
            tempIdx <- base::which.min(masses2)
            peakIdx2 <- peakIdx2[[tempIdx]]
          }
          
          if(peakIsProcessed[peakIdx2])
            ## next peak is already part of a isotope chain --> do not elongate chain
            break
          
          ## elongate chain
          candidateChain <- c(candidateChain, peakIdx2)
          ## set peak for next iteration
          peakIdx1 <- peakIdx2
        }#end of while loop for current peak chain
        if(base::length(candidateChain) == 1)
          ## a single peak
          next
        
        ## add candidate chain
        candidateChains[[base::length(candidateChains) + 1]] <- candidateChain
        candidateCharge[[base::length(candidateCharge) + 1]] <- chargeIdx
      }#end for peakIdx
    }#end for chargeIdx
    
    ## select the longest chain
    if(base::length(candidateChains) == 0){
      ## no more chains left -> stop searching
      goOnSearching <- FALSE
    } else {
      ## select the longest chain
      maxChainIdx <- base::which.max(base::sapply(candidateChains, function(x) base::length(x)))
      maxChain    <- candidateChains[[maxChainIdx]]
      maxCharge   <- candidateCharge[[maxChainIdx]]
      ## add chain
      resultChains[[base::length(resultChains) + 1]] <- maxChain
      chainCharges[[base::length(chainCharges) + 1]] <- maxCharge
      ## mark comprised peaks
      peakIsProcessed[maxChain] <- TRUE
    }
  }#end of while loop for peak chains
  
  if(base::length(resultChains) == 0)
    ## no isotope cluster found
    return(base::list())
  
  ###################################################################################################
  ## validate chains
  validatedResultChains <- base::list()
  validatedChainCharges <- base::list()
  ## for validation
  cpObj <- CAMERA::compoundQuantiles(compoundLibrary = "kegg")
  maximumIsotopeNumber <- base::max(cpObj@isotopeSet)
  
  ## available quantiles:
  ## 0.000005 0.999995 0.000010 0.999990 0.000050 0.999950 0.000100 0.999900 0.000500 0.999500 0.001000 0.999000
  ## 0.005000 0.995000 0.010000 0.990000 0.025000 0.975000 0.050000 0.950000 0.100000 0.900000 0.500000
  quantileLow  <- 0.00500
  quantileHigh <- 0.99500
  #quantileLow  <- 0.01
  #quantileHigh <- 0.99
  
  for(chainIdx in base::seq_len(length.out = base::length(resultChains))){
    ## get data
    charge  <- chainCharges[[chainIdx]]
    chain   <- resultChains[[chainIdx]]
    numberOfIsotopes <- base::length(chain)
    
    mass <- spectrum[[chain[[1]], "mz"]] * charge
    compoundMassInRange <- mass < cpObj@maxCompoundMass #the maximum mass allowed, default is 1000
    
    if(validateIsotopePatterns & compoundMassInRange){
      ## validate and decompose chain to segments
      monoisotopicPeakIdx <- 1
      isotopePatternMembership <- base::integer(length = base::length(chain))
      isotopePatternLabel <- 1
      isotopePatternMembership[[1]] <- isotopePatternLabel
      for(peakIdx in 2:base::length(chain)){
        ## check proportion
        #proportion <- spectrum[chain[[monoisotopicPeakIdx]], "intensity"] / spectrum[chain[[peakIdx]], "intensity"]
        
        
        proportionObservedMin <- spectrum[[chain[[monoisotopicPeakIdx]], "intensityMin"]] / spectrum[[chain[[peakIdx]], "intensityMax"]]
        proportionObservedMax <- spectrum[[chain[[monoisotopicPeakIdx]], "intensityMax"]] / spectrum[[chain[[peakIdx]], "intensityMin"]]
        
        isotopeNumber <- peakIdx - monoisotopicPeakIdx
        
        if(isotopeNumber > maximumIsotopeNumber){
          ## pattern is too long to be checkable
          proportionExpectedMin <- 0
          proportionExpectedMax <- Inf
        } else {
          ## fetch expected proportion interval
          proportionExpectedMin <- CAMERA::getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = isotopeNumber, mass = mass, quantile = quantileLow)
          proportionExpectedMax <- CAMERA::getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = isotopeNumber, mass = mass, quantile = quantileHigh)
        }
        
        centerObserved = (proportionObservedMin + proportionObservedMax) / 2;
        centerObserved = base::ifelse(base::is.nan(centerObserved),0,centerObserved);
        centerExpected  = (proportionExpectedMin + proportionExpectedMax) / 2;
        radiusObserved  = (proportionObservedMax - proportionObservedMin) / 2;
        radiusObserved  = base::ifelse(base::is.nan(radiusObserved),0,radiusObserved)
        radiusExpected  = (proportionExpectedMax - proportionExpectedMin) / 2;
        isotopeProportionFits <- base::abs(centerObserved - centerExpected) <= (radiusObserved + radiusExpected)
        
        if (isotopeProportionFits != TRUE & isotopeProportionFits != FALSE | base::is.na(isotopeProportionFits)) isotopeProportionFits <- FALSE
        
        if(!isotopeProportionFits){
          ## isotope proportion does not fit --> start new pattern
          isotopePatternLabel <- isotopePatternLabel + 1
          monoisotopicPeakIdx <- peakIdx
        }
        ## assign pattern membership
        isotopePatternMembership[[peakIdx]] <- isotopePatternLabel
      }## end for peakIdx
      
      # if(FALSE & base::length(base::unique(isotopePatternMembership)) > 1){
      #   if(exists("tableMM48"))
      #     for(labelIdx in unique(isotopePatternMembership)){
      #       chainPos <- which(isotopePatternMembership == labelIdx)[[1]]
      #       massHere <- spectrum[[chain[[chainPos]], "mz"]] * charge
      #       mm48 <- any(abs(tableMM48$Isotope.peak.0.Exact.mass - massHere) <= max(devppm * massHere, mzabs))
      #       cat(paste(", ", mm48))
      #     }
      #   cat(" - ")
      #   cat(paste(
      #     chainIdx, 
      #     numberOfIsotopes, isotopePatternLabel, mass, 
      #     "[", paste(isotopePatternMembership, collapse = ";"), "]", 
      #     "[", paste(spectrum[chain, "intensity"], collapse = ";"), "]",
      #     ifelse(test = isotopePatternLabel > 1, yes = 
      #              paste(spectrum[chain[[1]], "intensity"] / spectrum[chain[[(which(diff(isotopePatternMembership) == 1)[[1]]+1)]], "intensity"], "not in",
      #                    "[", CAMERA::getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = which(diff(isotopePatternMembership) == 1)[[1]], mass = mass, quantile = quantileLow),
      #                    CAMERA::getIsotopeProportion(object = cpObj, isotope1 = 0, isotope2 = which(diff(isotopePatternMembership) == 1)[[1]], mass = mass, quantile = quantileHigh), "]"), no = "n/a")
      #   ), "\n")
      # }
      
      ## extract chain segments
      chainSegments <- base::list()
      for(chainSegmentLabel in 1:isotopePatternLabel){
        chainSegment <- chain[base::which(isotopePatternMembership == chainSegmentLabel)]
        if(base::length(chainSegment) <= 1)
          next
        
        chainSegments[[chainSegmentLabel]] <- chainSegment
      }## end for chainSegment
    } else {
      chainSegments <- base::list()
      chainSegments[[1]] <- chain
    }
    
    ## box chain segments
    for(chainSegment in chainSegments){
      validatedResultChains[[base::length(validatedResultChains) + 1]] <- chainSegment
      validatedChainCharges[[base::length(validatedChainCharges) + 1]] <- charge
    }## end for chainSegment
  }## end for chain
  
  ###################################################################################################
  ## assemble list of isotope clusters as result
  isoMatrixForPS.list <- base::list()
  for(chainIdx in base::seq_len(length.out = base::length(validatedResultChains))){
    ## get data
    chain  <- validatedResultChains[[chainIdx]]
    charge <- validatedChainCharges[[chainIdx]]
    
    ## create and fill isotope cluster
    numberOfIsotopes <- base::length(chain)
    
    isoClusterMatrix <- base::matrix(nrow = numberOfIsotopes, ncol = 3)
    base::colnames(isoClusterMatrix) <- c("peak index", "isotope", "charge")
    for(isoIdx in 1:numberOfIsotopes){
      ## translate local peak index to global peak index
      peakIdx <- spectrum[chain[[isoIdx]], "peak index"]
      ## add isotope entry
      isoClusterMatrix[isoIdx, ] <- c(peakIdx, isoIdx, charge)
    }
    ## add isotope cluster to result list
    isoMatrixForPS.list[[base::length(isoMatrixForPS.list) + 1]] <- isoClusterMatrix
  }## end for chainSegment
  
  return(isoMatrixForPS.list)
}
