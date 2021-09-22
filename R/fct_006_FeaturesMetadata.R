

calculateMetadata <- function(obj,
                              ID = NULL) {
  
  pat <- obj@patdata
  
  if (!is.null(ID)) pat <- pat[, ID]
  
  #calculates qualities of each peak in each sample
  pat <- calculatePeakQualities(pat,
                                weights = NULL,
                                flatnessFactor = 0.05, #Passed to MetaClean as the flatness.factor argument to calculateJaggedness and calculateModality.
                                avgFunc = mean, # mean additional parameter for handling featureGroups
                                parallel = TRUE)
  
  
  
  
  
}