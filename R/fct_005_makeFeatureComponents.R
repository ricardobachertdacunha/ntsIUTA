


#' @title makeFeatureComponents
#' @description Uses the \pkg{CAMERA} package to find isotopes by grouping features over the retention time
#' and extracted ion chromotogram (EIC) for the given deviations. The isotopes are found for each replicate group.
#'
#' @param featData An \linkS4class{XCMSnExp} object containing features.
#' @param polarity The polarity of the replicate groups. Possible values are \code{positive} or \code{negative}.
#' @param sigma The multiplier of the standard deviation for grouping features by retention time.
#' @param perfwhm Percentage of the width of the FWHM.
#' @param cor_eic_th Minimum correlation index (from 0 to 1) for the EICs of features within the same sample.  
#' @param cor_exp_th Minimum correlation index (from 0 to 1) for the EICs of features across samples within the replicate group.
#' @param pval p-value threshold for testing correlation of significance.
#' @param calcCaS Logical, set to \code{TRUE} to calculate correlation accross samples. Deafault is \code{TRUE}.
#' @param calcIso Logical, set to \code{TRUE} to include isotope detection informationen for graph clustering. Deafault is \code{FALSE}.
#' @param ppmIsotopes The expected mass deviation to find isotopes.
#' @param mzabs The expected deviation of the \emph{m/z} to find isotopes.
#' @param noise The extimated intensity threshold for the noise level used for find isotopes.
#' @param validateIsotopePatterns Logical, set to \code{TRUE} for validating the annoatated isotopes with the kegg database. 
#' @param searchAdducts Logical, set to \code{TRUE} to screen for adducts after finding isotopes.
#' @param ppmAdducts The expected mass deviation to find adducts.
#' @param extendedList Logical, set to \code{TRUE} to use the extended list of adducts. The default (\code{FALSE} uses a shorter list.)
#' @param excludeBlanks Set to \code{TRUE} to not screen blank samples for isotopes and adducts.
#' This option is intended for saving processing time as the replicate blank groupes are expected to be removed before final feature list creation.
#' @param blankGroups A character vector with the name of the blank replicate groups in the given \code{featData} object.
#' For simplification use \code{base::unique(sampleInfo$blank)}.
#' @param save Logical, set to \code{TRUE} to save the generated \code{list} of \linkS4class{xsAnnotate} objects in the disk.
#' @param projPath The \code{projPath} directory as defined in the \code{setup} object.
#'
#' @return A \code{list} with an \linkS4class{xsAnnotate} object per replicate group as defined in the given \linkS4class{XCMSnExp} object.
#' 
#' @references
#' \insertRef{CAMERA}{ntsIUTA}
#' \insertRef{xcms1}{ntsIUTA}
#' \insertRef{xcms2}{ntsIUTA}
#' \insertRef{xcms3}{ntsIUTA}
#' 
#' @export
#' 
#' @importFrom methods as
#' @importFrom xcms filterMsLevel sampnames
#' @importFrom Biobase pData
#' @importFrom CAMERA xsAnnotate groupFWHM groupCorr findAdducts
#' @importFrom utils txtProgressBar setTxtProgressBar read.table
#'
#' @examples
#' 
#' 
#' 
makeFeatureComponents <- function(featData = featData, polarity = "positive",
                                  sigma = 5, perfwhm = 0.5, cor_eic_th = 0.85,
                                  cor_exp_th = 0.75, pval = 0.05,
                                  calcCaS = TRUE, calcIso = TRUE,
                                  ppmIsotopes = 50, mzabs = 0.01, noise = 350,
                                  validateIsotopePatterns = TRUE,
                                  searchAdducts = TRUE, ppmAdducts = 5, extendedList = FALSE,
                                  excludeBlanks = TRUE, blankGroups = "Blank",
                                  save = TRUE, projPath = setup$projPath) {
  
  #Examples
  # setup <- ntsIUTA::makeSetup(projPath = system.file(package = "ntsIUTA", dir = "extdata"), save = FALSE)
  # setup$sampleInfo[1:3,"group"] <- "Sample"
  # rawDataExample <- ntsIUTA::importRawData(setup$sampleInfo[1:3,], save = FALSE, centroidedData = TRUE)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = param, save = FALSE)
  # paramList <- ntsIUTA::paramListExample
  # instParam <- ntsIUTA::getInstParam(paramList)
  # peaksDataExample <- ntsIUTA::peakPicking(rawDataExample, param = instParam$PP, save = FALSE)
  # featDataExample <- ntsIUTA::makeFeatures(peaksDataExample, paramPreGrouping = instParam$preGrouping, paramAlignment = instParam$alignment, paramGrouping = instParam$grouping, save = FALSE)
  # featCompExample <- ntsIUTA::makeFeatureComponents(featData = featDataExample)
  
  
  # featData = featData
  # polarity = "positive"
  # sigma = 5
  # perfwhm = 0.45
  # cor_eic_th = 0.85
  # cor_exp_th = 0.85
  # pval = 0.05
  # ppmIsotopes = 50
  # mzabs = 0.01
  # noise = 350
  # searchAdducts = TRUE
  # ppmAdducts = 5
  # extendedList = FALSE
  # excludeBlanks = FALSE
  # calcCaS = TRUE
  # calcIso = TRUE
  
  #convert featData to xcmsSet class
  xSet <- methods::as(xcms::filterMsLevel(featData,  msLevel. = 1), "xcmsSet")
  xcms::sampnames(xSet) <- Biobase::pData(featData)$sample_name
  xcms::sampclass(xSet) <- Biobase::pData(featData)$sample_group
  #xSet <- xcms::fillPeaks(xSet) #not used as filling was already performed
  
  #Replicate sample group divider
  groups <- base::unique(featData$sample_group)
  
  
  if (excludeBlanks) groups <- groups[groups != blankGroups]
  
  #Holder for isotopes from each replicate group as list() and indices for which replicate group
  featComp <- base::list()
  
  xA <- CAMERA::xsAnnotate(xs = xSet, sample = c(1:base::length(featData$sample_group)),
                           polarity = polarity)
  
  xA <- CAMERA::groupFWHM(xA, sigma = sigma, #the multiplier of the standard deviation
                          perfwhm = perfwhm, #percentage of the width of the FWHM
                          intval = "maxo")
  
  xA <- CAMERA::groupCorr(xA, cor_eic_th = cor_eic_th, # Correlation threshold for EIC correlation
                          pval = pval, # p-value threshold for testing correlation of signiﬁcance
                          graphMethod = "hcs",
                          calcIso = calcIso, calcCiS = TRUE, calcCaS = calcCaS, psg_list = NULL, xraw = NULL,
                          cor_exp_th = cor_exp_th, # Threshold for intensity correlations across samples
                          intval = "maxo") #suppressWarnings()
  
  # test <- CAMERA::getPeaklist(xA)
  # test <- dplyr::filter(test, pcgroup == 8)
  # #testx <- as.data.frame(xcms::featureDefinitions(featData))
  # #testX <-base::row.names(xcms::featureDefinitions(featData))
  
  
  #Make peaksData by replicate groups as defined in the setup experiment
  pb <- utils::txtProgressBar(min = 0, max = 100, initial = 0, char = "=", width = 80, style = 3)
  
  for (rgidx in 1:base::length(groups)) {
    
    utils::setTxtProgressBar(pb, ((rgidx/base::length(groups))*100))
    
    sampleidxs <- base::which(featData$sample_group == groups[rgidx])
    
    xA_temp <- ntsIUTA::FindIsotopesWithValidationAltered(object = xA,
                                                          featData = featData,
                                                          sampleidxs = sampleidxs,
                                                          ppm = ppmIsotopes,
                                                          mzabs = mzabs,
                                                          noise = noise,
                                                          maxcharge = 3, #maxcharge set to match small molecules and lipids
                                                          intval = "maxo",
                                                          validateIsotopePatterns = validateIsotopePatterns)
    
    if (searchAdducts)
    {
      if (extendedList) {
        rules_pos <- base::system.file('rules/extended_adducts_pos.csv', package = "CAMERA")
        rules_neg <- base::system.file('rules/extended_adducts_neg.csv', package = "CAMERA")
      } else {
        rules_pos <- base::system.file('rules/primary_adducts_pos.csv', package = "CAMERA")
        rules_neg <- base::system.file('rules/primary_adducts_neg.csv', package = "CAMERA")
      }
      
      
      if (polarity == "positive")
      {rules <- utils::read.table(rules_pos, header = TRUE, sep = ",")}
      
      if (polarity == "negative")
      {rules <- utils::read.table(rules_neg, header = TRUE, sep = ",")}
      
      
      xA_temp <- CAMERA::findAdducts(xA_temp,
                                     ppm = ppmAdducts,
                                     mzabs = 0,
                                     multiplier = 2, # highest number(n) of allowed clusterion [nM+ion]
                                     polarity = polarity,
                                     rules = rules,
                                     max_peaks = 100, # If run in parralel mode, this number deﬁnes how much peaks will be calculated in every thread
                                     psg_list = NULL) # Vector of pseudospectra indices. The correlation analysis will be only done for those groups
    }
    
    featComp[[groups[rgidx]]] <- xA_temp
    
  }
  
  if (save)
  {
    rData <- base::paste0(projPath,"\\rData")
    if (!base::dir.exists(rData)) base::dir.create(rData)
    base::saveRDS(featComp, file = base::paste0(rData,"\\featComp.rds"))
  }
  
  return(featComp)
  
}




