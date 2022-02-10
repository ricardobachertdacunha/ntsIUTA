

### suspectList -----

#' @title suspectList
#'
#' @slot path The path of the csv file in disk.
#' @slot data The data.frame listing all the suspects of interest.
#' @slot length The number of compounds in the suspect list.
#' @slot comment Optional comment for the suspect list.
#'
#' @return An \linkS4class{suspectList} object to be used within screening
#' workflows of \pkg{ntsIUTA}.
#'
#' @export
#'
setClass("suspectList",
  slots = c(
    path = "character",
    data = "data.frame",
    length = "numeric",
    comment = "character"
  ),
  prototype = list(
    path = character(),
    data = data.frame(),
    length = numeric(),
    comment = character()
  )
)




### paramSet -----

#' @title paramSet
#'
#' @description A set of parameters for a data processing step of the basic workflow.
#' The \linkS4class{paramSet} object contains the algorithm to be used and the list of respective parameters. 
#' 
#' @slot algorithm A character string with the name of the algorithm to be used.
#' @slot param A list of parameters dependent on the algorithm used.
#'
#' @return A \linkS4class{paramSet} object to be used for
#' the basic workflow of \pkg{ntsIUTA}.
#'
#' @export
#'
setClass("paramSet",
  slots = c(
    algorithm = "character",
    param = "list"
  ),
  prototype = list(
    algorithm = NA_character_,
    param = list()
  )
)




### AlteredCameraParam -----

#' @title AlteredCameraParam
#'
#' @slot sigma The multiplier of the standard deviation for grouping features
#' by retention time.
#' @slot perfwhm Percentage of the overlapping FWHM to group features.
#' @slot cor_eic_th Threshold for feature EIC correlation in each sample.
#' @slot cor_exp_th Threshold for intensity correlation across samples.
#' @slot pval p-value threshold for testing correlation of significance.
#' @slot validateIsotopePatterns Logical, set to \code{TRUE} for validating
#' the annotated isotopes with the \emph{kegg} database.
#' @slot ppmIsotopes The expected mass deviation (in ppm) to find isotopes.
#' @slot noise numeric.
#' @slot searchAdducts Logical, set to \code{TRUE} to screen for adducts after finding isotopes.
#' @slot ppmAdducts The expected mass deviation (in ppm) to find adducts.
#' @slot extendedList Logical, set to \code{TRUE} to use the extended list of
#' adducts. The default is \code{FALSE}.
#'
#' @return An \linkS4class{AlteredCameraParam} object containing parameters for
#' annotation of features in an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("AlteredCameraParam",
  slots = c(
    sigma = "numeric",
    perfwhm = "numeric",
    cor_eic_th = "numeric",
    cor_exp_th = "numeric",
    pval = "numeric",
    validateIsotopePatterns = "logical",
    ppmIsotopes = "numeric",
    noise = "numeric",
    searchAdducts = "logical",
    ppmAdducts = "numeric",
    extendedList = "logical"
  )
)




### paramMS2 -----

#' @title paramMS2-class
#'
#' @slot maxMSRtWindow .
#' @slot precursorMzWindow .
#' @slot clusterMzWindow .
#' @slot topMost .
#' @slot minIntensityPre .
#' @slot minIntensityPost .
#'
#' @return An \linkS4class{paramMS2} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("paramMS2",
  slots = c(
    maxMSRtWindow = "numeric",
    precursorMzWindow = "numeric",
    clusterMzWindow = "numeric",
    topMost = "numeric",
    minIntensityPre = "numeric",
    minIntensityPost = "numeric"
  ),
  prototype = list(
    maxMSRtWindow = 10,
    precursorMzWindow = 1.3,
    clusterMzWindow = 0.003,
    topMost = 50,
    minIntensityPre = 10,
    minIntensityPost = 10
  )
)

# TODO Improve the extraction of MS2 data


### paramList -----

#' @title paramList
#'
#' @description A S4 class object to store the parameter sets for
#' each data processing step of the basic workflow.
#' Each parameter set is stored as \linkS4class{paramSet} object.
#' For MS2 data, the parameters are stored in a \linkS4class{paramMS2} object.
#' See \code{?"paramSet-class"} and \code{?"paramMS2-class"} for more information. 
#'
#' @slot peakPicking A paramSet object for performing peak picking in each ms file.
#' @slot peakGrouping A paramSet object for grouping and alignment of peaks withn replicate samples.
#' @slot fillMissing A paramSet object for performing recursive peak integration within replicates.
#' @slot annotateIsotopes A paramSet object for annotation of isotopes within replicates.
#' @slot annotateAdducts A paramSet object for annotation of adducts across unified features.
#' @slot MS2 An \linkS4class{paramMS2} object with settings for extraction of MS2 data.
#'
#' @return A \linkS4class{paramList} object to be used for the basic workflow of \pkg{ntsIUTA}.
#'
#' @export
#'
setClass("paramList",
  slots = c(
    peakPicking = "paramSet",
    peakGrouping = "paramSet",
    fillMissing = "paramSet",
    annotation = "paramSet",
    MS2 = "paramMS2"
  ),
  prototype = list(
    peakPicking = new("paramSet"),
    peakGrouping = new("paramSet"),
    fillMissing = new("paramSet"),
    annotation = new("paramSet"),
    MS2 = new("paramMS2")
  )
)




### QC -----

#' @title qcData
#'
#' @slot samples A data.frame listing the QC samples.
#' @slot targets A \linkS4class{suspectList} object with the QC target compounds.
#' @slot rtWindow The retention time window, in seconds, to look for target compounds.
#' @slot ppm The mass window, in ppm, to look for target compounds.
#' @slot MSnExp An \linkS4class{OnDiskMSnExp} object with QC raw data.
#' @slot patdata A \linkS4class{features} or \linkS4class{featureGroups}
#' object derived from the basic NTS workflow (peak picking, alignment and grouping).
#' @slot peaks A data.frame containing the peak information for each QC sample.
#' @slot features A data.frame with the features information.
#' @slot annotation A list with the annotation for the QC samples.
#' @slot results A data frame with the QC results summary.
#'
#' @return A \linkS4class{qcData} object to store the QC parameters, data and results
#' from an \linkS4class{ntsData} object.
#'
#' @export
#'
#' @importFrom data.table data.table
#'
setClass("qcData",
  slots = c(
    title = "character",
    description = "character",
    date = "Date",
    samples = "data.table",
    targets = "suspectList",
    rtWindow = "numeric",
    ppm = "numeric",
    pat = "list",
    peaks = "data.table",
    features = "data.table",
    unified = "data.table",
    results = "data.table"
  ),
  prototype = list(
    title = NA_character_,
    description = NA_character_,
    date = date,
    samples = data.table::data.table(
      file = character(),
      sample = character(),
      replicate = character(),
      blank = character(),
      polarity = character(),
      method = character()
    ),
    targets = new("suspectList"),
    rtWindow = 30,
    ppm = 15,
    pat = list(),
    peaks = data.table::data.table(),
    features = data.table::data.table(),
    unified = data.table::data.table(),
    results = data.table::data.table()
  )
)

# TODO make methods to get a ntsData out of a qcData for time analysis


### IS -----

#' @title isData
#'
#' @slot targets A \linkS4class{suspectList} object with the QC target compounds.
#' @slot ppm The mass window, in ppm, to look for target compounds.
#' @slot rtWindow The retention time window, in seconds, to look for target compounds.
#' @slot recoveryFrom A character object with the name of the sample replicate group
#' used for recovery evaluation.
#' @slot results A data frame with the QC results summary.
#'
#' @return A \linkS4class{qcData} object to store the QC parameters, data and results
#' from an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("isData",
  slots = c(
    targets = "suspectList",
    ppm = "numeric",
    rtWindow = "numeric",
    recoveryFrom = "character",
    results = "list"
  ),
  prototype = list(
    targets = new("suspectList"),
    ppm = 10,
    rtWindow = 30,
    recoveryFrom = NA_character_,
    results = list()
  )
)




### ntsData -----

#' @title ntsData
#'
#' @description A S4 class object to store project processed data within the \pkg{ntsIUTA} package.
#'
#' @slot title A character string with the project title.
#' @slot description A character string with the project description.
#' @slot date The \link{Date} of the project.
#' @slot path A character string with the project path.
#' @slot samples A data.table with six columns:
#' \enumerate{
#'  \item file: character string with the file path;
#'  \item sample: character string with the file name;
#'  \item replicate: character string with the assigned sample replicate group;
#'  \item blank: character string with associated blank replicate;
#'  \item polarity: the polarity mode of the respective file; possible values are "positive" or "negative");
#'  \item method: a character string with the method used to acquire the data file.
#' }
#' @slot metadata A data.table with the same number of rows
#' as the number of \code{samples}, containing metadata for each sample added as extra columns.
#' @slot parameters A \linkS4class{paramList} object containing process parameters for the basic workflow.
#' See \code{?"paranList-class"} for more information.
#' @slot QC A \linkS4class{qcData} object used for quality control. More information in \code{?"qcData-class"}.
#' @slot pat A length two list with:
#' \enumerate{
#'  \item A list of \linkS4class{features} or \linkS4class{featureGroups} objects
#' for each sample replicate derived from the basic workflow
#' (peak picking, alignment and grouping) using the \pkg{patRoon} package;
#'  \item A \linkS4class{featureGroups} object from the \pkg{patRoon} package with unified features.
#' }
#' @slot peaks A list of data.tables with peaks for each sample.
#' @slot features A list of data.tables with the features for each replicate.
#' @slot unified A data.table with unified features across different replicates.
#' @slot filters A list of applied filters to the unified features data.table.
#' @slot removed A data table with the removed unified features.
#' @slot workflows A list of objects inherent of downstream data processing steps, such as
#' suspect screening.
#'
#' @note The slot \code{pat} contains objects \linkS4class{features} or
#' \linkS4class{featureGroups} which can be used within native functions of the \pkg{patRoon} package.
#'
#' @export
#'
#' @importFrom data.table data.table
#'
setClass("ntsData",
  slots = c(
    title = "character",
    description = "character",
    date = "Date",
    path = "character",
    samples = "data.table",
    metadata = "data.table",
    parameters = "paramList",
    QC = "qcData",
    pat = "list",
    peaks = "list",
    features = "list",
    unified = "data.table",
    filters = "list",
    removed = "data.table",
    workflows = "list"
  ),
  prototype = list(
    title = NA_character_,
    description = NA_character_,
    date = Sys.Date(),
    path = NA_character_,
    samples = data.table::data.table(
      file = character(),
      sample = character(),
      replicate = character(),
      blank = character(),
      polarity = character(),
      method = character()
    ),
    metadata = data.table::data.table(replicate = character()),
    parameters = new("paramList"),
    QC = new("qcData"),
    pat = list(replicates = list(), unified = list()),
    peaks = list(),
    features = list(),
    unified = data.table::data.table(),
    filters = list(),
    removed = data.table::data.table(),
    workflows = list()
  )
)
