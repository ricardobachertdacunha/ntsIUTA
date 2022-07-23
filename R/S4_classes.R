

### ntsSettings -----

#' @title ntsSettings
#'
#' @description A set of settings for a data processing step of the basic workflow.
#' The \linkS4class{ntsSettings} object contains the algorithm to be used and the list of respective settings.
#'
#' @slot algorithm A character string with the name of the algorithm to be used.
#' @slot settings A list of settings dependent on the algorithm used.
#'
#' @return A \linkS4class{ntsSettings} object to be used for
#' the basic workflow of \pkg{ntsIUTA}.
#'
#' @export
#'
setClass("ntsSettings",
  slots = c(
    algorithm = "character",
    settings = "list"
  ),
  prototype = list(
    algorithm = NA_character_,
    settings = list()
  )
)


### ntsParameters -----

#' @title ntsParameters
#'
#' @description A S4 class object to store the parameter settings for
#' each data processing step of the basic workflow.
#' Each parameter set is stored as \linkS4class{ntsSettings} object.
#' See \code{?"ntsSettings-class"}.
#'
#' @slot picking Settings for peak picking.
#' @slot grouping Settings for grouping and alignment within replicates.
#' @slot filling Settings for filling peaks with recursive integration within replicates.
#' @slot annotation	Settings for annotation of isotopes and adducts.
#' @slot fragments Settings for extracting MS2 data.
#' @slot unification Settings for unifying features across replicates.
#'
#' @return A \linkS4class{ntsParameters} object to be used for the basic workflow of \pkg{ntsIUTA}.
#'
#' @export
#'
setClass("ntsParameters",
  slots = c(
    picking = "ntsSettings",
    grouping = "ntsSettings",
    filling = "ntsSettings",
    annotation = "ntsSettings",
    fragments = "ntsSettings",
    unification = "ntsSettings"
  ),
  prototype = list(
    picking = new("ntsSettings"),
    grouping = new("ntsSettings"),
    filling = new("ntsSettings"),
    annotation = new("ntsSettings"),
    fragments = new("ntsSettings"),
    unification = new("ntsSettings")
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


### QC -----

#' @title qcData
#'
#' @slot samples A data.frame listing the QC samples.
#' @slot scans A list with MS scans \link[data.table]{data.table} objects per sample.
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
#' @importFrom purrr quietly
#'
setClass("qcData",
  slots = c(
    title = "character",
    description = "character",
    date = "Date",
    polarity = "character",
    samples = "data.table",
    scans = "list",
    suspects = "suspectList",
    rtWindow = "numeric",
    ppm = "numeric",
    pat = "workflowStep",
    comp = "components",
    ms2 = "MSPeakLists",
    peaks = "data.table",
    features = "data.table",
    unified = "data.table",
    results = "data.table"
  ),
  prototype = list(
    title = NA_character_,
    description = NA_character_,
    date = date,
    polarity = NA_character_,
    samples = data.table::data.table(
      file = character(),
      sample = character(),
      replicate = character(),
      blank = character(),
      scans = numeric(),
      centroided = logical(),
      msLevels = character(),
      rtStart = numeric(),
      rtEnd = numeric(),
      mzLow = numeric(),
      mzHigh = numeric(),
      CE = character(),
      method = character()
    ),
    scans = list(),
    suspects = new("suspectList"),
    rtWindow = 30,
    ppm = 15,
    pat = new("featuresSIRIUS"),
    comp = new("componentsCliqueMS", fGroups = new("featureGroupsSIRIUS", features = new("featuresSIRIUS"))),
    ms2 = purrr::quietly(new)("MSPeakLists", algorithm = "mzr")$result,
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
#' The data within the slots is stored mostly as \link[data.table]{data.table} objects for a more
#' flexible and comprehensive processing and assembly of workflows. Only the slot workflows is a 
#' list composed of \link[data.table]{data.table} objects for each executed workflow step.
#'
#' @slot project A table with five columns:
#' \enumerate{
#'  \item title: a character string with the project title;
#'  \item date: the date of the project;
#'  \item path: a character string with the project path;
#'  \item polarity: character string with associated blank replicate;
#'  \item description: a character string with the project description.
#' }
#' @slot samples A table with six columns:
#' \enumerate{
#'  \item file: character string with the file path;
#'  \item sample: character string with the file name;
#'  \item replicate: character string with the assigned sample replicate group;
#'  \item blank: character string with associated blank replicate;
#'  \item datatype: character string with type of data (i.e., centroided, profile or chromatograms);
#'  \item size: a numeric value with the total number of scans/chromatograms;
#'  \item msLevels: a character string with the existing MS levels;
#'  \item timeStamp: a character string with the start time stamp of thesample/file;
#'  \item method: a character string with the method name or path used to acquire the data.
#' }
#' @slot scans A table with the MS information for each scan in each file/sample.
#' @slot ms1 A table with columns () with MS1 data for each file/sample.
#' @slot ms2 A table with MS2 data for each file/sample.
#' @slot chromsInfo A table with chromatographic info for each file/sample.
#' @slot chroms A table with chromatograms in each file/sample.
#' @slot metadata A table with the same number of rows
#' as the number of \code{samples} containing metadata for each sample added as extra columns.
#' @slot parameters A \linkS4class{ntsParameters} object containing process parameters for the basic workflow.
#' See \code{?"ntsParameters-class"} for more information.
#' @slot controlSamples A \linkS4class{qcData} object used for quality control. More information in \code{?"qcData-class"}.
#' @slot controlResults A \linkS4class{qcData} object used for quality control. More information in \code{?"qcData-class"}.
#' @slot pat An \linkS4class{workflowStep} object from the \pkg{patRoon} package.
#' @slot peaks A table with peaks for each sample.
#' @slot features A table with features (i.e., grouped peaks across samples in the project).
#' @slot unified A table with unified features
#' (i.e., features collapsed by their isotopes and adducts
#' as well as respective features acquired with a different ionization).
#' @slot filters A list of applied filters to the features.
#' @slot removed A table with the removed unified features.
#' @slot workflows A list of objects inherent of downstream data processing steps, such as
#' suspect screening, track of transformations and others.
#'
#' @note The slot \code{pat} contains the native \linkS4class{features} or
#' \linkS4class{featureGroups} object which can be used
#' within native functions of the \pkg{patRoon} package.
#'
#' @references \insertRef{Helmus2021}{ntsIUTA}
#'
#' @export
#'
#' @importFrom data.table data.table
#' @importFrom purrr quietly
#'
setClass("ntsData",
  slots = c(
    project = "data.table",
    samples = "data.table",
    scans = "data.table",
    ms1 = "data.table",
    ms2 = "data.table",
    chromsInfo = "data.table",
    chroms = "data.table",
    metadata = "data.table",
    parameters = "ntsParameters",
    controlSamples = "data.table",
    controlResults = "data.table",
    peaks = "data.table",
    features = "data.table",
    unified = "data.table",
    removed = "data.table",
    pat = "workflowStep",
    patComp = "components",
    patMS2 = "MSPeakLists",
    workflows = "list"
  ),
  prototype = list(
    project = data.table::data.table(
      title = NA_character_,
      date = Sys.Date(),
      path = NA_character_,
      polarity = NA_character_,
      description = NA_character_
    ),
    samples = data.table::data.table(
      file = character(),
      sample = character(),
      replicate = character(),
      blank = character(),
      datatype = character(),
      size = numeric(),
      msLevels = character(),
      timeStamp = character(),
      method = character()
    ),
    scans = data.table::data.table(),
    ms1 = data.table::data.table(
      sample = character(),
      rt = numeric(),
      mz = numeric(),
      intensity = numeric()
    ),
    ms2 = data.table::data.table(
      sample = character(),
      rt = numeric(),
      mz = numeric(),
      intensity = numeric(),
      precMZ = numeric(),
      voltage = numeric()
    ),
    chromsInfo = data.table::data.table(),
    chroms = data.table::data.table(),
    metadata = data.table::data.table(
      sample = character(),
      replicate = character()
    ),
    parameters = new("ntsParameters"),
    controlSamples = data.table::data.table(
      file = character(),
      sample = character(),
      replicate = character(),
      blank = character(),
      datatype = character(),
      size = numeric(),
      msLevels = character(),
      timeStamp = character(),
      method = character()
    ),
    controlResults = data.table::data.table(),
    peaks = data.table::data.table(),
    features = data.table::data.table(),
    unified = data.table::data.table(),
    removed = data.table::data.table(),
    pat = new("featuresSIRIUS"),
    patComp = new("componentsCliqueMS", fGroups = new("featureGroupsSIRIUS", features = new("featuresSIRIUS"))),
    patMS2 = purrr::quietly(new)("MSPeakLists", algorithm = "mzr")$result,
    workflows = list()
  )
)
