

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




### MS2param -----

#' @title MS2param-class
#'
#' @slot maxMSRtWindow .
#' @slot precursorMzWindow .
#' @slot clusterMzWindow .
#' @slot topMost .
#' @slot minIntensityPre .
#' @slot minIntensityPost .
#'
#' @return An \linkS4class{MS2param} object containing parameters for
#' extraction of MS2 spectra of given precursor ions in an \linkS4class{ntsData} object.
#'
#' @export
#'
setClass("MS2param",
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




### paramList -----

#' @title paramList
#'
#' @slot peakPicking A paramSet object for performing peak picking.
#' @slot peakGrouping A paramSet object for grouping peaks across samples.
#' @slot fillMissing A paramSet object for performing recursive peak integration.
#' @slot annotation A paramSet object for annotation of peaks and features in the \linkS4class{ntsData} object
#' @slot MS2 An \linkS4class{MS2param} object with settings for extraction of MS2 data.
#'
#' @return A \linkS4class{paramList} object to be used for
#' the basic workflow of \pkg{ntsIUTA}.
#'
#' @export
#'
setClass("paramList",
  slots = c(
    peakPicking = "paramSet",
    peakGrouping = "paramSet",
    fillMissing = "paramSet",
    annotation = "paramSet",
    MS2 = "MS2param"
  ),
  prototype = list(
    peakPicking = new("paramSet"),
    peakGrouping = new("paramSet"),
    fillMissing = new("paramSet"),
    annotation = new("paramSet"),
    MS2 = new("MS2param")
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
setClass("qcData",
  slots = c(
    samples = "data.frame",
    targets = "suspectList",
    rtWindow = "numeric",
    ppm = "numeric",
    MSnExp = "OnDiskMSnExp",
    patdata = "workflowStep",
    peaks = "data.frame",
    features = "data.frame",
    annotation = "list",
    results = "data.frame"
  ),
  prototype = list(
    samples = data.frame(
      file = character(),
      sample = character(),
      group = character(),
      blank = character()
    ),
    targets = new("suspectList"),
    rtWindow = 30,
    ppm = 15,
    MSnExp = new("OnDiskMSnExp"),
    patdata = new("featuresXCMS3"),
    peaks = data.frame(),
    features = data.frame(),
    annotation = list(),
    results = data.frame()
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
#' @description S4 class object to organize and store project processed data
#' within the \pkg{ntsIUTA} package.
#'
#' @slot title A character string with the project title.
#' @slot info A list to store project details.
#' The first position is for description of the project.
#' @slot path A character string with the project path.
#' as defined with \code{\link{setupProject}}.
#' @slot date The project date.
#' @slot polarity A string character defining the polarity mode of
#' the MS files in the \linkS4class{ntsData} object.
#' Possible values are "positive" or "negative".
#' @slot samples A \code{data.frame} as obtained by \code{\link{setupProject}}.
#' @slot metadata A \code{data.frame} with the same number of rows
#' as the \code{samples} containing metadata to support data interpretation.
#' @slot parameters A list containing process parameters for the workflow steps.
#' @slot QC Slot for assigning the QC sample replicate group and results.
#' Note that the assigned group is not included in the basic workflow
#' but only used for quality control, using the function \code{\link{checkQC}},
#' and for optimization of workflow parameters. A screening list with target substances
#' must be provided for running quality control.
#' @slot MSnExp An \linkS4class{OnDiskMSnExp} as generated by \code{\link{setupProject}}.
#' @slot patdata A \linkS4class{features} or \linkS4class{featureGroups}
#' object derived from the basic NTS workflow
#' (peak picking, alignment and grouping).
#' @slot peaks A list with a data.frame containing the peak information
#' for each sample in the slot \code{samples}.
#' @slot features A data.frame with the features information,
#' including feature quality data and filtering information.
#' @slot annotation Slot for storage of annotation of isotopes and adducts and
#' source fragments. Features with similar chromatographic behaviour
#' are grouped by components. Isotopologues are grouped by \code{isonumber} and
#' are created for each sample replicate group. Similarly, adducts are searched
#' for each sample replicate group and share the same molecular ion (\code{Mion}).
#' The functionality of the package \pkg{CAMERA} was adapted to work better for
#' environmental analyses, where sample replicate groups are, in principle,
#' very different from each other. When specified other functionalities
#' from the package \pkg{patRoon} can be used but the default
#' is the package \pkg{CAMERA} directly. The slot is composed of a list of length 2,
#' where the first entry is a data frame with summarized results for each feature
#' and the second is the raw data from the \pkg{CAMERA} package.
#' @slot IS Slot to store the list of internal standards (IS) and the results
#' from IS control of the samples in the \code{ntsData} object.
#' @slot filters A list of applied filters to the features data.frame.
#' @slot removed A data frame with the removed features.
#' @slot workflows A list of objects inherent of downstream workflows, such as
#' suspect screening.
#'
#' @return An S4 class object named \linkS4class{ntsData}.
#'
#' @note Slot \code{MSnExp} is used to save the \linkS4class{OnDiskMSnExp}
#' object obtained by \pkg{MSnbase} during the \code{\link{setupProject}}.
#' Slot \code{patdata} will contain the object \linkS4class{features} or
#' \linkS4class{featureGroups} generated from the \pkg{patRoon} package.
#' Both objects can be used within native functions of \pkg{MSnbase} and
#' \pkg{patRoon}.
#'
#' @export
#'
setClass("ntsData",
  slots = c(
    title = "character",
    info = "list",
    path = "character",
    date = "Date",
    polarity = "character",
    samples = "data.frame",
    metadata = "data.frame",
    parameters = "paramList",
    QC = "qcData",
    MSnExp = "OnDiskMSnExp",
    patdata = "workflowStep",
    peaks = "data.frame",
    features = "data.frame",
    annotation = "list",
    IS = "isData",
    filters = "list",
    removed = "data.frame",
    workflows = "list"
  ),
  prototype = list(
    title = "New Project",
    info = list(
      description = NA_character_,
      history = list()
    ),
    path = NA_character_,
    date = Sys.Date(),
    polarity = "positive",
    samples = data.frame(
      file = character(),
      sample = character(),
      group = character(),
      blank = character()
    ),
    metadata = data.frame(group = character()),
    parameters = new("paramList"),
    QC = new("qcData"),
    MSnExp = new("OnDiskMSnExp"),
    patdata = new("featuresXCMS3"),
    peaks = data.frame(),
    features = data.frame(),
    annotation = list(
      comp = data.frame(),
      raw = list()
    ),
    IS = new("isData"),
    filters = list(),
    removed = data.frame(),
    workflows = list()
  )
)
