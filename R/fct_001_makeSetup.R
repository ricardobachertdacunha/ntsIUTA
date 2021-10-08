

#' @title setupProject
#'
#' @description \code{setupProject} takes a given path to create a project for Non-target screening (NTS).
#' mzML or mzXML files in the path are directly added to the project.
#' Yet, the function \code{\link{addFiles}} can be used to add files or extra files to the project.
#' Setting \code{createRproject} to \code{TRUE} will create and open an R project in the selected or given folder.
#' When \code{convertFiles} is \code{TRUE} the function uses the function \code{\link{mzMLconverter}} to automatically convert
#' specified raw MS files in the path to mzML and add them to the project.
#'
#' @param path The directory for project setup. The default is the \code{utils::choose.dir()} to select or create a folder.
#' Note that the function will look into subfolders for mzML/mzXML files or specified raw files when \code{convertFiles} is set to \code{TRUE}.
#' @param title A character string with the project title.
#' @param description A character string with the project description.
#' @param date The project date.
#' @param polarity A vector specifying the polarity mode.
#' Possible values are \code{positive} or \code{negative} for positive and negative modes, respectively.
#' @param convertFiles Logical, set to \code{TRUE} for non mzML/mzXML files being converted to centroided mzML.
#' The default is \code{FALSE} for no conversion even if files are present.
#' @param convertFrom The name of the file format or vendor format to be found and converted.
#' Possible values are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.
#' @param convertToCentroid Logical, set to \code{TRUE} to convert profile data to centroided when converting raw MS files to mzML.
#' The default and recommended is \code{TRUE}.
#' @param groups A vector with the identifier/name of each sample replicate group. If \code{NULL},
#' a grouping tentative will be made using the name of the MS files found in the given \code{path}.
#' @param save Logical, set to \code{TRUE} to save the \linkS4class{ntsData} object in the \strong{rData} folder.
#' If \strong{rData} folder does not exist in given \code{path}, the folder will be created.
#' @param makeNewProject Logical, set to TRUE to create an R project in the given \code{path} and open a new R session.
#'
#' @return An S4 class object named \linkS4class{ntsData}
#' containing initialisation data for the NTS project.
#'
#' @details The \linkS4class{ntsData} object can be edited using the class methods.
#' See \linkS4class{ntsData} for more information about the slots and methods.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#' @importFrom methods new
#' @importFrom rstudioapi initializeProject openProject
#'
#' @examples
#' path <-  system.file(package = "ntsIUTA", dir = "extdata")
#' dt <- setupProject(path = path, save = FALSE)
#' dt
#'
setupProject <- function(path = utils::choose.dir(getwd(), "Select or create a project folder"),
                         title = NA_character_,
                         description = NA_character_,
                         date = Sys.Date(),
                         polarity = "positive",
                         convertFiles = FALSE,
                         convertFrom = NULL,
                         convertToCentroid = TRUE,
                         groups = NULL,
                         save = TRUE,
                         makeNewProject = FALSE) {

  obj <- new("ntsData")

  if (!is.na(title)) obj@title <- as.character(title)

  if (!is.na(description)) obj@info$description <- as.character(description)

  obj@path <- as.character(path)

  obj@date <- date

  if (testChoice(polarity, c("positive", "negative"))) obj@polarity <- polarity

  if (convertFiles) {

    if (testChoice(convertFrom, c("thermo", "bruker", "agilent", "ab", "waters"))) {

      mzMLconverter(path = path,
                    convertFrom = convertFrom,
                    centroidMethod = ifelse(convertToCentroid == TRUE, "vendor", NULL))

    } else {
      stop("Vendors should be specified for file recognition.
            Possible entries are: thermo, bruker, agilent,
            ab (from AB Sciex) and waters.")
    }
  }

  msFiles <- list.files(path = obj@path,
                        pattern = ".mzML|.mzXML",
                        recursive = TRUE,
                        full.names = TRUE,
                        no.. = TRUE)

  if (length(msFiles) == 0) {
    warning("No mzML or mzXML files were found in selected project path.
            You can use addFiles() to manually add files
            to the ntsData class object.")
  } else {

    obj <- addFiles(newFiles = msFiles,
                    obj = obj,
                    copyFiles = FALSE,
                    groups = groups)

  }

  obj@metadata <- obj@samples[, "sample", drop = FALSE]

  if (makeNewProject) {

    if (save) saveObject(obj = obj)

    #TODO add template for main script file
    sp <- file.path(obj@path, "script.R")
    cat(
"
\n
#Run the following code to load the project setup:\n
library(ntsIUTA)
\n
setup <- readRDS('rData/ntsData.rds')",

    file = sp, sep = "")

    if (!(is.na(Sys.getenv()["RSTUDIO"]))) {
      rstudioapi::initializeProject(obj@path)
      rstudioapi::openProject(obj@path, newSession = TRUE)
    } else {
      setwd(obj@path)
      print("When not using RStudio, it is recommended
          to open the project folder through the IDE
          so that the files and workspace can be loaded. \n
          Run  rstudioapi::navigateToFile('mainscript.R')
          in the new project to open the script.R file.")
      return(obj)
    }

  } else {

    if (save) saveObject(obj = obj)

    return(obj)

  }
}




#' @title addFiles
#' @description Adds mzML or mzXML files to an \linkS4class{ntsData} object.
#'
#' @param newFiles A list of paths for selected mzML or mzXML files.
#' The default is a UI to choose the files.
#' @param obj An \linkS4class{ntsData} object to add the selected files.
#' @param copyFiles Logical. Set to \code{TRUE} to copied the selected files
#' to the \code{path} of the \linkS4class{ntsData} object.
#' \code{FALSE} will add the original location of the files
#' to the \code{samples} table.
#' Note that if duplicate names exist, the process is aborted.
#' @param groups Optional vector with the same length as \code{newFiles}
#' with the identifier/name of each sample replicate group.
#' The default is \code{NULL} which attempts  to group the samples using
#' the name of the MS files. The groups can always be edited using
#' the class method \code{sampleGroups}.
#'
#' @return Returns the \linkS4class{ntsData} object including the added files
#' in the \code{samples} slot.
#'
#' @note When the \linkS4class{ntsData} has files, new files are added
#' below in the \code{samples}. If processed data is already present,
#' the new files will invalidate the structure.
#' Therefore, data should be processed again.
#'
#' @importFrom checkmate assertClass
#' @importFrom dplyr semi_join anti_join mutate
#' @importFrom utils choose.files askYesNo
#' @importFrom tools file_ext
#'
#' @export
#'
addFiles <- function(newFiles = utils::choose.files(),
                     obj = ntsData,
                     copyFiles = FALSE,
                     groups = NULL) {

  assertClass(obj, "ntsData")

  samples <- obj@samples

  if (length(newFiles) < 1) return(obj)

  samples2 <- data.frame(file = character(),
                         sample = character(),
                         group = character(),
                         blank = character())

  samples2[seq_len(length(newFiles)), ] <- NA
  samples2$file <- newFiles
  samples2$sample <- gsub(".mzML|.mzXML", "", basename(newFiles))

  if (is.null(groups)) {
    samples2$group <- gsub(".mzML|.mzXML", "", basename(newFiles))
    samples2 <- mutate(samples2, group = ifelse(grepl("A-r|qc|QC", samples2$sample), "QC", group))
    samples2 <- mutate(samples2, group = ifelse(grepl("Blank|blank|B-r", samples2$sample), "Blank", group))
  } else {
    if (nrow(samples2) == length(groups)) samples2$group <- groups
  }

  oldJoint <- rbind(samples, samples2)
  duplicates <- duplicated(oldJoint$sample)

  if (TRUE %in% duplicates) {
    warning("Process aborted, duplicate files names found! Rename files with the same name.")
    return(obj)
  } else {
    newJoint <- oldJoint
  }

  if (copyFiles) {
    msdir <- paste0(obj@path, "\\msFiles")
    if (!dir.exists(msdir)) dir.create(msdir)
  }

  for (i in seq_len(length(newFiles))) {
    i2 <- which(newJoint$file == newFiles[i])
    ext <- file_ext(newFiles[i])

    if (copyFiles) {
      dir <- msdir
    } else {
      dir <- dirname(newJoint$file[i2])
    }

    newfilepath <- paste(dir, "\\", newJoint$sample[i2], ".", ext, sep = "")

    if (copyFiles) {
      if (!(file.exists(newfilepath))) {
        file.copy(newFiles[i], newfilepath, overwrite = FALSE)
      } else {
        warning(paste("File ", newfilepath,
                " not added because already exists in given directory."))
        newfilepath <- NA
      }
    }

    newJoint$file[i2] <- newfilepath
  }

  newJoint <- newJoint[!is.na(newJoint$file), ]

  obj@samples <- newJoint

  if ("QC" %in% obj@samples$group) QC(obj) <- "QC"

  return(obj)
}




# TODO Adapt add metadata to ntsData-class

#' @title addMetadata
#' @description Adds additional information to the sampleInfo \code{data.frame}. For instance,
#' additional information could be concentration data or other known properties or classifiers for each sample in the sampleInfo \code{data.frame}.
#' The added metadata can then be used for inferential and statistical analyses.
#'
#' @param sampleInfo The sampleInfo \code{data.frame} object to be amended with metadata.
#' @param metadata A \code{data.frame} with the same number of rows as the sampleInfo contaninig 1 or more columns with metadata.
#' Note that the column names in the metadata \code{data.frame} will be used as classifiers.
#'
#' @return Returns an updated sampleInfo \code{data.frame} objects.
#'
#' @export
#'
addMetadata <- function(sampleInfo, metadata) {
  if (length(sampleInfo) == 0 & length(sampleInfo) != length(metadata)) {
    warning("Please make sure the metadata has the same dimensions as the sample data")
  } else {
    sampleInfo <- cbind(sampleInfo, metadata)
  }
  return(sampleInfo)
}




#' @title mzMLconverter
#' @description Converts raw MS files from different formats to mzML using the external tool \emph{msConvert}
#' from [ProtesWizard](http://proteowizard.sourceforge.net/) through the interface of \pkg{patRoon} \insertCite{Helmus2021}{ntsIUTA}.
#' Possible formats are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.
#'
#' @param path The \code{path} where the files to be converted can be found. The function will look into subfolders.
#' @param files Alternatevely, a list with complete file paths for conversion.
#' Note that this will overwrite the \code{path} argument to look for files. The \code{convertFrom} argument must be anyway given.
#' @param convertFrom The format/vendor of files to be found and converted.
#' Possible entries are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.
#' The default is \emph{agilent}.
#' @param centroidMethod The method for centroiding the data. Possible values are \emph{vendor} (the default) to use the vendor algorithm,
#' \emph{cwt} to use the wavelet algorithm or \code{NULL} for skiping dada centroiding. See [ProtesWizard](http://proteowizard.sourceforge.net/) for more information.
#' @param outPath The path to store the mzML files. When \code{NULL} the \code{path} is used instead.
#' @param overWrite Logical, set to \code{TRUE} to overwrite existing files in the defined storing location (either \code{path} or \code{outPath}).
#'
#' @return The converted files will be saved as mzML in the given \code{path} or \code{outPath}.
#'
#' @references \insertRef{Helmus2021}{ntsIUTA}
#' @references \insertRef{proteo}{ntsIUTA}
#'
#' @export
#'
#' @importFrom patRoon convertMSFiles
#'
mzMLconverter <- function(path = getwd(),
                          files = NULL,
                          convertFrom = "agilent",
                          centroidMethod = "vendor",
                          outPath = NULL,
                          overWrite = FALSE) {

  msFileFilters <- data.frame(type = c("agilent", "thermo", "ab", "waters"),
                              ext =  c(".\\.d$", ".\\.RAW$", ".\\.wiff$", ".\\.RAW$"))

  if (convertFrom %in% msFileFilters$type) {

    if (is.null(files)) {
      fileType <- msFileFilters[msFileFilters$type == convertFrom, "ext", drop = TRUE]
      if (length(fileType) == 1) {
        rawFiles <- list.files(path = path, pattern = fileType,
                                    recursive = TRUE, full.names = TRUE,
                                    no.. = FALSE, include.dirs = TRUE, ignore.case = TRUE)
      } else {
        return(cat("Warning: files with the specified format not found in the given path.
        See ?ntsIUTA::mzMLconverter for information."))
      }

    } else {
      matchfiles <- grepl(pattern = msFileFilters[msFileFilters$type == convertFrom, "ext", drop = TRUE], files)
      if (unique(matchfiles)) {
        rawFiles <- files
      } else {
        return(cat("Warning: one or more files do not match with the specified format.
        See ?ntsIUTA::mzMLconverter for information."))
      }
    }

    if (length(rawFiles) > 1) {
      if (!is.null(centroidMethod) & centroidMethod %in% c("vendor", "cwt")) {
        patRoon::convertMSFiles(files = rawFiles,
                                outPath = ifelse(is.null(outPath), path, outPath),
                                dirs = FALSE,
                                anaInfo = NULL,
                                from = convertFrom,
                                to = "mzML",
                                overWrite = overWrite,
                                algorithm = "pwiz",
                                centroid = centroidMethod,
                                filters = c("msLevel 1-2"),
                                extraOpts = NULL,
                                PWizBatchSize = 1)
        return(cat("done"))

      } else {
        patRoon::convertMSFiles(files = rawFiles,
                                outPath = ifelse(is.null(outPath), path, outPath),
                                dirs = FALSE,
                                anaInfo = NULL,
                                from = convertFrom,
                                to = "mzML",
                                overWrite = overWrite,
                                algorithm = "pwiz",
                                centroid = NULL,
                                filters = NULL,
                                extraOpts = NULL,
                                PWizBatchSize = 1)
        return(cat("done"))
      }
    } else {
      return(cat("No files found with the given format or vendor.
      See ?ntsIUTA::mzMLconverter for information."))
    }

  } else {
    return(cat("Warning: the format given in convertFrom is not recognized.
    See ?ntsIUTA::mzMLconverter for possible data formats."))
  }
}
