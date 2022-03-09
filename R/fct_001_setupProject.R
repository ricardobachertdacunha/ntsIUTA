

#' @title setupProject
#'
#' @description \link{setupProject} takes a given path to create a project for Non-target screening (NTS).
#' mzML and/or mzXML files in the path are directly added to the project.
#' Yet, the function \code{\link{addFiles}} can be used to add files or extra files to the project.
#' Setting \code{createRproject} to \code{TRUE} will create and open an R project in the selected or given folder.
#' When \code{convertFiles} is \code{TRUE} the function uses the function \code{\link{mzMLconverter}} to automatically convert
#' specified MS files in the path to mzML and add them to the project.
#'
#' @param path The directory for project setup. The default is the \code{utils::choose.dir()} to select or create a folder.
#' Note that the function will look into the folder or subfolders for mzML and/or mzXML files
#' or specified raw files when \code{convertFiles} is set to \code{TRUE}.
#' @param title A character string with the project title.
#' @param description A character string with the project description.
#' @param date The project date.
#' @param convertFiles Logical, set to \code{TRUE} for converting vendor MS files to mzML.
#' The default is \code{FALSE} for no conversion even if vendor MS files are present.
#' @param convertFrom The name of the vendor format to be found and converted.
#' Possible values are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.
#' @param convertToCentroid Logical, set to \code{TRUE} to centroid profile data when converting to mzML. The default is \code{TRUE}.
#' @param replicates Optional character vector with the same length as the number of MS files in the project folder
#' with the identifier/name of the sample replicate group for each file.
#' The default is \code{NULL} which attempts  to group the samples using
#' the name of the MS files. The groups can always be edited using
#' the class method \link{sampleGroups<-}.
#' See \code{?"sampleGroups<-"} for more information.
#' @param polarity A vector specifying the ionization polarity of the files.
#' Possible values are \emph{positive} (the default) or \emph{negative} for positive or negative modes, respectively.
#' @param method A character vector with the method used for acquiring the MS data.
#' When different methods were used, the method should be given per file.
#' Note that the character vector must have the same length as the number of mzML/mzXML files in the project folder.
#' @param save Logical, set to \code{TRUE} to save the \linkS4class{ntsData} object in the \strong{rData} folder.
#' If \strong{rData} folder does not exist in given \code{path}, the folder will be created.
#' @param makeNewProject Logical, set to TRUE to create an R project in the given \code{path} and open a new R session.
#'
#' @return An S4 class object named \linkS4class{ntsData} containing initialisation data for the NTS project.
#'
#' @details The \linkS4class{ntsData} object can be edited using the class methods.
#' See \code{?"ntsData-class"} for more information about the slots and methods.
#'
#' @export
#'
#' @importFrom checkmate testChoice
#' @importFrom methods new
#' @importFrom rstudioapi initializeProject openProject
#'
#' @examples
#' path <-  system.file(package = "ntsIUTA", dir = "extdata")
#' object <- setupProject(path = path, save = FALSE)
#' object
#'
setupProject <- function(path = utils::choose.dir(getwd(), "Select or create a project folder"),
                         title = NA_character_,
                         description = NA_character_,
                         date = Sys.Date(),
                         convertFiles = FALSE,
                         convertFrom = NULL,
                         convertToCentroid = TRUE,
                         replicates = NULL,
                         polarity = "positive",
                         method = NA_character_,
                         save = TRUE,
                         makeNewProject = FALSE) {

  object <- new("ntsData")
  object@title <- as.character(title)
  object@description <- as.character(description)
  object@path <- as.character(path)
  object@date <- date

  if (convertFiles) {
    if (testChoice(convertFrom, c("thermo", "bruker", "agilent", "ab", "waters"))) {
      mzMLconverter(
        path = path,
        convertFrom = convertFrom,
        centroidMethod = ifelse(convertToCentroid == TRUE, "vendor", NULL)
      )
    } else {
      stop("Vendors should be specified for file recognition.
            Possible entries are: thermo, bruker, agilent, ab (from AB Sciex) and waters.")
    }
  }

  files <- list.files(
    path = object@path,
    pattern = ".mzML|.mzXML",
    recursive = TRUE,
    full.names = TRUE,
    no.. = TRUE
  )

  if (length(files) == 0) {
    warning("No mzML or mzXML files were found in selected project path.
            You can use the addFiles funtion to manually add files
            to the project (see ?addFiles for more information).")
  } else {
    object <- addFiles(
      files = files,
      object = object,
      copy = FALSE,
      replicates = replicates,
      polarity = polarity,
      method = method
    )
  }

  object@metadata <- object@samples[, .(sample, replicate)]

  if (makeNewProject) {

    if (save) saveObject(object = object)

    #TODO add template for main script file
    sp <- file.path(object@path, "script.R")
    cat(
"
\n
#Run the following code to load the project setup ntsData object:\n
library(ntsIUTA)
\n
setup <- readRDS('rData/ntsData.rds')",

    file = sp, sep = "")

    if (!(is.na(Sys.getenv()["RSTUDIO"]))) {
      rstudioapi::initializeProject(object@path)
      rstudioapi::openProject(object@path, newSession = TRUE)
    } else {
      setwd(object@path)
      print("When not using RStudio, it is recommended
          to open the project folder through the IDE
          so that the files and workspace can be loaded. \n
          Run  rstudioapi::navigateToFile('mainscript.R')
          in the new project to open the script.R file.")
      return(object)
    }
  } else {
    if (save) saveObject(object = object)
    return(object)
  }
}




#' @title addFiles
#'
#' @description Adds mzML and/or mzXML files to an \linkS4class{ntsData} object.
#'
#' @param files A list of paths for selected mzML and/or mzXML files.
#' The default is a UI to choose the files.
#' @param object An \linkS4class{ntsData} object to add the selected files.
#' @param copy Logical. Set to \code{TRUE} to copied the selected files
#' to the \code{path} of the \linkS4class{ntsData} object.
#' \code{FALSE} will add the original location of the files
#' to the \code{samples} table.
#' Note that if duplicate names exist, the process is aborted.
#' @param replicates Optional character vector with the same length as the number of MS files to add
#' with the identifier/name of the sample replicate group for each file.
#' The default is \code{NULL} which attempts  to group the replicate samples using
#' the characters before a last \emph{-} (hyphen) or \emph{_} (underscore) of the MS file names.
#' The replicate names can always be edited using the class method \link{sampleGroups<-}.
#' See \code{?"sampleGroups<-"} for more information.
#' @param polarity A vector specifying the ionization polarity of the files.
#' Possible values are \emph{positive} (the default) or \emph{negative} for positive or negative modes, respectively.
#' Note that it will overwrite the polarity of the \linkS4class{ntsData} object.
#' @param method A character vector with the method used for acquiring the MS data.
#' When different methods were used, the method should be given per file.
#' So the character vector must have the same length as the number of MS files to add.
#'
#' @return Returns the \linkS4class{ntsData} object including the added files
#' in the \code{samples} slot.
#'
#' @note When the \linkS4class{ntsData} has files, new files are added
#' below in the \code{samples}. If processed data is already present,
#' the new files will invalidate the structure.
#' Therefore, data should be processed again.
#'
#' @importFrom checkmate assertClass testChoice
#' @importFrom data.table data.table
#' @importFrom stringr str_replace
#' @importFrom tools file_ext
#'
#' @export
#'
addFiles <- function(files = utils::choose.files(),
                     object = ntsData,
                     copy = FALSE,
                     replicates = NULL,
                     polarity = "positive",
                     method = NULL) {

  assertClass(object, "ntsData")
  samples <- object@samples
  files <- sort(files)

  if (length(files) < 1) return(object)

  samples2 <- data.table(
    file = character(),
    sample = character(),
    replicate = character(),
    blank = character(),
    method = character(),
    scans = numeric(),
    centroided = logical(),
    msLevels = character(),
    rtStart = numeric(),
    rtEnd = numeric(),
    mzLow = numeric(),
    mzHigh = numeric(),
    CE = character()
  )

  samples2 <- samples2[seq_len(length(files)), ]
  samples2$file <- files
  samples2$sample <- gsub(".mzML|.mzXML", "", basename(files))

  if (is.null(replicates)) {
    samples2$replicate <- samples2$sample
    samples2$replicate <- str_replace(samples2$replicate, "-", "_")
    samples2$replicate <- sub("_[^_]+$", "", samples2$replicate)
  } else {
    if (nrow(samples2) == length(replicates)) samples2$replicates <- replicates
  }

  if (testChoice(polarity, c("positive", "negative"))) object@polarity <- polarity

  if (length(method) == 1 | length(method) == nrow(samples2)) {
    samples2$method <- method
  }

  samples2 <- rbind(samples, samples2)
  duplicates <- duplicated(samples2$sample)

  if (TRUE %in% duplicates) {
    warning("Files not added because duplicate file names found! Rename files with the same name.")
    return(object)
  }

  if (copy) {
    msdir <- paste0(object@path, "/msfiles")
    if (!dir.exists(msdir)) dir.create(msdir)

    for (i in seq_len(length(files))) {
      i2 <- which(samples2$file == files[i])
      ext <- file_ext(files[i])

      newfilepath <- paste(msdir, "/", samples2$sample[i2], ".", ext, sep = "")

      if (!file.exists(newfilepath)) {
        file.copy(files[i], newfilepath, overwrite = FALSE)
        samples2$file[i2] <- newfilepath
      }
    }
  }

  object@samples <- samples2

  object <- addRawInfo(object)

  return(object)
}




#' @title addMetadata
#'
#' @description Adds metadata for each sample replicate in an \linkS4class{ntsData} object.
#' If metadata already exists in the \linkS4class{ntsData} object,
#' the new variable/s are added as new column/s.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param var A \code{data.frame}, \code{data.table} or a \code{vector} with the same
#' row number or length as the number of sample replicates in the \linkS4class{ntsData}.
#' For \code{data.frame} or \code{data.table}, multiple columns can be added.
#' @param varname A character vector to name the column in the metadata
#' \code{data.table} when metadata is a \code{vector}. The default is \emph{var}.
#'
#' @return An \linkS4class{ntsData} with metadata added to the slot metadata.
#'
#' @importFrom checkmate assertClass
#' @importFrom data.table as.data.table data.table is.data.table
#'
#' @export
#'
addMetadata <- function(object, var = NULL, varname = "var") {

  assertClass(object, "ntsData")

  if (is.data.frame(var) | is.data.table(var)) {
    if (nrow(var) == nrow(object@metadata)) {
      var <- as.data.table(var)
      object@metadata <- cbind(object@metadata, var)
    } else {
      warning("var does not have the same number of rows as sample replicates, metadata not added.")
      return(object)
    }
  }

  if (is.vector(var)) {
    if (length(var) == nrow(object@metadata)) {
      var <- data.table(var = var)
      colnames(var) <- varname
      object@metadata <- cbind(object@metadata, var)
    } else {
      warning("var does not have the same length as sample replicates, metadata not added.")
      return(object)
    }
  }
  return(object)
}




#' @title removeMetadata
#'
#' @description Removes metadata in an \linkS4class{ntsData} object.
#'
#' @param object An \linkS4class{ntsData} object.
#' @param varname A character vector with column metadata name/s
#' to be removed from the object. If \code{NULL}, the default,
#' all metadata is removed.
#'
#' @return An \linkS4class{ntsData} with specified metadata removed for the slot metadata.
#'
#' @importFrom checkmate assertClass
#'
#' @export
#'
removeMetadata <- function(object, varname = NULL) {

  assertClass(object, "ntsData")

  if (is.null(varname)) {
    object@metadata <- object@metadata[, .(sample, replicate)]
    return(object)
  }

  if (FALSE %in% (varname %in% colnames(metadata(object)))) {
    warning("Given varname/s not found in the metadata of the object!")
    return(object)
  }

  varname <- varname[!varname %in% c("sample", "replicate")]
  object@metadata <- object@metadata[, -varname, with = FALSE]
  return(object)
}




#' @title mzMLconverter
#'
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
#' @param overwrite Logical, set to \code{TRUE} to overwrite existing files in the defined storing location (either \code{path} or \code{outPath}).
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
                          overwrite = FALSE) {

  msFileFilters <- data.frame(type = c("agilent", "thermo", "ab", "waters"),
                              ext =  c(".\\.d$", ".\\.RAW$", ".\\.wiff$", ".\\.RAW$"))

  if (convertFrom %in% msFileFilters$type) {

    if (is.null(files)) {
      fileType <- msFileFilters[msFileFilters$type == convertFrom, "ext", drop = TRUE]
      if (length(fileType) == 1) {
        rawFiles <- list.files(
          path = path,
          pattern = fileType,
          recursive = TRUE,
          full.names = TRUE,
          no.. = FALSE,
          include.dirs = TRUE,
          ignore.case = TRUE
        )
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

    if (length(rawFiles) >= 1) {
      if (!is.null(centroidMethod) & centroidMethod %in% c("vendor", "cwt")) {
        patRoon::convertMSFiles(files = rawFiles,
                                outPath = ifelse(is.null(outPath), path, outPath),
                                dirs = FALSE,
                                anaInfo = NULL,
                                from = convertFrom,
                                to = "mzML",
                                overWrite = overwrite,
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
                                overWrite = overwrite,
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
