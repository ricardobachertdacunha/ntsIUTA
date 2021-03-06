% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_001_makeSetup.R
\name{setupProject}
\alias{setupProject}
\title{setupProject}
\usage{
setupProject(
  path = utils::choose.dir(getwd(), "Select or create a project folder"),
  title = NA_character_,
  description = NA_character_,
  date = Sys.Date(),
  polarity = "positive",
  convertFiles = FALSE,
  convertFrom = NULL,
  convertToCentroid = TRUE,
  groups = NULL,
  save = TRUE,
  makeNewProject = FALSE
)
}
\arguments{
\item{path}{The directory for project setup. The default is the \code{utils::choose.dir()} to select or create a folder.
Note that the function will look into subfolders for mzML/mzXML files or specified raw files when \code{convertFiles} is set to \code{TRUE}.}

\item{title}{A character string with the project title.}

\item{description}{A character string with the project description.}

\item{date}{The project date.}

\item{polarity}{A vector specifying the polarity mode.
Possible values are \code{positive} or \code{negative} for positive and negative modes, respectively.}

\item{convertFiles}{Logical, set to \code{TRUE} for non mzML/mzXML files being converted to centroided mzML.
The default is \code{FALSE} for no conversion even if files are present.}

\item{convertFrom}{The name of the file format or vendor format to be found and converted.
Possible values are: \emph{thermo}, \emph{bruker}, \emph{agilent}, \emph{ab} (from AB Sciex) and \emph{waters}.}

\item{convertToCentroid}{Logical, set to \code{TRUE} to convert profile data to centroided when converting raw MS files to mzML.
The default and recommended is \code{TRUE}.}

\item{groups}{A vector with the identifier/name of each sample replicate group. If \code{NULL},
a grouping tentative will be made using the name of the MS files found in the given \code{path}.}

\item{save}{Logical, set to \code{TRUE} to save the \linkS4class{ntsData} object in the \strong{rData} folder.
If \strong{rData} folder does not exist in given \code{path}, the folder will be created.}

\item{makeNewProject}{Logical, set to TRUE to create an R project in the given \code{path} and open a new R session.}
}
\value{
An S4 class object named \linkS4class{ntsData}
containing initialisation data for the NTS project.
}
\description{
\code{setupProject} takes a given path to create a project for Non-target screening (NTS).
mzML or mzXML files in the path are directly added to the project.
Yet, the function \code{\link{addFiles}} can be used to add files or extra files to the project.
Setting \code{createRproject} to \code{TRUE} will create and open an R project in the selected or given folder.
When \code{convertFiles} is \code{TRUE} the function uses the function \code{\link{mzMLconverter}} to automatically convert
specified raw MS files in the path to mzML and add them to the project.
}
\details{
The \linkS4class{ntsData} object can be edited using the class methods.
See \linkS4class{ntsData} for more information about the slots and methods.
}
\examples{
path <-  system.file(package = "ntsIUTA", dir = "extdata")
dt <- setupProject(path = path, save = FALSE)
dt

}
