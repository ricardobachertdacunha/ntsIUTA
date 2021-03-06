% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_003_getPatData.R
\name{getPatData}
\alias{getPatData}
\title{getPatData}
\usage{
getPatData(featData, sampleInfo, save = FALSE, projPath = setup$projPath)
}
\arguments{
\item{featData}{A \linkS4class{XCMSnExp} object with feature data to be converted.}

\item{sampleInfo}{The \code{sampleInfo} obtained by the \code{\link{setupProject}} function, matching the \code{featData}.}

\item{save}{Logical, set to \code{TRUE} to save the generated \linkS4class{featureGroups} object in the disk.}

\item{projPath}{The \code{projPath} in the \code{setup} object, or the location where to save the generated R object.}
}
\value{
A \linkS4class{featureGroups} object corresponding matching the given \linkS4class{XCMSnExp} object.
}
\description{
getPatData
}
