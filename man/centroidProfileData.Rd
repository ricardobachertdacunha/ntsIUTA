% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_002_importRawData.R
\name{centroidProfileData}
\alias{centroidProfileData}
\title{centroidProfileData}
\usage{
centroidProfileData(
  obj,
  halfwindow = 2,
  SNR = 0,
  noiseMethod = "MAD",
  smoothing = FALSE,
  methodSmoothing = "SavitzkyGolay",
  methodRefineMz = "kNeighbors",
  k = 1,
  signalPercentage = 10,
  stopAtTwo = TRUE,
  save = FALSE,
  ...
)
}
\arguments{
\item{obj}{A \linkS4class{ntsData} object with profile data for centroiding.}

\item{halfwindow}{Sets the window size for centroiding as \code{2 * halfwindow + 1}.
The \code{halfwindow} should be slightly larger than the full width
at half maximum of the profile peak.}

\item{SNR}{The signal-to-noise ratio to consider a local maximum as peak.}

\item{noiseMethod}{The method to estimate the noise level.
Possible methods are "MAD" (the default) and "SuperSmoother".
See \code{?MSnbase::pickPeaks} for more information.}

\item{smoothing}{Logical, set to \code{TRUE} for applying smothing
to the profile data before centroiding. The default is FALSE.}

\item{methodSmoothing}{Method for data smoothing.
The possible methods are "SavitzkyGolay" (the default) and "MovingAverage".
See \code{?MSnbase::smooth} for more information and arguments,
which are passed by \code{...}.}

\item{methodRefineMz}{Method for refinement.
Possible methods are "none" (the default, for not applying \emph{m/z} refinement),
"kNeighbors" and "descendPeak". See \code{?MSnbase::pickPeaks} for more information.}

\item{k}{When refine method is "kNeighbors",
\code{k} is number of closest signals to the centroid.}

\item{signalPercentage}{When refine method is "descendPeak",
\code{signalPercentage} is the minimum signal percentage of centroids to refine \emph{m/z}.}

\item{stopAtTwo}{Logical, when refine method is "descendPeak",
set to \code{TRUE} for allowing two consecutive equal or higher signals.
\code{FALSE} will stop when one equal or higher centroid is found.}

\item{save}{Logical, set to \code{TRUE} to replace
the original files by the centroided files in disk.
The location is taken from the originbal file paths.}

\item{...}{Arguments for selected smoothing method.
See \code{?MSnbase::smooth} for possible arguments for each method.}
}
\value{
Centroided \linkS4class{ntsData} object.
When \code{save} is set to TRUE, the profile data in the original
mzML or mzXML files is replaced by the centroided data.
}
\description{
Centroiding of profile data with additional possibility
for data smoothing before centroiding and \emph{m/z} refinement.
The \code{centroidProfileData} function combines functions \code{smooth}
and \code{pickPeaks} from the \code{MSnbase} package, see references.
}
\references{
\insertRef{MSnbase2}{ntsIUTA}
}
