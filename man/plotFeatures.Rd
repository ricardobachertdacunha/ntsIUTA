% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_plotFeatures.R
\name{plotFeatures}
\alias{plotFeatures}
\title{plotFeatures}
\usage{
plotFeatures(
  obj,
  samples = NULL,
  ID = NULL,
  mz = NULL,
  ppm = 20,
  rt = NULL,
  rtWindow = NULL,
  rtUnit = "sec",
  msLevel = 1,
  colorBy = "features",
  interactive = TRUE
)
}
\arguments{
\item{obj}{An \linkS4class{ntsData} object.}

\item{samples}{The index or name of the sample/s.
The default is \code{NULL} and all samples are used.}

\item{ID}{The identifier of the features of interest.
When not \code{NULL}, overwrites any given \code{mz} and \code{rt} value.}

\item{mz}{Optional target \emph{m/z} to find features using
the mass deviation specified by the \code{ppm} parameter.}

\item{ppm}{The mass deviation to extract the features
when \code{mz} is specified.}

\item{rt}{The retention time in minutes or seconds,
depending on the defined \code{rtUnit}, see below.
Only used when \code{mz} is specified.}

\item{rtWindow}{The time deviation to collect features.
The time unit is the defined by \code{rtUnit}.
A time interval can be given with a length 2 vector,
defining the minimum and maximum retention time.}

\item{rtUnit}{Possible entries are \code{min} or \code{sec}.
The default is \code{min}.}

\item{msLevel}{The MS level to extract the data.
For the moment, only 1 is possible.}

\item{colorBy}{Possible values are \code{"features"},
\code{"samples"} or \code{sampleGroups} (the default),
for colouring by samples or sample replicate groups respectively.}

\item{interactive}{Logical, set to \code{TRUE} to use
the \pkg{plotly} instead of \pkg{ggplot2}. The default is \code{TRUE}.}
}
\value{
A feature/s plot produced through \pkg{gglpot2} for interactive
\code{FALSE} or an interactive plot through the \pkg{plotly}.
}
\description{
Method for plotting features from a \linkS4class{ntsData} object.
}
