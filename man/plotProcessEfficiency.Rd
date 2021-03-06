% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fct_plotMonitoring.R
\name{plotProcessEfficiency}
\alias{plotProcessEfficiency}
\title{plotProcessEfficiency}
\usage{
plotProcessEfficiency(obj, sequences = NULL)
}
\arguments{
\item{obj}{An \linkS4class{ntsMonitoringData} object with categorized features
to evaluate the process efficiency.}

\item{sequences}{Optional vector with sequence names from the \linkS4class{ntsMonitoringData} object to plot.
When \code{NULL}, all sequences are plotted.}
}
\value{
A bar plot with the number of features
for each category for each step in the respective process sequence.
}
\description{
plotProcessEfficiency
}
