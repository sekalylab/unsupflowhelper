% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/controls.R
\name{get_reference_matrix}
\alias{get_reference_matrix}
\title{Get reference matrix}
\usage{
get_reference_matrix(flow_object, add_sample_id = TRUE)
}
\arguments{
\item{flow_object}{A Flow Object}

\item{add_sample_id}{Add a SampleID column to the matrix. Defaults to TRUE.}
}
\value{
A data.frame or matrix.
}
\description{
Extract fluorescence values from the reference control FlowSet.
}
