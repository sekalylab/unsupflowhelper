% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics.R
\name{get_cluster_frequencies}
\alias{get_cluster_frequencies}
\title{Get cluster frequencies}
\usage{
get_cluster_frequencies(
  flow_object,
  louvain_select = NULL,
  annotation = "louvain",
  format = "long",
  output = "relative"
)
}
\arguments{
\item{flow_object}{A flow_object}

\item{louvain_select}{A filter for a subset of louvain clusters for which you want frequencies. NULL by default.}

\item{annotation}{Add additional annotations to the data, including clusters, cluster annotations and sample metadata, and dimensions reductions. Set to 'louvain' by default.}

\item{format}{The format of the output. Can be set to 'long' or 'wide'. Annotations can only be added in a long format. Defaults to 'long'.}

\item{output}{Retrieve relative or absolute frequencies. If the downsampling was NOT equalized at flow object creation, relative frequencies are recommended.}
}
\value{
A data.frame
}
\description{
Extract cluster frequencies from a flow object. If the object contains cluster annotations, the function will calculate frequencies for the indicated annotation.
}
