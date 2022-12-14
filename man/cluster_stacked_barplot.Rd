% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{cluster_stacked_barplot}
\alias{cluster_stacked_barplot}
\title{Plot Stacked barplot}
\usage{
cluster_stacked_barplot(flow_object, annotation = "louvain", group = NULL)
}
\arguments{
\item{flow_object}{A Flow Object}

\item{annotation}{Cluster annotation to display. Defaults to `louvain`.}

\item{group}{Pass a metadata variable to compare frequencies across. Defaults to `NULL`.}
}
\value{
A ggplot object
}
\description{
Graphs a histogram, similar to a 1D FlowJo plot, of the expression of a selected marker for each cluster in a Flow Object.
}
