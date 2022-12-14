% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudotime.R
\name{plot_grouped_pseudotime}
\alias{plot_grouped_pseudotime}
\title{Plot a pseudotime density/relative frequency plot}
\usage{
plot_grouped_pseudotime(
  flow_object,
  trajectory,
  group.by,
  adjust = 1,
  y = "freq"
)
}
\arguments{
\item{flow_object}{A Flow Object}

\item{trajectory}{Name of the trajectory to plot.}

\item{group.by}{Sample or cluster annotation to separate cells by.}

\item{adjust}{Bin size adjustment to compute normalization over. Defaults to 1.}

\item{y}{Parameter to be displayed on the Y-axis. Accepted values are 'freq' or 'density'. Defaults to 'freq'.}
}
\value{
A ggplot2 object
}
\description{
Display a density plot for pseudotime on the X-axis, and relative density on the Y-axis. Can be grouped by metadata annotations.
}
