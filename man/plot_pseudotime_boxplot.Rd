% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pseudotime.R
\name{plot_pseudotime_boxplot}
\alias{plot_pseudotime_boxplot}
\title{Plot a pseudotime boxplot}
\usage{
plot_pseudotime_boxplot(
  flow_object,
  trajectory,
  metric = "mean",
  group.by,
  stats = NULL
)
}
\arguments{
\item{flow_object}{A Flow Object}

\item{trajectory}{Name of the trajectory to plot.}

\item{metric}{Summary statistic to display. Accepts `mean` or `median`.}

\item{group.by}{Metadata variable to compare pseudotime.}

\item{stats}{Selection of a statistical test to compare groups (if using the `groups` parameter). Accepts `wilcox`, `t.test` or `none`. Defaults to `wilcox`.}
}
\value{
A ggplot boxplot
}
\description{
Display mean or median pseudotime values within a trajectory across metadata groups
}
